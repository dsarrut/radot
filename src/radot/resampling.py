from pathlib import Path
from typing import Tuple, Optional
import SimpleITK as sitk

import math


def resample_image(
    image: sitk.Image,
    spacing: Optional[Tuple[float, float, float]] = None,
    image_type: str = "ct",
    preserve_total: bool = False,
    gaussian: bool = True,
    reference_image: Optional[sitk.Image] = None,
) -> sitk.Image:
    """
    Resample a SimpleITK image to a target spacing or to match a reference image grid.
    
    Parameters:
    - image: sitk.Image, input image.
    - spacing: target spacing (spacing_x, spacing_y, spacing_z). Optional if reference_image is provided.
    - image_type: "ct", "labels", "pet", or "spect".
    - preserve_total: if True, scales intensities by (new_voxel_volume / old_voxel_volume) to preserve the sum of voxel values.
    - gaussian: if True, applies a Gaussian anti-aliasing filter along the downsampled axes.
    - reference_image: sitk.Image, optional reference image to align size, spacing, origin, and direction.
    """
    if image_type != "labels":
        image = sitk.Cast(image, sitk.sitkFloat32)

    original_spacing = image.GetSpacing()
    original_size = image.GetSize()
    
    if reference_image is not None:
        target_spacing = reference_image.GetSpacing()
        new_size = reference_image.GetSize()
        target_origin = reference_image.GetOrigin()
        target_direction = reference_image.GetDirection()
    else:
        if spacing is None:
            raise ValueError("Either spacing or reference_image must be specified.")
        target_spacing = spacing
        new_size = [
            max(1, int(round(original_size[i] * original_spacing[i] / target_spacing[i])))
            for i in range(3)
        ]
        target_origin = image.GetOrigin()
        target_direction = image.GetDirection()
    
    # If gaussian smoothing is requested and this is not a label image, apply it before resampling
    if gaussian and image_type != "labels":
        smoothed_image = sitk.Cast(image, sitk.sitkFloat32)
        smoothing_needed = False
        for i in range(3):
            if target_spacing[i] > original_spacing[i]:
                # Calculate sigma in physical units (FWHM mapped to spacing change)
                sigma = math.sqrt(target_spacing[i] ** 2 - original_spacing[i] ** 2) / 2.355
                if sigma > 0.0:
                    gaussian_filter = sitk.RecursiveGaussianImageFilter()
                    gaussian_filter.SetDirection(i)
                    gaussian_filter.SetSigma(sigma)
                    smoothed_image = gaussian_filter.Execute(smoothed_image)
                    smoothing_needed = True
        if smoothing_needed:
            image = smoothed_image

    resampler = sitk.ResampleImageFilter()
    resampler.SetSize(new_size)
    resampler.SetOutputSpacing(target_spacing)
    resampler.SetOutputOrigin(target_origin)
    resampler.SetOutputDirection(target_direction)
    
    # Choose interpolator
    if image_type == "labels":
        resampler.SetInterpolator(sitk.sitkNearestNeighbor)
    else:
        resampler.SetInterpolator(sitk.sitkLinear)


    # Set default pixel value (background value)
    if image_type == "ct":
        resampler.SetDefaultPixelValue(-1000.0)
    else:
        resampler.SetDefaultPixelValue(0.0)

    resampled_image = resampler.Execute(image)

    # Scale intensities if preserve_total is True (for absolute activity Bq)
    if preserve_total:
        old_volume = original_spacing[0] * original_spacing[1] * original_spacing[2]
        new_volume = target_spacing[0] * target_spacing[1] * target_spacing[2]
        volume_ratio = new_volume / old_volume

        # Cast to float to preserve decimals
        if resampled_image.GetPixelID() not in [sitk.sitkFloat32, sitk.sitkFloat64]:
            resampled_image = sitk.Cast(resampled_image, sitk.sitkFloat32)
        resampled_image = resampled_image * volume_ratio

    return resampled_image


def resample_image_file(
    image_path: Path,
    spacing: Tuple[float, float, float],
    image_type: Optional[str] = None,
    output_path: Optional[Path] = None,
    unit: Optional[str] = None,
    preserve_total: bool = False,
    gaussian: bool = True,
) -> None:
    """
    Load an image from path, resample it, and write it to output path.
    """
    if not image_path.exists():
        raise FileNotFoundError(f"Input image not found: {image_path}")

    if output_path is None:
        raise ValueError("output_path must be specified.")

    print(f"Reading image from: {image_path}")
    image = sitk.ReadImage(str(image_path))

    from radot.meta import get_metadata_case_insensitive

    # 1. Resolve image_type (modality) from metadata
    if image_type is None:
        image_type = get_metadata_case_insensitive(image, "modality", image_path)
        if image_type is None:
            raise ValueError("Image modality type is not specified and was not found in the image metadata.")
    image_type = image_type.lower()

    # 2. Resolve unit from metadata
    if unit is None:
        unit = get_metadata_case_insensitive(image, "unit", image_path)
        if unit is None:
            # Default fallbacks
            if image_type == "ct":
                unit = "hu"
            elif image_type in ("pet", "spect"):
                unit = "suv"
            else:
                unit = "none"
    else:
        unit = unit.lower()

    # Validate unit vs image_type
    if unit == "hu" and image_type != "ct":
        raise ValueError("HU units can only be used with CT images.")
    if unit in ["suv", "bqml", "bq"] and image_type not in ["pet", "spect"]:
        raise ValueError("SUV/Bqml/Bq units can only be used with PET or SPECT images.")
    if unit != "none" and image_type == "labels":
        raise ValueError("Labels cannot have a unit option.")

    # Force preserve_total if unit is absolute activity (bq)
    if unit == "bq":
        preserve_total = True

    d = 0
    if image_type == "ct":
        d = -1000
    print(
        f"Resampling (type: {image_type}, unit: {unit}, preserve_total: {preserve_total}, gaussian: {gaussian}), default pixel {d}"
    )
    resampled = resample_image(
        image=image,
        spacing=spacing,
        image_type=image_type,
        preserve_total=preserve_total,
        gaussian=gaussian,
    )

    print(f"Saving resampled image to: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(resampled, str(output_path))

    # Copy sidecar JSON if it exists
    from radot.meta import get_sidecar_path
    import shutil
    in_sidecar = get_sidecar_path(image_path)
    if in_sidecar.exists():
        out_sidecar = get_sidecar_path(output_path)
        shutil.copy2(in_sidecar, out_sidecar)

    print("Resampling completed successfully.")
