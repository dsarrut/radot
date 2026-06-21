from pathlib import Path
from typing import Tuple, Optional
import SimpleITK as sitk

import math


def resample_image(
    image: sitk.Image,
    spacing: Tuple[float, float, float],
    image_type: str = "ct",
    preserve_total: bool = False,
    gaussian: bool = True,
) -> sitk.Image:
    """
    Resample a SimpleITK image to a target spacing.

    Parameters:
    - image: sitk.Image, input image.
    - spacing: target spacing (spacing_x, spacing_y, spacing_z).
    - image_type: "ct", "labels", "pet", or "spect".
    - preserve_total: if True, scales intensities by (new_voxel_volume / old_voxel_volume) to preserve the sum of voxel values.
    - gaussian: if True, applies a Gaussian anti-aliasing filter along the downsampled axes.
    """
    original_spacing = image.GetSpacing()
    original_size = image.GetSize()

    # If gaussian smoothing is requested and this is not a label image, apply it before resampling
    if gaussian and image_type != "labels":
        smoothed_image = sitk.Cast(image, sitk.sitkFloat32)
        smoothing_needed = False
        for i in range(3):
            if spacing[i] > original_spacing[i]:
                # Calculate sigma in physical units (FWHM mapped to spacing change)
                sigma = math.sqrt(spacing[i] ** 2 - original_spacing[i] ** 2) / 2.355
                if sigma > 0.0:
                    gaussian_filter = sitk.RecursiveGaussianImageFilter()
                    gaussian_filter.SetDirection(i)
                    gaussian_filter.SetSigma(sigma)
                    smoothed_image = gaussian_filter.Execute(smoothed_image)
                    smoothing_needed = True
        if smoothing_needed:
            image = smoothed_image

    # Calculate new size based on physical dimensions and target spacing
    new_size = [
        max(1, int(round(original_size[i] * original_spacing[i] / spacing[i])))
        for i in range(3)
    ]

    resampler = sitk.ResampleImageFilter()
    resampler.SetSize(new_size)
    resampler.SetOutputSpacing(spacing)
    resampler.SetOutputOrigin(image.GetOrigin())
    resampler.SetOutputDirection(image.GetDirection())

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
        new_volume = spacing[0] * spacing[1] * spacing[2]
        volume_ratio = new_volume / old_volume

        # Cast to float to preserve decimals
        if resampled_image.GetPixelID() not in [sitk.sitkFloat32, sitk.sitkFloat64]:
            resampled_image = sitk.Cast(resampled_image, sitk.sitkFloat32)
        resampled_image = resampled_image * volume_ratio

    return resampled_image


def resample_image_file(
    image_path: Path,
    spacing: Tuple[float, float, float],
    image_type: str = "ct",
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
    print("Resampling completed successfully.")
