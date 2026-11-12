from pathlib import Path
from typing import Tuple, Optional
import logging

logger = logging.getLogger("radot")
import SimpleITK as sitk

import math

# Per-voxel quantities (activity Bq, time-integrated activity Bq.s): their total is kept when resampled
PER_VOXEL_UNITS = ("bq", "bqs")

# Image types: interpolation and padding depend on it (labels: nearest neighbour, no smoothing)
IMAGE_TYPES = ("ct", "pet", "spect", "dose", "labels")


def check_image_type(image_type: Optional[str]) -> str:
    """Lower-cased image type; an unknown or missing type is an error (never smoothed by default)."""
    if image_type is None:
        raise ValueError(f"Image type is required: one of {', '.join(IMAGE_TYPES)}.")
    t = image_type.lower()
    if t not in IMAGE_TYPES:
        raise ValueError(
            f"Unknown image type '{image_type}': expected one of {', '.join(IMAGE_TYPES)}. "
            f"A label map or a mask is 'labels' (nearest neighbour, no smoothing)."
        )
    return t


def resolve_type_and_unit(
    image: sitk.Image, image_path: Path, image_type: Optional[str], unit: Optional[str]
) -> Tuple[str, str]:
    """Image type and unit for resampling: the option, else the sidecar. An option that contradicts
    the sidecar is an error (resampling never relabels nor converts)."""
    from radot.meta import get_metadata_case_insensitive

    meta_type = get_metadata_case_insensitive(image, "modality", image_path)
    meta_type = meta_type.lower() if meta_type else None
    if image_type is not None:
        t = image_type.lower()
        if meta_type is not None and meta_type != t:
            raise ValueError(
                f"Image type mismatch: the sidecar says '{meta_type}', the option says '{t}'. Fix the "
                f"sidecar (radot meta IMAGE --modality {t}) or drop the option."
            )
    else:
        if meta_type is None:
            raise ValueError("Image modality type is not specified and was not found in the image metadata.")
        t = meta_type
    t = check_image_type(t)

    meta_unit = get_metadata_case_insensitive(image, "unit", image_path)
    meta_unit = meta_unit.lower() if meta_unit else None
    if unit is not None:
        u = unit.lower()
        if meta_unit is not None and meta_unit != u:
            raise ValueError(
                f"Unit mismatch: the sidecar says '{meta_unit}', the option says '{u}'. Resampling does "
                f"not convert units: use radot convert-units first, or fix the sidecar."
            )
    elif meta_unit is not None:
        u = meta_unit
    elif t == "ct":
        u = "hu"
    elif t in ("pet", "spect"):
        raise ValueError("The unit option is required for PET and SPECT modalities (not found in image metadata).")
    else:
        u = "none"

    if u == "hu" and t != "ct":
        raise ValueError("HU units can only be used with CT images.")
    if u in ["suv", "bqml", "bq"] and t not in ["pet", "spect"]:
        raise ValueError("SUV/Bqml/Bq units can only be used with PET or SPECT images.")
    if u != "none" and t == "labels":
        raise ValueError("Labels cannot have a unit option.")
    return t, u


def _write_output_sidecar(image_path: Path, output_path: Path, image_type: str, unit: str) -> None:
    """Write the sidecar of a resampled image: the source keys, the resolved modality and unit, and the
    SHA-256 of the new file. A label map sidecar ({"1": "Liver", ...}) is copied as is."""
    import json
    from radot.meta import get_sidecar_path, compute_file_sha256

    in_sidecar = get_sidecar_path(image_path)
    meta = {}
    if in_sidecar.exists():
        meta = json.loads(in_sidecar.read_text())
    if image_type != "labels":
        meta["modality"] = image_type
        if unit != "none":
            meta["unit"] = unit
        meta.pop("sha256", None)
        sha = compute_file_sha256(output_path)
        if sha:
            meta["sha256"] = sha
    elif not meta:
        return
    get_sidecar_path(output_path).write_text(json.dumps(meta, indent=4))


def resample_image(
    image: sitk.Image,
    spacing: Optional[Tuple[float, float, float]] = None,
    image_type: Optional[str] = None,
    preserve_total: bool = False,
    gaussian: bool = True,
    reference_image: Optional[sitk.Image] = None,
) -> sitk.Image:
    """
    Resample a SimpleITK image to a target spacing or to match a reference image grid.
    
    Parameters:
    - image: sitk.Image, input image.
    - spacing: target spacing (spacing_x, spacing_y, spacing_z). Optional if reference_image is provided.
    - image_type: required, one of IMAGE_TYPES ("ct", "pet", "spect", "dose", "labels").
    - preserve_total: if True, scales intensities by (new_voxel_volume / old_voxel_volume) to preserve the sum of voxel values.
    - gaussian: if True, applies a Gaussian anti-aliasing filter along the downsampled axes.
    - reference_image: sitk.Image, optional reference image to align size, spacing, origin, and direction.
    """
    image_type = check_image_type(image_type)
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
        # Same physical extent: the first voxel edge stays in place (origin = first voxel centre)
        target_direction = image.GetDirection()
        shift = [(target_spacing[i] - original_spacing[i]) / 2.0 for i in range(3)]
        d = target_direction
        target_origin = tuple(
            image.GetOrigin()[r] + sum(d[3 * r + c] * shift[c] for c in range(3)) for r in range(3)
        )
    
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

    logger.info(f"Reading image from: {image_path}")
    image = sitk.ReadImage(str(image_path))

    image_type, unit = resolve_type_and_unit(image, image_path, image_type, unit)

    # Force preserve_total for a per-voxel unit (bq, bqs)
    if unit in PER_VOXEL_UNITS:
        preserve_total = True

    d = 0
    if image_type == "ct":
        d = -1000
    logger.info(
        f"Resampling (type: {image_type}, unit: {unit}, preserve_total: {preserve_total}, gaussian: {gaussian}), default pixel {d}"
    )
    resampled = resample_image(
        image=image,
        spacing=spacing,
        image_type=image_type,
        preserve_total=preserve_total,
        gaussian=gaussian,
    )

    logger.info(f"Saving resampled image to: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(resampled, str(output_path))

    _write_output_sidecar(image_path, output_path, image_type, unit)

    logger.info("Resampling completed successfully.")


def resample_image_file_like(
    image_path: Path,
    like_image_path: Path,
    output_path: Path,
    type: Optional[str] = None,
    unit: Optional[str] = None,
    preserve_total: bool = False,
    gaussian: bool = True,
) -> None:
    """
    Load an image from image_path, resample it to match the grid of the image at like_image_path,
    and save to output_path.
    """
    image_path = Path(image_path)
    like_image_path = Path(like_image_path)
    output_path = Path(output_path)

    if not image_path.exists():
        raise FileNotFoundError(f"Input image not found: {image_path}")
    if not like_image_path.exists():
        raise FileNotFoundError(f"Reference ('like') image not found: {like_image_path}")

    logger.info(f"Reading image from: {image_path}")
    image = sitk.ReadImage(str(image_path))

    logger.info(f"Reading reference image from: {like_image_path}")
    like_image = sitk.ReadImage(str(like_image_path))

    resolved_type, unit = resolve_type_and_unit(image, image_path, type, unit)

    # Force preserve_total for a per-voxel unit (bq, bqs)
    if unit in PER_VOXEL_UNITS:
        preserve_total = True

    logger.info(
        f"Resampling like-grid (type: {resolved_type}, unit: {unit}, preserve_total: {preserve_total}, gaussian: {gaussian})"
    )

    resampled = resample_image(
        image=image,
        image_type=resolved_type,
        preserve_total=preserve_total,
        gaussian=gaussian,
        reference_image=like_image,
    )

    logger.info(f"Saving resampled image to: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(resampled, str(output_path))

    _write_output_sidecar(image_path, output_path, resolved_type, unit)

    logger.info("Resampling completed successfully.")
