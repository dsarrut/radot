from pathlib import Path
from typing import Optional
import SimpleITK as sitk

def convert_image_units(
    image: sitk.Image,
    to_unit: str,
    from_unit: Optional[str] = None,
    weight: Optional[float] = None,
    activity: Optional[float] = None,
    image_path: Optional[Path] = None,
) -> sitk.Image:
    """
    Convert image voxel values between units (bq, bqml, suv).
    
    Parameters:
    - image: sitk.Image, input image.
    - to_unit: str, target unit ("bq", "bqml", "suv").
    - from_unit: Optional[str], source unit ("bq", "bqml", "suv"). If None, read from metadata.
    - weight: Optional[float], patient weight in kg (required for SUV conversions). If None, read from metadata.
    - activity: Optional[float], injected activity in MBq (required for SUV conversions). If None, read from metadata.
    - image_path: Optional[Path], optional file path of the image for sidecar resolution.
    """
    from radot.meta import get_metadata_case_insensitive

    to_unit = to_unit.lower()

    # 1. Resolve from_unit
    if from_unit is None:
        from_unit = get_metadata_case_insensitive(image, "unit", image_path)
        if from_unit is None:
            raise ValueError("Source unit is not specified and was not found in the image metadata.")
    else:
        from_unit = from_unit.lower()
    
    if from_unit == to_unit:
        return image

    # 2. Resolve weight and activity if SUV is involved
    if "suv" in (from_unit, to_unit):
        if weight is None:
            w_str = get_metadata_case_insensitive(image, "weight", image_path)
            if w_str is not None:
                try:
                    weight = float(w_str)
                except ValueError:
                    raise ValueError(f"Invalid weight format in image metadata: {w_str}")
            else:
                raise ValueError("A positive --weight (in kg) is required for SUV conversion (not found in image metadata).")
        if weight <= 0:
            raise ValueError("Weight must be positive.")

        if activity is None:
            act_str = get_metadata_case_insensitive(image, "activity", image_path)
            if act_str is not None:
                try:
                    activity = float(act_str)
                except ValueError:
                    raise ValueError(f"Invalid activity format in image metadata: {act_str}")
            else:
                raise ValueError("A positive --activity (in MBq) is required for SUV conversion (not found in image metadata).")
        if activity <= 0:
            raise ValueError("Activity must be positive.")

    # Calculate voxel volume in mL (1 mL = 1000 mm^3)
    spacing = image.GetSpacing()
    voxel_volume_ml = (spacing[0] * spacing[1] * spacing[2]) / 1000.0

    # Step 1: Convert input to bqml (Becquerels per mL)
    # We cast to Float32 to ensure floating point division works properly
    converted = sitk.Cast(image, sitk.sitkFloat32)
    
    if from_unit == "bq":
        # bq -> bqml
        converted = converted / voxel_volume_ml
    elif from_unit == "suv":
        # suv -> bqml
        # SUV = bqml / (activity_bq / weight_g)
        # bqml = SUV * (activity_bq / weight_g)
        activity_bq = activity * 1_000_000.0
        weight_g = weight * 1000.0
        suv_factor = activity_bq / weight_g
        converted = converted * suv_factor
    elif from_unit != "bqml":
        raise ValueError(f"Unsupported source unit: {from_unit}")

    # Step 2: Convert bqml to target unit
    if to_unit == "bq":
        # bqml -> bq
        converted = converted * voxel_volume_ml
    elif to_unit == "suv":
        # bqml -> suv
        activity_bq = activity * 1_000_000.0
        weight_g = weight * 1000.0
        suv_factor = activity_bq / weight_g
        converted = converted / suv_factor
    elif to_unit != "bqml":
        raise ValueError(f"Unsupported target unit: {to_unit}")

    return sitk.Cast(converted, sitk.sitkFloat32)


def convert_image_units_file(
    image_path: Path,
    to_unit: str,
    from_unit: Optional[str] = None,
    output_path: Optional[Path] = None,
    weight: Optional[float] = None,
    activity: Optional[float] = None,
) -> None:
    """
    Load image from path, convert its units, and save to output path.
    """
    if not image_path.exists():
        raise FileNotFoundError(f"Input image not found: {image_path}")

    to_unit = to_unit.lower()

    if output_path is None:
        suffix = "".join(image_path.suffixes)
        base = image_path.name[:-len(suffix)] if suffix else image_path.name
        output_path = image_path.parent / f"{base}_{to_unit}{suffix}"

    print(f"Reading image from: {image_path}")
    image = sitk.ReadImage(str(image_path))

    # Resolve from_unit if None for print info
    resolved_from_unit = from_unit
    if resolved_from_unit is None:
        from radot.meta import get_metadata_case_insensitive
        resolved_from_unit = get_metadata_case_insensitive(image, "unit", image_path)
        if resolved_from_unit is None:
            raise ValueError("Source unit is not specified and was not found in the image metadata.")
    
    print(f"Converting units from {resolved_from_unit} to {to_unit}...")
    converted_image = convert_image_units(
        image=image,
        from_unit=from_unit,
        to_unit=to_unit,
        weight=weight,
        activity=activity,
        image_path=image_path
    )

    print(f"Saving converted image to: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(converted_image, str(output_path))

    # Write updated metadata to output sidecar JSON
    from radot.meta import get_sidecar_path, get_metadata_case_insensitive
    import json
    
    modality = get_metadata_case_insensitive(image, "modality", image_path)
    weight_val = get_metadata_case_insensitive(image, "weight", image_path)
    activity_val = get_metadata_case_insensitive(image, "activity", image_path)
    injection_date = get_metadata_case_insensitive(image, "injection_date", image_path)
    acquisition_date = get_metadata_case_insensitive(image, "acquisition_date", image_path)

    new_meta = {}
    if modality:
        new_meta["modality"] = modality
    new_meta["unit"] = to_unit
    if weight_val:
        new_meta["weight"] = float(weight_val)
    if activity_val:
        new_meta["activity"] = float(activity_val)
    if injection_date:
        new_meta["injection_date"] = injection_date
    if acquisition_date:
        new_meta["acquisition_date"] = acquisition_date

    out_sidecar = get_sidecar_path(output_path)
    with open(out_sidecar, "w") as f:
        json.dump(new_meta, f, indent=4)

    print("Unit conversion completed successfully.")
