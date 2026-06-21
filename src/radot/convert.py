from pathlib import Path
from typing import Optional
import SimpleITK as sitk

def convert_image_units(
    image: sitk.Image,
    from_unit: str,
    to_unit: str,
    weight: Optional[float] = None,
    activity: Optional[float] = None,
) -> sitk.Image:
    """
    Convert image voxel values between units (bq, bqml, suv).
    
    Parameters:
    - image: sitk.Image, input image.
    - from_unit: str, source unit ("bq", "bqml", "suv").
    - to_unit: str, target unit ("bq", "bqml", "suv").
    - weight: float, patient weight in kg (required for SUV conversions).
    - activity: float, injected activity in MBq (required for SUV conversions).
    """
    from_unit = from_unit.lower()
    to_unit = to_unit.lower()
    
    if from_unit == to_unit:
        return image

    # Validate weight and activity if SUV is involved
    if "suv" in (from_unit, to_unit):
        if weight is None or weight <= 0:
            raise ValueError("A positive --weight (in kg) is required when converting to or from SUV.")
        if activity is None or activity <= 0:
            raise ValueError("A positive --activity (in MBq) is required when converting to or from SUV.")

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

    return converted


def convert_image_units_file(
    image_path: Path,
    from_unit: str,
    to_unit: str,
    output_path: Optional[Path] = None,
    weight: Optional[float] = None,
    activity: Optional[float] = None,
) -> None:
    """
    Load image from path, convert its units, and save to output path.
    """
    if not image_path.exists():
        raise FileNotFoundError(f"Input image not found: {image_path}")

    from_unit = from_unit.lower()
    to_unit = to_unit.lower()

    if output_path is None:
        suffix = "".join(image_path.suffixes)
        base = image_path.name[:-len(suffix)] if suffix else image_path.name
        output_path = image_path.parent / f"{base}_{to_unit}{suffix}"

    print(f"Reading image from: {image_path}")
    image = sitk.ReadImage(str(image_path))

    print(f"Converting units from {from_unit} to {to_unit}...")
    converted_image = convert_image_units(
        image=image,
        from_unit=from_unit,
        to_unit=to_unit,
        weight=weight,
        activity=activity,
    )

    print(f"Saving converted image to: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(converted_image, str(output_path))
    print("Unit conversion completed successfully.")
