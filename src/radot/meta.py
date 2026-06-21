import sys
import json
from pathlib import Path
from typing import Optional
import SimpleITK as sitk
from datetime import datetime

METADATA_KEYS = {
    "modality": "modality",
    "unit": "unit",
    "weight": "weight",
    "activity": "activity",
    "injection_date": "injection_date",
    "acquisition_date": "acquisition_date",
}

def validate_date(date_str: str) -> bool:
    """Validate that date_str is in 'YYYY-MM-DD HH:MM:SS' format."""
    try:
        datetime.strptime(date_str, "%Y-%m-%d %H:%M:%S")
        return True
    except ValueError:
        return False

def get_sidecar_path(image_path: Path) -> Path:
    """Get the path to the JSON sidecar file matching the image file."""
    name_str = image_path.name
    # Strip common medical image extensions
    for ext in ('.nii.gz', '.nii', '.nrrd', '.mhd', '.mha'):
        if name_str.lower().endswith(ext):
            base_name = name_str[:-len(ext)]
            return image_path.parent / f"{base_name}.json"
    return image_path.with_suffix(".json")

def parse_description_metadata(desc: str) -> dict:
    """Parse key=value pairs from radot: description metadata string."""
    metadata = {}
    
    key_mapping = {
        "m": "modality",
        "u": "unit",
        "w": "weight",
        "a": "activity",
        "i": "injection_date",
        "d": "acquisition_date"
    }
    
    if desc.startswith("radot:"):
        parts = desc[len("radot:") :].split(";")
        for part in parts:
            if "=" in part:
                k, v = part.split("=", 1)
                k = k.strip().lower()
                v = v.strip()
                if v:
                    resolved_key = key_mapping.get(k, k)
                    metadata[resolved_key] = v
    return metadata

def get_metadata_case_insensitive(
    image: Optional[sitk.Image],
    key: str,
    image_path: Optional[Path] = None
) -> Optional[str]:
    """Retrieve metadata value case-insensitively, preferring JSON sidecar and falling back to image headers."""
    lower_key = key.lower()

    # 1. Try sidecar JSON first
    if image_path is not None:
        sidecar_path = get_sidecar_path(image_path)
        if sidecar_path.exists():
            try:
                with open(sidecar_path, "r") as f:
                    data = json.load(f)
                # Case-insensitive lookup in JSON dict
                for k, v in data.items():
                    if k.lower() == lower_key:
                        return str(v) if v is not None else None
            except Exception as e:
                print(f"Warning: Failed to read metadata from sidecar {sidecar_path}: {e}")

    # 2. Try loading image if not provided but image_path is
    if image is None and image_path is not None and image_path.exists():
        try:
            image = sitk.ReadImage(str(image_path))
        except Exception:
            pass

    # 3. Fallback: Parse from loaded image object headers
    if image is not None:
        keys = image.GetMetaDataKeys()
        
        # Try case-insensitive lookup of individual keys
        for k in keys:
            if k.lower() == lower_key:
                return image.GetMetaData(k)
                
        # Try description field fallback
        desc_val = None
        for k in keys:
            if k.lower() in ("descrip", "itk_filenotes"):
                desc_val = image.GetMetaData(k)
                break
                
        if desc_val:
            parsed = parse_description_metadata(desc_val)
            return parsed.get(lower_key)

    return None

def manage_metadata(
    image_path: Path,
    modality: Optional[str] = None,
    unit: Optional[str] = None,
    weight: Optional[float] = None,
    activity: Optional[float] = None,
    injection_date: Optional[str] = None,
    acquisition_date: Optional[str] = None,
) -> None:
    """
    Read or write metadata fields in the JSON sidecar file.
    If no new metadata is specified, prints existing metadata to console.
    """
    if not image_path.exists():
        raise FileNotFoundError(f"Image file not found: {image_path}")

    # Validate dates if provided
    if injection_date is not None and not validate_date(injection_date):
        raise ValueError(
            f"Invalid injection_date format: {injection_date}. Expected 'YYYY-MM-DD HH:MM:SS'"
        )
    if acquisition_date is not None and not validate_date(acquisition_date):
        raise ValueError(
            f"Invalid acquisition_date format: {acquisition_date}. Expected 'YYYY-MM-DD HH:MM:SS'"
        )

    # Check if we are writing metadata
    is_writing = any(
        v is not None
        for v in (modality, unit, weight, activity, injection_date, acquisition_date)
    )

    if is_writing:
        print("Writing metadata...")
        # Read existing metadata to merge (using our helper that checks sidecars and falls back to headers)
        try:
            image = sitk.ReadImage(str(image_path))
        except Exception:
            image = None

        existing_modality = get_metadata_case_insensitive(image, "modality", image_path)
        existing_unit = get_metadata_case_insensitive(image, "unit", image_path)
        existing_weight = get_metadata_case_insensitive(image, "weight", image_path)
        existing_activity = get_metadata_case_insensitive(image, "activity", image_path)
        existing_injection_date = get_metadata_case_insensitive(image, "injection_date", image_path)
        existing_acquisition_date = get_metadata_case_insensitive(image, "acquisition_date", image_path)

        # Merge with new values
        final_modality = modality.lower() if modality is not None else existing_modality
        final_unit = unit.lower() if unit is not None else existing_unit
        final_weight = weight if weight is not None else (float(existing_weight) if existing_weight else None)
        final_activity = activity if activity is not None else (float(existing_activity) if existing_activity else None)
        final_injection_date = injection_date if injection_date is not None else existing_injection_date
        final_acquisition_date = acquisition_date if acquisition_date is not None else existing_acquisition_date

        # Write to JSON sidecar
        meta_dict = {}
        if final_modality is not None:
            meta_dict["modality"] = final_modality
        if final_unit is not None:
            meta_dict["unit"] = final_unit
        if final_weight is not None:
            meta_dict["weight"] = final_weight
        if final_activity is not None:
            meta_dict["activity"] = final_activity
        if final_injection_date is not None:
            meta_dict["injection_date"] = final_injection_date
        if final_acquisition_date is not None:
            meta_dict["acquisition_date"] = final_acquisition_date

        sidecar_path = get_sidecar_path(image_path)
        print(f"Saving metadata to sidecar: {sidecar_path}")
        with open(sidecar_path, "w") as f:
            json.dump(meta_dict, f, indent=4)
        print("Metadata saved successfully.")
    else:
        # Just reading metadata
        print("\n--- Image Metadata ---")
        try:
            image = sitk.ReadImage(str(image_path))
        except Exception:
            image = None

        mod_val = get_metadata_case_insensitive(image, "modality", image_path)
        if mod_val:
            print(f"Modality: {mod_val.upper()}")
        else:
            print("Modality: (not set)")

        unit_val = get_metadata_case_insensitive(image, "unit", image_path)
        if unit_val:
            print(f"Unit: {unit_val.lower()}")
        else:
            print("Unit: (not set)")

        weight_val = get_metadata_case_insensitive(image, "weight", image_path)
        if weight_val:
            print(f"Weight: {weight_val} kg")
        else:
            print("Weight: (not set)")

        activity_val = get_metadata_case_insensitive(image, "activity", image_path)
        if activity_val:
            print(f"Activity: {activity_val} MBq")
        else:
            print("Activity: (not set)")

        inj_val = get_metadata_case_insensitive(image, "injection_date", image_path)
        if inj_val:
            print(f"Injection Date: {inj_val}")
        else:
            print("Injection Date: (not set)")

        acq_val = get_metadata_case_insensitive(image, "acquisition_date", image_path)
        if acq_val:
            print(f"Acquisition Date: {acq_val}")
        else:
            print("Acquisition Date: (not set)")

        # Time from Injection (computed)
        if inj_val and acq_val:
            try:
                inj_dt = datetime.strptime(inj_val, "%Y-%m-%d %H:%M:%S")
                acq_dt = datetime.strptime(acq_val, "%Y-%m-%d %H:%M:%S")
                time_diff = acq_dt - inj_dt
                print(f"Time from injection: {time_diff}")
            except Exception:
                pass

        print("----------------------\n")
