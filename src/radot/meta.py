import sys
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


def get_metadata_case_insensitive(image: sitk.Image, key: str) -> Optional[str]:
    """Retrieve metadata value case-insensitively, falling back to NIfTI description field."""
    keys = image.GetMetaDataKeys()
    lower_key = key.lower()

    # 1. Try case-insensitive lookup of individual keys
    for k in keys:
        if k.lower() == lower_key:
            return image.GetMetaData(k)

    # 2. Fallback: Parse from NIfTI description fields (keys 'descrip' or 'itk_filenotes')
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
    Read or write metadata fields in the given image.
    If no new metadata is specified, prints existing metadata to console.
    """
    if not image_path.exists():
        raise FileNotFoundError(f"Image file not found: {image_path}")

    # Read image
    print(f"Reading image: {image_path}")
    image = sitk.ReadImage(str(image_path))

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
        # Read existing metadata to merge
        existing_modality = get_metadata_case_insensitive(image, "modality")
        existing_unit = get_metadata_case_insensitive(image, "unit")
        existing_weight = get_metadata_case_insensitive(image, "weight")
        existing_activity = get_metadata_case_insensitive(image, "activity")
        existing_injection_date = get_metadata_case_insensitive(image, "injection_date")
        existing_acquisition_date = get_metadata_case_insensitive(
            image, "acquisition_date"
        )

        # Merge with new values
        final_modality = modality.lower() if modality is not None else existing_modality
        final_unit = unit.lower() if unit is not None else existing_unit
        final_weight = str(weight) if weight is not None else existing_weight
        final_activity = str(activity) if activity is not None else existing_activity
        final_injection_date = (
            injection_date if injection_date is not None else existing_injection_date
        )
        final_acquisition_date = (
            acquisition_date
            if acquisition_date is not None
            else existing_acquisition_date
        )

        # Write individual keys (for NRRD / MHD)
        if final_modality is not None:
            image.SetMetaData("modality", final_modality)
        if final_unit is not None:
            image.SetMetaData("unit", final_unit)
        if final_weight is not None:
            image.SetMetaData("weight", final_weight)
        if final_activity is not None:
            image.SetMetaData("activity", final_activity)
        if final_injection_date is not None:
            image.SetMetaData("injection_date", final_injection_date)
        if final_acquisition_date is not None:
            image.SetMetaData("acquisition_date", final_acquisition_date)

        # Write to both descrip and ITK_FileNotes (for NIfTI)
        desc_parts = []
        if final_modality is not None:
            desc_parts.append(f"m={final_modality}")
        if final_unit is not None:
            desc_parts.append(f"u={final_unit}")
        if final_weight is not None:
            desc_parts.append(f"w={final_weight}")
        if final_activity is not None:
            desc_parts.append(f"a={final_activity}")
        if final_injection_date is not None:
            desc_parts.append(f"i={final_injection_date}")
        if final_acquisition_date is not None:
            desc_parts.append(f"d={final_acquisition_date}")

        is_nifti = image_path.name.lower().endswith((".nii", ".nii.gz"))
        desc_str = "radot:" + ";".join(desc_parts)
        if is_nifti:
            if len(desc_str) > 79:
                print(
                    f"Warning: Metadata string length ({len(desc_str)}) exceeds 79 characters. It will be truncated in NIfTI headers to prevent writer errors."
                )
                desc_str = desc_str[:79]
        image.SetMetaData("descrip", desc_str)
        image.SetMetaData("ITK_FileNotes", desc_str)

        # Save back to same file
        sitk.WriteImage(image, str(image_path))
        print("Metadata saved successfully.")
    else:
        # Just reading metadata
        print("\n--- Image Metadata ---")
        found_any = False

        # Modality
        mod_val = get_metadata_case_insensitive(image, "modality")
        if mod_val:
            print(f"Modality: {mod_val.upper()}")
            found_any = True
        else:
            print("Modality: (not set)")

        # Unit
        unit_val = get_metadata_case_insensitive(image, "unit")
        if unit_val:
            print(f"Unit: {unit_val.lower()}")
            found_any = True
        else:
            print("Unit: (not set)")

        # Weight
        weight_val = get_metadata_case_insensitive(image, "weight")
        if weight_val:
            print(f"Weight: {weight_val} kg")
            found_any = True
        else:
            print("Weight: (not set)")

        # Activity
        activity_val = get_metadata_case_insensitive(image, "activity")
        if activity_val:
            print(f"Activity: {activity_val} MBq")
            found_any = True
        else:
            print("Activity: (not set)")

        # Injection Date
        inj_val = get_metadata_case_insensitive(image, "injection_date")
        if inj_val:
            print(f"Injection Date: {inj_val}")
            found_any = True
        else:
            print("Injection Date: (not set)")

        # Acquisition Date
        acq_val = get_metadata_case_insensitive(image, "acquisition_date")
        if acq_val:
            print(f"Acquisition Date: {acq_val}")
            found_any = True
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
