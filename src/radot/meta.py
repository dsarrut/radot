import sys
from pathlib import Path
from typing import Optional
import SimpleITK as sitk

METADATA_KEYS = {
    "modality": "modality",
    "unit": "unit",
    "weight": "weight",
    "activity": "activity"
}

def parse_description_metadata(desc: str) -> dict:
    """Parse key=value pairs from radot: description metadata string."""
    metadata = {}
    if desc.startswith("radot:"):
        parts = desc[len("radot:"):].split(";")
        for part in parts:
            if "=" in part:
                k, v = part.split("=", 1)
                k = k.strip().lower()
                v = v.strip()
                if v:
                    metadata[k] = v
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

    # Check if we are writing metadata
    is_writing = any(v is not None for v in (modality, unit, weight, activity))

    if is_writing:
        print("Writing metadata...")
        # Read existing metadata to merge
        existing_modality = get_metadata_case_insensitive(image, "modality")
        existing_unit = get_metadata_case_insensitive(image, "unit")
        existing_weight = get_metadata_case_insensitive(image, "weight")
        existing_activity = get_metadata_case_insensitive(image, "activity")
        
        # Merge with new values
        final_modality = modality.lower() if modality is not None else existing_modality
        final_unit = unit.lower() if unit is not None else existing_unit
        final_weight = str(weight) if weight is not None else existing_weight
        final_activity = str(activity) if activity is not None else existing_activity
        
        # Write individual keys (for NRRD / MHD)
        if final_modality is not None:
            image.SetMetaData("modality", final_modality)
        if final_unit is not None:
            image.SetMetaData("unit", final_unit)
        if final_weight is not None:
            image.SetMetaData("weight", final_weight)
        if final_activity is not None:
            image.SetMetaData("activity", final_activity)
            
        # Write to both descrip and ITK_FileNotes (for NIfTI)
        desc_parts = []
        if final_modality is not None:
            desc_parts.append(f"modality={final_modality}")
        if final_unit is not None:
            desc_parts.append(f"unit={final_unit}")
        if final_weight is not None:
            desc_parts.append(f"weight={final_weight}")
        if final_activity is not None:
            desc_parts.append(f"activity={final_activity}")
            
        desc_str = "radot:" + ";".join(desc_parts)
        if len(desc_str) > 80:
            print(f"Warning: Metadata string length ({len(desc_str)}) exceeds 80 characters. It will be truncated in NIfTI headers.")
            desc_str = desc_str[:80]
            
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
            
        print("----------------------\n")
