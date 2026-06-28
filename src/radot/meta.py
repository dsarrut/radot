import sys
import json
from pathlib import Path
import logging

logger = logging.getLogger("radot")
from typing import Optional
import SimpleITK as sitk
from datetime import datetime, timedelta

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


def parse_time_from_injection(t_str: str) -> timedelta:
    """Parse a time string formatted as 'day:hours:min' (e.g. '2:15:30') into a timedelta."""
    parts = t_str.split(":")
    if len(parts) != 3:
        raise ValueError(
            f"Invalid time_from_injection format: {t_str}. Expected 'day:hours:min' (e.g. 2:15:30)"
        )
    try:
        days = int(parts[0])
        hours = int(parts[1])
        minutes = int(parts[2])
        return timedelta(days=days, hours=hours, minutes=minutes)
    except ValueError:
        raise ValueError(
            f"Invalid time_from_injection format: {t_str}. Expected integers 'day:hours:min'"
        )


def get_sidecar_path(image_path: Path) -> Path:
    """Get the path to the JSON sidecar file matching the image file."""
    name_str = image_path.name
    # Strip common medical image extensions
    for ext in (".nii.gz", ".nii", ".nrrd", ".mhd", ".mha"):
        if name_str.lower().endswith(ext):
            base_name = name_str[: -len(ext)]
            return image_path.parent / f"{base_name}.json"
    return image_path.with_suffix(".json")


def get_metadata_case_insensitive(
    image: Optional[sitk.Image], key: str, image_path: Optional[Path] = None
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
                logger.warning(
                    f"Failed to read metadata from sidecar {sidecar_path}: {e}"
                )

    # 2. Try parsing headers using ImageFileReader (no pixel loading)
    if image_path is not None and image_path.exists():
        try:
            reader = sitk.ImageFileReader()
            reader.SetFileName(str(image_path))
            reader.ReadImageInformation()
            keys = reader.GetMetaDataKeys()
            for k in keys:
                if k.lower() == lower_key:
                    return reader.GetMetaData(k)
        except Exception:
            pass

    # 3. Fallback: Parse from loaded image object headers
    if image is not None:
        keys = image.GetMetaDataKeys()

        # Try case-insensitive lookup of individual keys
        for k in keys:
            if k.lower() == lower_key:
                return image.GetMetaData(k)

    return None


def get_pixel_size_bytes(pixel_id: int) -> int:
    """Return the byte size of a single component of the given pixel ID."""
    # Mapping of common SimpleITK pixel types to byte sizes
    mapping = {
        sitk.sitkInt8: 1,
        sitk.sitkUInt8: 1,
        sitk.sitkInt16: 2,
        sitk.sitkUInt16: 2,
        sitk.sitkInt32: 4,
        sitk.sitkUInt32: 4,
        sitk.sitkInt64: 8,
        sitk.sitkUInt64: 8,
        sitk.sitkFloat32: 4,
        sitk.sitkFloat64: 8,
    }
    return mapping.get(pixel_id, 1)


def format_memory_size(size_in_bytes: int) -> str:
    """Format memory size in a human-readable string."""
    val = float(size_in_bytes)
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if val < 1024.0:
            return f"{val:.2f} {unit}"
        val /= 1024.0
    return f"{val:.2f} PB"


def read_metadata(image_path: Path) -> dict:
    """
    Read metadata and image grid information from sidecar and file headers.
    Returns a dictionary structured as:
    {
        "image_info": { ... },
        "metadata": { ... }
    }
    """
    image_path = Path(image_path)
    if not image_path.exists():
        raise FileNotFoundError(f"Image file not found: {image_path}")

    # Gather image information (headers only)
    image_info = {}
    try:
        reader = sitk.ImageFileReader()
        reader.SetFileName(str(image_path))
        reader.ReadImageInformation()
        
        size = reader.GetSize()
        spacing = reader.GetSpacing()
        origin = reader.GetOrigin()
        dir_matrix = reader.GetDirection()
        
        pixel_type = sitk.GetPixelIDValueAsString(reader.GetPixelID())
        
        num_voxels = 1
        for s in size:
            num_voxels *= s

        pixel_size = get_pixel_size_bytes(reader.GetPixelID()) * reader.GetNumberOfComponents()
        mem_size = num_voxels * pixel_size
        mem_size_str = format_memory_size(mem_size)

        image_info = {
            "size": size,
            "spacing": spacing,
            "origin": origin,
            "direction": dir_matrix,
            "pixel_type": pixel_type,
            "num_voxels": num_voxels,
            "memory_size": mem_size_str
        }
    except Exception:
        pass

    # Gather therapy metadata
    modality = get_metadata_case_insensitive(None, "modality", image_path)
    unit = get_metadata_case_insensitive(None, "unit", image_path)
    weight = get_metadata_case_insensitive(None, "weight", image_path)
    activity = get_metadata_case_insensitive(None, "activity", image_path)
    injection_date = get_metadata_case_insensitive(None, "injection_date", image_path)
    acquisition_date = get_metadata_case_insensitive(None, "acquisition_date", image_path)

    # Convert numeric metadata types back if possible
    if weight is not None:
        try:
            weight = float(weight)
        except ValueError:
            pass
    if activity is not None:
        try:
            activity = float(activity)
        except ValueError:
            pass

    metadata = {
        "modality": modality,
        "unit": unit,
        "weight": weight,
        "activity": activity,
        "injection_date": injection_date,
        "acquisition_date": acquisition_date,
    }

    # Compute time from injection if dates are set
    if injection_date and acquisition_date:
        try:
            inj_dt = datetime.strptime(injection_date, "%Y-%m-%d %H:%M:%S")
            acq_dt = datetime.strptime(acquisition_date, "%Y-%m-%d %H:%M:%S")
            time_diff = acq_dt - inj_dt
            metadata["time_from_injection"] = time_diff.total_seconds() / 3600.0
        except Exception:
            pass

    return {
        "image_info": image_info,
        "metadata": metadata
    }


def write_metadata(
    image_path: Path,
    modality: Optional[str] = None,
    unit: Optional[str] = None,
    weight: Optional[float] = None,
    activity: Optional[float] = None,
    injection_date: Optional[str] = None,
    acquisition_date: Optional[str] = None,
    time_from_injection: Optional[str] = None,
) -> None:
    """
    Write or merge metadata fields in the JSON sidecar file.
    """
    image_path = Path(image_path)
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

    # Read existing metadata to merge
    existing_modality = get_metadata_case_insensitive(None, "modality", image_path)
    existing_unit = get_metadata_case_insensitive(None, "unit", image_path)
    existing_weight = get_metadata_case_insensitive(None, "weight", image_path)
    existing_activity = get_metadata_case_insensitive(None, "activity", image_path)
    existing_injection_date = get_metadata_case_insensitive(None, "injection_date", image_path)
    existing_acquisition_date = get_metadata_case_insensitive(None, "acquisition_date", image_path)

    # Merge with new values
    final_modality = modality.lower() if modality is not None else existing_modality
    final_unit = unit.lower() if unit is not None else existing_unit
    if final_modality == "ct" and final_unit is None:
        final_unit = "hu"
    if final_modality == "dose" and final_unit is None:
        final_unit = "gy"
    final_weight = (
        weight
        if weight is not None
        else (float(existing_weight) if existing_weight else None)
    )
    final_activity = (
        activity
        if activity is not None
        else (float(existing_activity) if existing_activity else None)
    )
    final_injection_date = (
        injection_date if injection_date is not None else existing_injection_date
    )
    final_acquisition_date = (
        acquisition_date
        if acquisition_date is not None
        else existing_acquisition_date
    )

    # Apply time_from_injection logic
    if time_from_injection is not None:
        delta = parse_time_from_injection(time_from_injection)
        if final_injection_date is not None and final_acquisition_date is not None:
            raise ValueError("Cannot specify --time-from-injection when both injection and acquisition dates are present.")
        elif final_injection_date is not None:
            inj_dt = datetime.strptime(final_injection_date, "%Y-%m-%d %H:%M:%S")
            acq_dt = inj_dt + delta
            final_acquisition_date = acq_dt.strftime("%Y-%m-%d %H:%M:%S")
        elif final_acquisition_date is not None:
            acq_dt = datetime.strptime(final_acquisition_date, "%Y-%m-%d %H:%M:%S")
            inj_dt = acq_dt - delta
            final_injection_date = inj_dt.strftime("%Y-%m-%d %H:%M:%S")
        else:
            final_injection_date = "2000-01-01 12:00:00"
            inj_dt = datetime.strptime(final_injection_date, "%Y-%m-%d %H:%M:%S")
            acq_dt = inj_dt + delta
            final_acquisition_date = acq_dt.strftime("%Y-%m-%d %H:%M:%S")

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
    logger.info(f"Saving metadata to sidecar: {sidecar_path}")
    with open(sidecar_path, "w") as f:
        json.dump(meta_dict, f, indent=4)
    logger.info("Metadata saved successfully.")


def print_metadata(info_dict: dict) -> None:
    """
    Print gathered metadata and image information to the console.
    """
    import sys
    image_info = info_dict.get("image_info", {})
    if image_info:
        size = image_info.get("size")
        sys.stdout.write(f"Image Size: {size}\n")
        sys.stdout.write(f"Spacing: {image_info.get('spacing')}\n")
        sys.stdout.write(f"Origin: {image_info.get('origin')}\n")
        
        dir_matrix = image_info.get("direction")
        if dir_matrix:
            if len(dir_matrix) == 9:
                matrix_str = (
                    f"[{dir_matrix[0]:.4f}, {dir_matrix[1]:.4f}, {dir_matrix[2]:.4f}]\n"
                    f"                    [{dir_matrix[3]:.4f}, {dir_matrix[4]:.4f}, {dir_matrix[5]:.4f}]\n"
                    f"                    [{dir_matrix[6]:.4f}, {dir_matrix[7]:.4f}, {dir_matrix[8]:.4f}]"
                )
            else:
                matrix_str = str(dir_matrix)
            sys.stdout.write(f"Orientation Matrix: {matrix_str}\n")
            
        sys.stdout.write(f"Pixel Type: {image_info.get('pixel_type')}\n")
        sys.stdout.write(f"Number of Voxels: {image_info.get('num_voxels')}\n")
        sys.stdout.write(f"Memory Size: {image_info.get('memory_size')}\n")

    metadata = info_dict.get("metadata", {})
    
    mod_val = metadata.get("modality")
    sys.stdout.write(f"Modality: {mod_val.upper() if mod_val else '(not set)'}\n")

    unit_val = metadata.get("unit")
    sys.stdout.write(f"Unit: {unit_val.lower() if unit_val else '(not set)'}\n")

    weight_val = metadata.get("weight")
    sys.stdout.write(f"Weight: {f'{weight_val} kg' if weight_val is not None else '(not set)'}\n")

    activity_val = metadata.get("activity")
    sys.stdout.write(f"Activity: {f'{activity_val} MBq' if activity_val is not None else '(not set)'}\n")

    inj_val = metadata.get("injection_date")
    sys.stdout.write(f"Injection Date: {inj_val if inj_val else '(not set)'}\n")

    acq_val = metadata.get("acquisition_date")
    sys.stdout.write(f"Acquisition Date: {acq_val if acq_val else '(not set)'}\n")

    time_from_injection = metadata.get("time_from_injection")
    if time_from_injection is not None:
        sys.stdout.write(f"Time from injection: {time_from_injection:.2f} hours\n")


def manage_metadata(
    image_path: Path,
    modality: Optional[str] = None,
    unit: Optional[str] = None,
    weight: Optional[float] = None,
    activity: Optional[float] = None,
    injection_date: Optional[str] = None,
    acquisition_date: Optional[str] = None,
    time_from_injection: Optional[str] = None,
) -> None:
    """
    Read or write metadata fields in the JSON sidecar file (wrapper for backward compatibility).
    """
    is_writing = any(
        v is not None
        for v in (modality, unit, weight, activity, injection_date, acquisition_date, time_from_injection)
    )
    if is_writing:
        write_metadata(
            image_path=image_path,
            modality=modality,
            unit=unit,
            weight=weight,
            activity=activity,
            injection_date=injection_date,
            acquisition_date=acquisition_date,
            time_from_injection=time_from_injection,
        )
    else:
        info_dict = read_metadata(image_path)
        print_metadata(info_dict)


def get_time_from_injection_h(
    image: sitk.Image, image_path: Optional[Path] = None
) -> float:
    """Retrieve or compute the time from injection in hours."""
    # 1. Try reading directly
    t_str = get_metadata_case_insensitive(image, "time_from_injection", image_path)
    if t_str is not None:
        try:
            return float(t_str)
        except ValueError:
            pass

    # 2. Try calculating from acquisition and injection dates
    inj_val = get_metadata_case_insensitive(image, "injection_date", image_path)
    acq_val = get_metadata_case_insensitive(image, "acquisition_date", image_path)
    if inj_val and acq_val:
        try:
            inj_dt = datetime.strptime(inj_val, "%Y-%m-%d %H:%M:%S")
            acq_dt = datetime.strptime(acq_val, "%Y-%m-%d %H:%M:%S")
            time_diff = acq_dt - inj_dt
            return time_diff.total_seconds() / 3600.0
        except Exception:
            pass

    raise ValueError(
        "Time from injection is required but was not found in the image metadata."
    )
