import sys
import json
from pathlib import Path
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
                print(
                    f"Warning: Failed to read metadata from sidecar {sidecar_path}: {e}"
                )

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
        for v in (modality, unit, weight, activity, injection_date, acquisition_date, time_from_injection)
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
        existing_injection_date = get_metadata_case_insensitive(
            image, "injection_date", image_path
        )
        existing_acquisition_date = get_metadata_case_insensitive(
            image, "acquisition_date", image_path
        )

        # Merge with new values
        final_modality = modality.lower() if modality is not None else existing_modality
        final_unit = unit.lower() if unit is not None else existing_unit
        if final_modality == "ct" and final_unit is None:
            final_unit = "hu"
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
        print(f"Saving metadata to sidecar: {sidecar_path}")
        with open(sidecar_path, "w") as f:
            json.dump(meta_dict, f, indent=4)
        print("Metadata saved successfully.")
    else:
        # Just reading metadata
        try:
            image = sitk.ReadImage(str(image_path))
        except Exception:
            image = None

        if image is not None:
            size = image.GetSize()
            spacing = image.GetSpacing()
            origin = image.GetOrigin()
            dir_matrix = image.GetDirection()
            if len(dir_matrix) == 9:
                matrix_str = (
                    f"[{dir_matrix[0]:.4f}, {dir_matrix[1]:.4f}, {dir_matrix[2]:.4f}]\n"
                    f"                    [{dir_matrix[3]:.4f}, {dir_matrix[4]:.4f}, {dir_matrix[5]:.4f}]\n"
                    f"                    [{dir_matrix[6]:.4f}, {dir_matrix[7]:.4f}, {dir_matrix[8]:.4f}]"
                )
            else:
                matrix_str = str(dir_matrix)

            pixel_type = image.GetPixelIDTypeAsString()
            
            # Number of voxels
            num_voxels = 1
            for s in size:
                num_voxels *= s

            pixel_size = get_pixel_size_bytes(image.GetPixelID()) * image.GetNumberOfComponentsPerPixel()
            mem_size = num_voxels * pixel_size
            mem_size_str = format_memory_size(mem_size)

            print(f"Image Size: {size}")
            print(f"Spacing: {spacing}")
            print(f"Origin: {origin}")
            print(f"Orientation Matrix: {matrix_str}")
            print(f"Pixel Type: {pixel_type}")
            print(f"Number of Voxels: {num_voxels}")
            print(f"Memory Size: {mem_size_str}")

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
