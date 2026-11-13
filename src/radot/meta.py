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
    "height": "height",
    "sex": "sex",
    "radionuclide": "radionuclide",
    "isotope": "radionuclide",
}


def compute_file_sha256(path: Path, chunk_size: int = 65536) -> Optional[str]:
    """Compute SHA-256 checksum of a file in chunks."""
    import hashlib

    p = Path(path)
    if not p.exists() or not p.is_file():
        return None
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while chunk := f.read(chunk_size):
            h.update(chunk)
    return h.hexdigest()


# The only date format of radot (sidecars, options, configs): no day-first, no date without time
DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def parse_date(value, name: str = "date") -> datetime:
    """A radot date: 'YYYY-MM-DD HH:MM:SS', optional fraction of second (or a datetime). Any other
    format (day first, slashes, date without time) is an error."""
    if isinstance(value, datetime):
        return value
    text = str(value).strip()
    for fmt in (DATE_FORMAT, DATE_FORMAT + ".%f"):
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            continue
    raise ValueError(f"{name} '{value}' is not in the format 'YYYY-MM-DD HH:MM:SS'.")


def metadata_number(value, key: str) -> Optional[float]:
    """A numeric sidecar value (weight, activity, height): None if absent, else a float; a value that is
    not a number (e.g. '70kg') is an error naming the key, never ignored nor kept as text."""
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        raise ValueError(
            f"{key} '{value}' in the metadata is not a number. Fix it (radot meta IMAGE --{key} ...)."
        ) from None


def validate_date(date_str: str) -> bool:
    """Validate that date_str is in 'YYYY-MM-DD HH:MM:SS' format."""
    try:
        parse_date(date_str)
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


def _read_sidecar(sidecar_path: Path) -> dict:
    """Content of a JSON sidecar. An unreadable sidecar is an error: its values are never replaced by
    the header ones or dropped when it is rewritten."""
    try:
        data = json.loads(Path(sidecar_path).read_text())
    except Exception as e:
        raise ValueError(
            f"Sidecar {sidecar_path} cannot be read ({e}). Fix it, or delete it to start again."
        ) from e
    if not isinstance(data, dict):
        raise ValueError(f"Sidecar {sidecar_path} is not a JSON object ({{\"key\": value, ...}}).")
    return data


def get_metadata_case_insensitive(
    image: Optional[sitk.Image], key: str, image_path: Optional[Path] = None
) -> Optional[str]:
    """Retrieve metadata value case-insensitively, preferring JSON sidecar and falling back to image headers."""
    lower_key = key.lower()

    # 1. Try sidecar JSON first
    if image_path is not None:
        sidecar_path = get_sidecar_path(image_path)
        if sidecar_path.exists():
            data = _read_sidecar(sidecar_path)
            # Case-insensitive lookup in JSON dict
            for k, v in data.items():
                if k.lower() == lower_key:
                    return str(v) if v is not None else None

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
    for unit in ["B", "KB", "MB", "GB", "TB"]:
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

        pixel_size = (
            get_pixel_size_bytes(reader.GetPixelID()) * reader.GetNumberOfComponents()
        )
        mem_size = num_voxels * pixel_size
        mem_size_str = format_memory_size(mem_size)

        image_info = {
            "size": size,
            "spacing": spacing,
            "origin": origin,
            "direction": dir_matrix,
            "pixel_type": pixel_type,
            "num_voxels": num_voxels,
            "memory_size": mem_size_str,
        }
    except Exception:
        pass

    # Gather therapy metadata
    modality = get_metadata_case_insensitive(None, "modality", image_path)
    unit = get_metadata_case_insensitive(None, "unit", image_path)
    weight = get_metadata_case_insensitive(None, "weight", image_path)
    activity = get_metadata_case_insensitive(None, "activity", image_path)
    injection_date = get_metadata_case_insensitive(None, "injection_date", image_path)
    acquisition_date = get_metadata_case_insensitive(
        None, "acquisition_date", image_path
    )
    height = get_metadata_case_insensitive(None, "height", image_path)
    sex = get_metadata_case_insensitive(None, "sex", image_path)
    radionuclide = get_metadata_case_insensitive(None, "radionuclide", image_path)
    decay_correction = get_metadata_case_insensitive(None, "decay_correction", image_path)

    # Convert numeric metadata types back if possible
    # Numbers; a value that is not a number is kept as text for display, with a WARNING (radot check
    # reports it, the tools that use it raise)
    def _num(key, value):
        try:
            return metadata_number(value, key)
        except ValueError as e:
            logger.warning(f"{image_path.name}: {e}")
            return value

    weight = _num("weight", weight)
    activity = _num("activity", activity)
    height = _num("height", height)
    if sex is not None:
        sex = sex.lower().strip()

    metadata = {
        "modality": modality,
        "unit": unit,
        "weight": weight,
        "activity": activity,
        "injection_date": injection_date,
        "acquisition_date": acquisition_date,
        "height": height,
        "sex": sex,
        "radionuclide": radionuclide,
        "decay_correction": decay_correction,
    }

    # Compute time from injection if dates are set
    if injection_date and acquisition_date:
        try:
            inj_dt = parse_date(injection_date)
            acq_dt = parse_date(acquisition_date)
            time_diff = acq_dt - inj_dt
            metadata["time_from_injection"] = time_diff.total_seconds() / 3600.0
        except Exception:
            pass

    return {"image_info": image_info, "metadata": metadata}


def write_metadata(
    image_path: Path,
    modality: Optional[str] = None,
    unit: Optional[str] = None,
    weight: Optional[float] = None,
    activity: Optional[float] = None,
    injection_date: Optional[str] = None,
    acquisition_date: Optional[str] = None,
    time_from_injection: Optional[str] = None,
    height: Optional[float] = None,
    sex: Optional[str] = None,
    radionuclide: Optional[str] = None,
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

    # Validate sex if provided
    if sex is not None:
        norm_sex = sex.lower().strip()
        if norm_sex not in ("male", "female"):
            raise ValueError(f"Invalid sex value: {sex}. Expected 'male' or 'female'.")
        sex = norm_sex

    # Validate radionuclide if provided (known to ICRP 107)
    if radionuclide is not None:
        from radot.db import get_physical_half_life_h

        radionuclide = radionuclide.strip()
        get_physical_half_life_h(radionuclide)  # raises ValueError if unknown

    # Read existing metadata to merge
    existing_modality = get_metadata_case_insensitive(None, "modality", image_path)
    existing_unit = get_metadata_case_insensitive(None, "unit", image_path)
    existing_weight = get_metadata_case_insensitive(None, "weight", image_path)
    existing_activity = get_metadata_case_insensitive(None, "activity", image_path)
    existing_injection_date = get_metadata_case_insensitive(
        None, "injection_date", image_path
    )
    existing_acquisition_date = get_metadata_case_insensitive(
        None, "acquisition_date", image_path
    )
    existing_height = get_metadata_case_insensitive(None, "height", image_path)
    existing_sex = get_metadata_case_insensitive(None, "sex", image_path)

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
        else metadata_number(existing_weight, "weight")
    )
    final_activity = (
        activity
        if activity is not None
        else metadata_number(existing_activity, "activity")
    )
    final_injection_date = (
        injection_date if injection_date is not None else existing_injection_date
    )
    final_acquisition_date = (
        acquisition_date if acquisition_date is not None else existing_acquisition_date
    )
    final_height = (
        height
        if height is not None
        else metadata_number(existing_height, "height")
    )
    final_sex = (
        sex.lower().strip()
        if sex is not None
        else (existing_sex.lower().strip() if existing_sex else None)
    )

    # Apply time_from_injection logic
    if time_from_injection is not None:
        delta = parse_time_from_injection(time_from_injection)
        if final_injection_date is not None and final_acquisition_date is not None:
            raise ValueError(
                "Cannot specify --time-from-injection when both injection and acquisition dates are present."
            )
        elif final_injection_date is not None:
            inj_dt = parse_date(final_injection_date)
            acq_dt = inj_dt + delta
            final_acquisition_date = acq_dt.strftime("%Y-%m-%d %H:%M:%S")
        elif final_acquisition_date is not None:
            acq_dt = parse_date(final_acquisition_date)
            inj_dt = acq_dt - delta
            final_injection_date = inj_dt.strftime("%Y-%m-%d %H:%M:%S")
        else:
            final_injection_date = "2000-01-01 12:00:00"
            inj_dt = parse_date(final_injection_date)
            acq_dt = inj_dt + delta
            final_acquisition_date = acq_dt.strftime("%Y-%m-%d %H:%M:%S")

    # Write to JSON sidecar, keeping the keys not managed here (e.g. radionuclide, decay_correction)
    sidecar_path = get_sidecar_path(image_path)
    meta_dict = _read_sidecar(sidecar_path) if sidecar_path.exists() else {}
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
    if final_height is not None:
        meta_dict["height"] = final_height
    if final_sex is not None:
        meta_dict["sex"] = final_sex
    if radionuclide is not None:
        meta_dict["radionuclide"] = radionuclide

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
    sys.stdout.write(
        f"Weight: {f'{weight_val} kg' if weight_val is not None else '(not set)'}\n"
    )

    activity_val = metadata.get("activity")
    sys.stdout.write(
        f"Activity: {f'{activity_val} MBq' if activity_val is not None else '(not set)'}\n"
    )

    inj_val = metadata.get("injection_date")
    sys.stdout.write(f"Injection Date: {inj_val if inj_val else '(not set)'}\n")

    acq_val = metadata.get("acquisition_date")
    sys.stdout.write(f"Acquisition Date: {acq_val if acq_val else '(not set)'}\n")

    time_from_injection = metadata.get("time_from_injection")
    if time_from_injection is not None:
        sys.stdout.write(f"Time from injection: {time_from_injection:.2f} hours\n")

    height_val = metadata.get("height")
    sys.stdout.write(
        f"Height: {f'{height_val} cm' if height_val is not None else '(not set)'}\n"
    )

    sex_val = metadata.get("sex")
    sys.stdout.write(f"Sex: {sex_val if sex_val else '(not set)'}\n")

    nuclide_val = metadata.get("radionuclide")
    sys.stdout.write(f"Radionuclide: {nuclide_val if nuclide_val else '(not set)'}\n")

    decay_val = metadata.get("decay_correction")
    if decay_val:
        sys.stdout.write(f"Decay Correction: {decay_val}\n")


def manage_metadata(
    image_path: Path,
    modality: Optional[str] = None,
    unit: Optional[str] = None,
    weight: Optional[float] = None,
    activity: Optional[float] = None,
    injection_date: Optional[str] = None,
    acquisition_date: Optional[str] = None,
    time_from_injection: Optional[str] = None,
    height: Optional[float] = None,
    sex: Optional[str] = None,
    radionuclide: Optional[str] = None,
) -> None:
    """
    Read or write metadata fields in the JSON sidecar file (wrapper for backward compatibility).
    """
    is_writing = any(
        v is not None
        for v in (
            modality,
            unit,
            weight,
            activity,
            injection_date,
            acquisition_date,
            time_from_injection,
            height,
            sex,
            radionuclide,
        )
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
            height=height,
            sex=sex,
            radionuclide=radionuclide,
        )
    else:
        info_dict = read_metadata(image_path)
        print_metadata(info_dict)


def get_time_from_injection_h(
    image: sitk.Image, image_path: Optional[Path] = None
) -> float:
    """Time from injection in hours: from the injection and acquisition dates, or from a
    `time_from_injection` key (hours). Both present and different, or one malformed, is an error."""
    t_str = get_metadata_case_insensitive(image, "time_from_injection", image_path)
    t_key = None
    if t_str is not None:
        try:
            t_key = float(t_str)
        except ValueError:
            raise ValueError(
                f"time_from_injection '{t_str}' is not a number of hours."
            ) from None

    t_dates = None
    inj_val = get_metadata_case_insensitive(image, "injection_date", image_path)
    acq_val = get_metadata_case_insensitive(image, "acquisition_date", image_path)
    if inj_val and acq_val:
        inj_dt = parse_date(inj_val, "injection_date")
        acq_dt = parse_date(acq_val, "acquisition_date")
        t_dates = (acq_dt - inj_dt).total_seconds() / 3600.0

    if t_key is not None and t_dates is not None and abs(t_key - t_dates) > 1.0 / 60.0:
        raise ValueError(
            f"time_from_injection ({t_key:g} h) differs from the dates ({t_dates:.4g} h: injection "
            f"{inj_val}, acquisition {acq_val}). Remove one of them."
        )
    if t_dates is not None:
        return t_dates
    if t_key is not None:
        return t_key

    raise ValueError(
        "Time from injection is required but was not found in the image metadata."
    )
