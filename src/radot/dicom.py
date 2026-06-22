import json
from pathlib import Path
from typing import Optional
import SimpleITK as sitk

def clean_tag(val: Optional[str]) -> Optional[str]:
    if val is None:
        return None
    val = val.strip()
    return val if val else None

def parse_date(date_str: str) -> Optional[str]:
    date_str = date_str.strip()
    if len(date_str) >= 8:
        return f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:8]}"
    return None

def parse_time(time_str: str) -> Optional[str]:
    time_str = time_str.strip()
    if len(time_str) >= 6:
        return f"{time_str[:2]}:{time_str[2:4]}:{time_str[4:6]}"
    elif len(time_str) >= 4:
        return f"{time_str[:2]}:{time_str[2:4]}:00"
    return None

def parse_datetime(dt_str: str) -> Optional[str]:
    dt_str = dt_str.strip()
    if len(dt_str) >= 14:
        return f"{dt_str[:4]}-{dt_str[4:6]}-{dt_str[6:8]} {dt_str[8:10]}:{dt_str[10:12]}:{dt_str[12:14]}"
    elif len(dt_str) >= 12:
        return f"{dt_str[:4]}-{dt_str[4:6]}-{dt_str[6:8]} {dt_str[8:10]}:{dt_str[10:12]}:00"
    elif len(dt_str) >= 8:
        date_part = parse_date(dt_str[:8])
        if date_part:
            return f"{date_part} 00:00:00"
    return None

def combine_date_time(date_str: Optional[str], time_str: Optional[str]) -> Optional[str]:
    if not date_str:
        return None
    d = parse_date(date_str)
    if not d:
        return None
    t = parse_time(time_str) if time_str else "00:00:00"
    if not t:
        t = "00:00:00"
    return f"{d} {t}"

def convert_dicom(dicom_path: Path, output_path: Optional[Path] = None) -> None:
    """
    Convert a DICOM directory (series) or single DICOM file into another format (e.g. NIfTI).
    Extracts DICOM metadata and writes it to a matching JSON sidecar file.
    If output_path is None, prints the extracted DICOM metadata to the console.
    """
    if not dicom_path.exists():
        raise FileNotFoundError(f"DICOM path not found: {dicom_path}")

    # Determine if directory (series) or single file
    is_dir = dicom_path.is_dir()

    if is_dir:
        print(f"Scanning directory for DICOM series: {dicom_path}")
        reader = sitk.ImageSeriesReader()
        dicom_names = reader.GetGDCMSeriesFileNames(str(dicom_path))
        if not dicom_names:
            raise ValueError(f"No DICOM series found in directory: {dicom_path}")
        header_file = dicom_names[0]
    else:
        print(f"Reading single DICOM file: {dicom_path}")
        header_file = str(dicom_path)

    # Read DICOM tags using ImageFileReader to be robust and header-only
    header_reader = sitk.ImageFileReader()
    header_reader.SetFileName(header_file)
    header_reader.ReadImageInformation()
    
    def get_tag(tag_str: str) -> Optional[str]:
        if header_reader.HasMetaDataKey(tag_str):
            return clean_tag(header_reader.GetMetaData(tag_str))
        return None

    # Parse modality
    modality_raw = get_tag("0008|0060")
    modality = None
    if modality_raw:
        m_lower = modality_raw.lower()
        if m_lower in ("ct", "pt", "pet", "nm", "spect"):
            if m_lower == "pt":
                modality = "pet"
            elif m_lower == "nm":
                modality = "spect"
            else:
                modality = m_lower
        else:
            modality = m_lower

    # Parse weight
    weight_raw = get_tag("0010|1030")
    weight = None
    if weight_raw:
        try:
            weight = float(weight_raw)
        except ValueError:
            pass

    # Parse activity (Radionuclide Total Dose 0018|1074 is in Bq, we store in MBq)
    activity_raw = get_tag("0018|1074")
    activity = None
    if activity_raw:
        try:
            activity = float(activity_raw) / 1_000_000.0
        except ValueError:
            pass

    # Parse unit (Units 0054|1001)
    unit_raw = get_tag("0054|1001")
    unit = None
    if unit_raw:
        u_upper = unit_raw.upper().strip()
        if "BQML" in u_upper:
            unit = "bqml"
        elif "GML" in u_upper:
            unit = "suv"
        else:
            unit = u_upper.lower()

    # Default unit fallbacks based on modality
    if unit is None:
        if modality == "ct":
            unit = "hu"
        elif modality in ("pet", "spect"):
            # If weight and activity are available, default to suv, otherwise bqml
            if weight is not None and activity is not None:
                unit = "suv"
            else:
                unit = "bqml"

    # Study/Series/Acquisition Dates and Times
    study_date = get_tag("0008|0020")
    study_time = get_tag("0008|0030")
    series_date = get_tag("0008|0021")
    series_time = get_tag("0008|0031")
    acq_date = get_tag("0008|0022")
    acq_time = get_tag("0008|0032")
    acq_datetime = get_tag("0008|002a")

    # Resolve acquisition_date
    acquisition_date = None
    if acq_datetime:
        acquisition_date = parse_datetime(acq_datetime)
    if not acquisition_date:
        acquisition_date = combine_date_time(acq_date, acq_time)
    if not acquisition_date:
        acquisition_date = combine_date_time(series_date, series_time)
    if not acquisition_date:
        acquisition_date = combine_date_time(study_date, study_time)

    # Radiopharmaceutical details for injection_date
    rad_datetime = get_tag("0018|1078")  # Radiopharmaceutical Start DateTime
    rad_time = get_tag("0018|1072")      # Radiopharmaceutical Start Time
    
    injection_date = None
    if rad_datetime:
        injection_date = parse_datetime(rad_datetime)
    if not injection_date and rad_time:
        # Fall back to matching acquisition/series/study date with rad_time
        for d in (acq_date, series_date, study_date):
            if d:
                injection_date = combine_date_time(d, rad_time)
                if injection_date:
                    break

    if output_path is None:
        tags_to_print = [
            ("0010|0010", "Patient's Name"),
            ("0010|0020", "Patient ID"),
            ("0010|1030", "Patient's Weight"),
            ("0008|0060", "Modality"),
            ("0008|1030", "Study Description"),
            ("0008|103e", "Series Description"),
            ("0008|0070", "Manufacturer"),
            ("0008|1090", "Manufacturer's Model Name"),
            ("0020|000d", "Study Instance UID"),
            ("0020|000e", "Series Instance UID"),
            ("0054|1001", "Units"),
            ("0018|1074", "Radionuclide Total Dose"),
            ("0018|1078", "Radiopharmaceutical Start DateTime"),
            ("0018|1072", "Radiopharmaceutical Start Time"),
            ("0008|0020", "Study Date"),
            ("0008|0030", "Study Time"),
            ("0008|0021", "Series Date"),
            ("0008|0031", "Series Time"),
            ("0008|0022", "Acquisition Date"),
            ("0008|0032", "Acquisition Time"),
            ("0008|002a", "Acquisition DateTime"),
        ]
        for key, name in tags_to_print:
            val = get_tag(key)
            tag_display = f"({key.replace('|', ',')})"
            print(f"{tag_display} {name}: {val if val is not None else '(not set)'}")
        return

    # Load full image only if saving to path
    if is_dir:
        reader.SetFileNames(dicom_names)
        reader.MetaDataDictionaryArrayUpdateOn()
        print(f"Reading DICOM series ({len(dicom_names)} slices)...")
        image = reader.Execute()
    else:
        image = sitk.ReadImage(str(dicom_path))

    # Build sidecar dictionary
    meta_dict = {}
    if modality:
        meta_dict["modality"] = modality
    if unit:
        meta_dict["unit"] = unit
    if weight is not None:
        meta_dict["weight"] = weight
    if activity is not None:
        meta_dict["activity"] = activity
    if injection_date:
        meta_dict["injection_date"] = injection_date
    if acquisition_date:
        meta_dict["acquisition_date"] = acquisition_date

    # Write output image
    print(f"Saving converted image to: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(image, str(output_path))

    # Write output sidecar JSON
    from radot.meta import get_sidecar_path
    sidecar_path = get_sidecar_path(output_path)
    print(f"Saving metadata sidecar to: {sidecar_path}")
    with open(sidecar_path, "w") as f:
        json.dump(meta_dict, f, indent=4)

    print("DICOM conversion completed successfully.")
