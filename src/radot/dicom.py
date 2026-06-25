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

def convert_dicom(
    dicom_path: Path,
    output_path: Optional[Path] = None,
    all_tags: bool = False
) -> None:
    """
    Convert a DICOM directory (series) or single DICOM file into another format (e.g. NIfTI).
    Extracts DICOM metadata and writes it to a matching JSON sidecar file.
    If output_path is None, prints the extracted DICOM metadata to the console.
    """
    if not dicom_path.exists():
        raise FileNotFoundError(f"DICOM path not found: {dicom_path}")

    if output_path is not None:
        if output_path.suffix.lower() in (".dcm", ".dicom"):
            raise ValueError("Writing back to DICOM format is not supported. Please use NIfTI (.nii.gz), NRRD (.nrrd), or MetaImage (.mhd) for the output file.")

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

    # Parse image type for sinogram detection
    image_type_raw = get_tag("0008|0008")
    is_sinogram = False
    if image_type_raw:
        parts = [p.strip().upper() for p in image_type_raw.replace("\\", "/").replace(",", "/").split("/")]
        if "TOMO" in parts and not any("RECON" in p for p in parts):
            is_sinogram = True
        elif "PROJECTION" in parts or "SINOGRAM" in parts:
            is_sinogram = True

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
        if all_tags:
            keys = sorted(header_reader.GetMetaDataKeys())
            for key in keys:
                val = clean_tag(header_reader.GetMetaData(key))
                tag_display = f"({key.replace('|', ',')})"
                print(f"{tag_display}: {val if val is not None else '(not set)'}")
            return

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
        print(f"Is Sinogram: {'yes' if is_sinogram else 'no'}")
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
    if is_sinogram:
        meta_dict["sinogram"] = True

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


def list_dicom_series(folder_path: Path) -> None:
    if not folder_path.exists():
        raise FileNotFoundError(f"Folder not found: {folder_path}")
    if not folder_path.is_dir():
        raise ValueError(f"Path is not a directory: {folder_path}")

    from tqdm import tqdm

    # Collect subdirectories recursively and include the root folder itself
    subdirs = [p for p in folder_path.rglob("*") if p.is_dir()]
    subdirs.append(folder_path)
    subdirs = sorted(list(set(subdirs)))

    print(f"Scanning folder recursively for DICOM series: {folder_path}\n")
    series_list = []

    # Disable SimpleITK global warnings during the scan to avoid "No Series were found" warnings
    was_warning_on = sitk.ProcessObject.GetGlobalWarningDisplay()
    sitk.ProcessObject.SetGlobalWarningDisplay(False)

    try:
        # Walk recursively with tqdm progress bar
        for path in tqdm(subdirs, desc="Scanning folders", unit="dir"):
            try:
                series_ids = sitk.ImageSeriesReader.GetGDCMSeriesIDs(str(path))
            except Exception:
                continue
            if not series_ids:
                continue
            for sid in series_ids:
                file_names = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(str(path), sid)
                if file_names:
                    first_file = file_names[0]
                    try:
                        reader = sitk.ImageFileReader()
                        reader.SetFileName(first_file)
                        reader.ReadImageInformation()
                        
                        def get_tag(tag_str: str) -> Optional[str]:
                            if reader.HasMetaDataKey(tag_str):
                                return clean_tag(reader.GetMetaData(tag_str))
                            return None

                        modality = get_tag("0008|0060") or "Unknown"
                        patient_id = get_tag("0010|0020") or "Unknown"
                        series_desc = get_tag("0008|103e") or get_tag("0008|1030") or ""
                        
                        study_date = get_tag("0008|0020")
                        study_time = get_tag("0008|0030")
                        series_date = get_tag("0008|0021")
                        series_time = get_tag("0008|0031")
                        acq_date = get_tag("0008|0022")
                        acq_time = get_tag("0008|0032")
                        acq_datetime = get_tag("0008|002a")

                        study_uid = get_tag("0020|000d") or "UnknownStudy"
                        frame_of_ref = get_tag("0020|0052")

                        inst_creation_date = get_tag("0008|0012")
                        inst_creation_time = get_tag("0008|0013")

                        best_dt = None
                        if acq_datetime:
                            best_dt = parse_datetime(acq_datetime)
                        if not best_dt:
                            best_dt = combine_date_time(acq_date, acq_time)
                        if not best_dt:
                            best_dt = combine_date_time(series_date, series_time)
                        if not best_dt:
                            best_dt = combine_date_time(study_date, study_time)
                        
                        if not best_dt:
                            best_dt = "(no date)"

                        recon_dt = combine_date_time(inst_creation_date, inst_creation_time) or best_dt

                        dim = reader.GetDimension()
                        size = reader.GetSize()
                        if dim >= 3 and len(file_names) == 1:
                            slices = size[2]
                        else:
                            slices = len(file_names)

                        # Clean/standardize modality name
                        m_lower = modality.lower()
                        if m_lower == "pt":
                            modality_str = "PET"
                        elif m_lower == "nm":
                            modality_str = "SPECT"
                        else:
                            modality_str = modality.upper()

                        # Check if it's a sinogram
                        is_sinogram = False
                        image_type_raw = get_tag("0008|0008")
                        if image_type_raw:
                            parts = [p.strip().upper() for p in image_type_raw.replace("\\", "/").replace(",", "/").split("/")]
                            if "TOMO" in parts and not any("RECON" in p for p in parts):
                                    is_sinogram = True
                            elif "PROJECTION" in parts or "SINOGRAM" in parts:
                                is_sinogram = True
                        
                        if is_sinogram:
                            modality_str += " (Sino)"

                        try:
                            rel_folder = path.relative_to(folder_path)
                            if str(rel_folder) == ".":
                                folder_str = "./"
                            else:
                                folder_str = f"./{rel_folder}"
                        except ValueError:
                            folder_str = str(path)

                        series_list.append({
                            "datetime": best_dt,
                            "recon_datetime": recon_dt,
                            "modality": modality_str,
                            "patient_id": patient_id,
                            "slices": slices,
                            "description": series_desc,
                            "folder": folder_str,
                            "uid": sid,
                            "study_uid": study_uid,
                            "frame_of_ref": frame_of_ref
                        })
                    except Exception:
                        continue
    finally:
        sitk.ProcessObject.SetGlobalWarningDisplay(was_warning_on)

    if not series_list:
        print("No DICOM series found.")
        return

    # Identify groups by (study_uid, frame_of_ref, date)
    from collections import defaultdict
    group_counts = defaultdict(int)
    for s in series_list:
        study_uid = s["study_uid"]
        ref = s["frame_of_ref"]
        day = s["datetime"][:10] if len(s["datetime"]) >= 10 else s["datetime"]
        if ref:
            group_counts[(study_uid, ref, day)] += 1

    groups = defaultdict(list)
    single_series = []
    for s in series_list:
        study_uid = s["study_uid"]
        ref = s["frame_of_ref"]
        day = s["datetime"][:10] if len(s["datetime"]) >= 10 else s["datetime"]
        if ref and group_counts[(study_uid, ref, day)] > 1:
            groups[(study_uid, ref, day)].append(s)
        else:
            single_series.append(s)

    blocks = []
    for key, group in groups.items():
        # Sort internally within the group by reconstruction date/time, then description/modality
        group.sort(key=lambda x: (x["recon_datetime"], x["datetime"], x["description"], x["modality"]))
        # Only indent subsequent series in the group (the first one is the group head and not indented)
        group[0]["indented"] = False
        for s in group[1:]:
            s["indented"] = True
        # Block datetime is the datetime of the first series (group head)
        block_datetime = group[0]["datetime"]
        blocks.append((block_datetime, group))

    for s in single_series:
        s["indented"] = False
        blocks.append((s["datetime"], [s]))

    # Sort blocks by representative datetime
    blocks.sort(key=lambda x: x[0])

    # Flatten back to list
    sorted_series_list = []
    for _, block in blocks:
        sorted_series_list.extend(block)

    headers = ["Date/Time", "Modality", "Patient ID", "Slices", "Description", "Folder"]
    rows = []
    for s in sorted_series_list:
        dt_str = s["datetime"]
        if s["indented"]:
            dt_str = f"  {dt_str}"
        rows.append([
            dt_str,
            s["modality"],
            s["patient_id"],
            str(s["slices"]),
            s["description"],
            s["folder"]
        ])

    # Compute column widths
    col_widths = [len(h) for h in headers]
    for row in rows:
        for i, val in enumerate(row):
            col_widths[i] = max(col_widths[i], len(val))

    def format_row(row_data):
        parts = []
        for i, val in enumerate(row_data):
            width = col_widths[i]
            if i in (0, 1, 2, 4, 5):  # left aligned
                parts.append(val.ljust(width))
            else:  # right aligned (Slices)
                parts.append(val.rjust(width))
        return " | ".join(parts)

    header_line = format_row(headers)
    separator = "-+-".join("-" * w for w in col_widths)

    print(header_line)
    print(separator)
    for row in rows:
        print(format_row(row))
    print()
