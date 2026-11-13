import json
from pathlib import Path
from typing import Optional
import logging

logger = logging.getLogger("radot")
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
    # A date without time is not a time point (no midnight default)
    return None

def combine_date_time(date_str: Optional[str], time_str: Optional[str]) -> Optional[str]:
    """'YYYY-MM-DD HH:MM:SS' from DICOM date and time; None when the time is missing (no midnight
    default)."""
    if not date_str:
        return None
    d = parse_date(date_str)
    if not d:
        return None
    t = parse_time(time_str) if time_str else None
    if not t:
        return None
    return f"{d} {t}"

# Element names and symbols of the radionuclides used in nuclear medicine
ELEMENT_SYMBOLS = {
    "actinium": "Ac", "astatine": "At", "bismuth": "Bi", "carbon": "C", "copper": "Cu",
    "fluorine": "F", "gallium": "Ga", "holmium": "Ho", "indium": "In", "iodine": "I",
    "lead": "Pb", "lutetium": "Lu", "nitrogen": "N", "oxygen": "O", "radium": "Ra",
    "rubidium": "Rb", "samarium": "Sm", "technetium": "Tc", "terbium": "Tb",
    "thorium": "Th", "yttrium": "Y", "zirconium": "Zr",
}

def parse_radionuclide(text: Optional[str]) -> Optional[str]:
    """
    Parse a DICOM radionuclide description into the radot name (e.g. 'Lu177', 'Tc99m').
    Accepts code meanings such as '177 Lutetium', '^177^Lutetium', 'Fluorine^18^[^18^Fluorodeoxyglucose]',
    and symbols such as 'Lu-177' or '177Lu'. Returns None if no radionuclide is recognized.
    """
    if not text:
        return None
    import re
    tokens = re.findall(r"[A-Za-z]+|\d+m?", text)
    symbols = set(ELEMENT_SYMBOLS.values())
    for i, tok in enumerate(tokens):
        symbol = ELEMENT_SYMBOLS.get(tok.lower())
        if symbol is None and tok in symbols:
            symbol = tok
        if symbol is None:
            continue
        # Mass number just after or just before the element
        for j in (i + 1, i - 1):
            if 0 <= j < len(tokens) and tokens[j][0].isdigit():
                return f"{symbol}{tokens[j]}"
    return None

def read_radiopharmaceutical(dicom_file: str) -> dict:
    """
    Read the Radiopharmaceutical Information Sequence (0054,0016) of a PET/SPECT DICOM file:
    injected activity (Bq), injection date/time and radionuclide. These tags are inside the
    sequence, where the SimpleITK header reader does not see them.
    """
    import pydicom

    try:
        ds = pydicom.dcmread(dicom_file, stop_before_pixels=True)
    except Exception as e:
        logger.debug(f"Could not read DICOM header with pydicom: {e}")
        return {}
    seq = ds.get("RadiopharmaceuticalInformationSequence")
    if not seq:
        return {}
    item = seq[0]
    info = {}
    for key in ("RadionuclideTotalDose", "RadiopharmaceuticalStartDateTime", "RadiopharmaceuticalStartTime"):
        val = item.get(key)
        if val not in (None, ""):
            info[key] = str(val)
    nuclide = None
    code_seq = item.get("RadionuclideCodeSequence")
    if code_seq:
        nuclide = parse_radionuclide(str(code_seq[0].get("CodeMeaning", "")))
    if nuclide is None:
        nuclide = parse_radionuclide(str(item.get("Radiopharmaceutical", "")))
    if nuclide is not None:
        info["radionuclide"] = nuclide
    return info

def select_series(dicom_path: Path, series: Optional[str] = None) -> str:
    """Series UID to read in a folder: the only one, or the one matching `series` (a unique part
    of the Series Instance UID). Several series without `series` is an error listing them."""
    ids = list(sitk.ImageSeriesReader.GetGDCMSeriesIDs(str(dicom_path)))
    if not ids:
        raise ValueError(f"No DICOM series found in directory: {dicom_path}")
    if series is not None:
        match = [sid for sid in ids if series in sid]
        if len(match) == 1:
            return match[0]
        what = "matches no series" if not match else f"matches {len(match)} series"
        raise ValueError(f"--series '{series}' {what} in {dicom_path}.\n" + _describe_series(dicom_path, ids))
    if len(ids) == 1:
        return ids[0]
    raise ValueError(
        f"{dicom_path} contains {len(ids)} DICOM series: choose one with --series "
        f"(a unique part of its UID).\n" + _describe_series(dicom_path, ids)
    )


def _describe_series(dicom_path: Path, ids) -> str:
    lines = []
    for sid in ids:
        files = sitk.ImageSeriesReader.GetGDCMSeriesFileNames(str(dicom_path), sid)
        desc = modality = ""
        try:
            r = sitk.ImageFileReader()
            r.SetFileName(files[0])
            r.ReadImageInformation()
            modality = clean_tag(r.GetMetaData("0008|0060")) if r.HasMetaDataKey("0008|0060") else ""
            desc = clean_tag(r.GetMetaData("0008|103e")) if r.HasMetaDataKey("0008|103e") else ""
        except Exception:
            pass
        lines.append(f"  {sid}  {modality or '?'}  {len(files)} file(s)  {desc or ''}".rstrip())
    return "\n".join(lines)


def convert_dicom(
    dicom_path: Path,
    output_path: Optional[Path] = None,
    all_tags: bool = False,
    series: Optional[str] = None,
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
        logger.info(f"Scanning directory for DICOM series: {dicom_path}")
        reader = sitk.ImageSeriesReader()
        series_uid = select_series(dicom_path, series)
        dicom_names = reader.GetGDCMSeriesFileNames(str(dicom_path), series_uid)
        if not dicom_names:
            raise ValueError(f"No DICOM series found in directory: {dicom_path}")
        header_file = dicom_names[0]
    else:
        if series is not None:
            raise ValueError("--series selects a series in a folder, not in a single DICOM file.")
        logger.info(f"Reading single DICOM file: {dicom_path}")
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
            logger.warning(f"Patient's Weight (0010,1030) '{weight_raw}' is not a number: weight not written.")

    # Radiopharmaceutical tags are in the sequence (0054,0016); top-level tags are a fallback
    radio = read_radiopharmaceutical(header_file)

    # Parse height (Patient's Size 0010|1020 is in m, we store in cm)
    height_raw = get_tag("0010|1020")
    height = None
    if height_raw:
        try:
            height = round(float(height_raw) * 100.0, 2)
        except ValueError:
            logger.warning(f"Patient's Size (0010,1020) '{height_raw}' is not a number: height not written.")
        if height is not None and height <= 0:
            height = None

    # Parse sex (Patient's Sex 0010|0040: M, F or O)
    sex_raw = get_tag("0010|0040")
    sex = {"M": "male", "F": "female"}.get(sex_raw.upper()) if sex_raw else None

    # Decay correction reference (0054|1102: NONE, START or ADMIN)
    decay_raw = get_tag("0054|1102")
    decay_correction = decay_raw.lower() if decay_raw else None

    radionuclide = radio.get("radionuclide")

    # Parse activity (Radionuclide Total Dose 0018|1074 is in Bq, we store in MBq)
    activity_raw = radio.get("RadionuclideTotalDose") or get_tag("0018|1074")
    activity = None
    if activity_raw:
        try:
            activity = float(activity_raw) / 1_000_000.0
        except ValueError:
            logger.warning(f"Radionuclide Total Dose (0018,1074) '{activity_raw}' is not a number: activity not written.")

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

    # CT values are HU by definition. PET/SPECT without Units (0054,1001): not guessed (it may be
    # counts); the unit stays unset and the tools that need it raise an explicit error
    if unit is None:
        if modality == "ct":
            unit = "hu"
        elif modality in ("pet", "spect") and output_path is not None:
            logger.warning(
                "No Units tag (0054,1001): the unit is not written to the sidecar. Set it with "
                "radot meta IMAGE --unit bqml|bq|suv|..."
            )

    # Study/Series/Acquisition Dates and Times
    study_date = get_tag("0008|0020")
    study_time = get_tag("0008|0030")
    series_date = get_tag("0008|0021")
    series_time = get_tag("0008|0031")
    acq_date = get_tag("0008|0022")
    acq_time = get_tag("0008|0032")
    acq_datetime = get_tag("0008|002a")

    # Resolve acquisition_date (date and time from the same tag pair, never a default time)
    acquisition_date = None
    acquisition_source = None
    for source, value in (
        ("Acquisition DateTime", parse_datetime(acq_datetime) if acq_datetime else None),
        ("Acquisition Date/Time", combine_date_time(acq_date, acq_time)),
        ("Series Date/Time", combine_date_time(series_date, series_time)),
        ("Study Date/Time", combine_date_time(study_date, study_time)),
    ):
        if value:
            acquisition_date, acquisition_source = value, source
            break
    if acquisition_source not in (None, "Acquisition DateTime", "Acquisition Date/Time"):
        logger.warning(f"acquisition_date from the {acquisition_source} tags (no acquisition date and time).")

    # Radiopharmaceutical details for injection_date
    rad_datetime = radio.get("RadiopharmaceuticalStartDateTime") or get_tag("0018|1078")
    rad_time = radio.get("RadiopharmaceuticalStartTime") or get_tag("0018|1072")
    
    injection_date = None
    if rad_datetime:
        injection_date = parse_datetime(rad_datetime)
    if not injection_date and rad_time and output_path is not None:
        # The injection date is not guessed from the acquisition date: a late time point is
        # acquired days after the injection
        t = parse_time(rad_time)
        logger.warning(
            f"Injection time {t} without date (no Radiopharmaceutical Start DateTime 0018,1078): "
            f"injection_date is not written. Set it with radot meta IMAGE --injection-date "
            f"'YYYY-MM-DD {t}'."
        )

    if output_path is None:
        if all_tags:
            keys = sorted(header_reader.GetMetaDataKeys())
            import sys
            for key in keys:
                val = clean_tag(header_reader.GetMetaData(key))
                tag_display = f"({key.replace('|', ',')})"
                sys.stdout.write(f"{tag_display}: {val if val is not None else '(not set)'}\n")
            return

        tags_to_print = [
            ("0010|0010", "Patient's Name"),
            ("0010|0020", "Patient ID"),
            ("0010|1030", "Patient's Weight"),
            ("0010|1020", "Patient's Size"),
            ("0010|0040", "Patient's Sex"),
            ("0054|1102", "Decay Correction"),
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
        import sys
        # Tags read from the Radiopharmaceutical Information Sequence
        seq_values = {
            "0018|1074": radio.get("RadionuclideTotalDose"),
            "0018|1078": radio.get("RadiopharmaceuticalStartDateTime"),
            "0018|1072": radio.get("RadiopharmaceuticalStartTime"),
        }
        for key, name in tags_to_print:
            val = seq_values.get(key) or get_tag(key)
            tag_display = f"({key.replace('|', ',')})"
            sys.stdout.write(f"{tag_display} {name}: {val if val is not None else '(not set)'}\n")
        sys.stdout.write(f"Radionuclide: {radionuclide if radionuclide else '(not set)'}\n")
        sys.stdout.write(f"Is Sinogram: {'yes' if is_sinogram else 'no'}\n")
        return

    # Load full image only if saving to path
    if is_dir:
        reader.SetFileNames(dicom_names)
        reader.MetaDataDictionaryArrayUpdateOn()
        logger.info(f"Reading DICOM series ({len(dicom_names)} slices)...")
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
    if height is not None:
        meta_dict["height"] = height
    if sex:
        meta_dict["sex"] = sex
    if radionuclide:
        meta_dict["radionuclide"] = radionuclide
    if decay_correction:
        meta_dict["decay_correction"] = decay_correction
    if is_sinogram:
        meta_dict["sinogram"] = True

    # Write output image
    logger.info(f"Saving converted image to: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(image, str(output_path))

    # Write output sidecar JSON
    from radot.meta import get_sidecar_path
    sidecar_path = get_sidecar_path(output_path)
    logger.info(f"Saving metadata sidecar to: {sidecar_path}")
    with open(sidecar_path, "w") as f:
        json.dump(meta_dict, f, indent=4)

    logger.info("DICOM conversion completed successfully.")


def list_dicom_series(folder_path: Path) -> list:
    if not folder_path.exists():
        raise FileNotFoundError(f"Folder not found: {folder_path}")
    if not folder_path.is_dir():
        raise ValueError(f"Path is not a directory: {folder_path}")

    from tqdm import tqdm

    # Collect subdirectories recursively and include the root folder itself
    subdirs = [p for p in folder_path.rglob("*") if p.is_dir()]
    subdirs.append(folder_path)
    subdirs = sorted(list(set(subdirs)))
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
                            # date only (display: no time shown, none invented)
                            d = next((x for x in (acq_date, series_date, study_date) if x), None)
                            best_dt = (parse_date(d) if d else None) or "(no date)"

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
        logger.warning("No DICOM series found.")
        return []

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

    return sorted_series_list


def print_dicom_series(series_list: list) -> None:
    import sys
    if not series_list:
        sys.stdout.write("No DICOM series found.\n")
        return

    headers = ["Date/Time", "Modality", "Patient ID", "Slices", "Description", "Folder"]
    rows = []
    for s in series_list:
        dt_str = s["datetime"]
        if s.get("indented", False):
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

    sys.stdout.write(header_line + "\n")
    sys.stdout.write(separator + "\n")
    for row in rows:
        sys.stdout.write(format_row(row) + "\n")
    sys.stdout.write("\n")
