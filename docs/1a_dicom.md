# DICOM Module Documentation

The `radot` DICOM module provides tools to scan directories for DICOM series, inspect headers, and convert them to therapy dosimetry-ready continuous formats (e.g. NIfTI, NRRD) with automatic metadata extraction.

---

## 1. CLI Commands

### Scanning Directories: `radot dicom-list`
Scans a folder and all its subfolders recursively to list all available DICOM series.

```bash
radot dicom-list <folder_path>
```

#### Features:
- Recursive Scanning: Automatically discovers all series in deeply nested subdirectories.
- Visual Grouping: Series belonging to the same Study and Frame of Reference are automatically grouped, sorted chronologically, and indented for easy correlation.
- Modality Identification: Detects and standardizes CT, PET, SPECT (derived from Nuclear Medicine `NM`), and sinogram projection datasets.

#### Example Output:
```text
Scanning folder recursively for DICOM series: /data/patient_01

Date/Time           | Modality     | Patient ID | Slices | Description              | Folder    
--------------------+--------------+------------+--------+--------------------------+-----------
2026-06-25 08:30:00 | CT           | PAT001     |    120 | CT Thorax                | ./ct      
2026-06-25 09:00:00 | SPECT        | PAT001     |     64 | NM SPECT Activity        | ./spect   
  2026-06-25 09:15:00 | SPECT (Sino) | PAT001     |    120 | NM Sinogram Projection   | ./spect/raw
```

---

### Converting DICOM: `radot dicom`
Converts a single DICOM file or an entire DICOM series directory into standard continuous formats (e.g., `.nii.gz`, `.nrrd`, `.mhd`).

```bash
radot dicom <input_path> -o <output_path> [--all-tags]
```

#### Options:
* `-o, --output <output_path>`: The output image file path. 
  > [!WARNING]
  > Writing back to DICOM format is not supported. Use formats like NIfTI (`.nii.gz`), NRRD (`.nrrd`), or MetaImage (`.mhd`).
* `--all-tags`: Print all DICOM tags present in the header (only works when the `-o` / `--output` option is omitted).

#### Behavior:
1. Header Inspection: Inspects the primary header file.
2. Metadata Extraction: Extracts key metadata:
   - Modality (e.g. `ct`, `pet`, `spect`)
   - Unit from the DICOM Units tag (e.g. `hu` for CT; `bqml` or `suv` for PET/SPECT). Without a Units tag, PET/SPECT images are assumed to be in `bqml`.
   - Patient's weight (in kg), height (Patient's Size, converted from m to cm) and sex (`male` / `female`)
   - Radionuclide injected activity (converted from Bq to MBq)
   - Radionuclide (e.g. `Lu177`), parsed from the Radionuclide Code Sequence (e.g. "177 Lutetium") or the Radiopharmaceutical name
   - Acquisition and Injection date/time stamps
   - Decay correction reference (`none`, `start` or `admin`, DICOM tag Decay Correction)

   The injected activity, injection date/time and radionuclide are read from the Radiopharmaceutical Information Sequence (0054,0016), where PET and SPECT scanners store them; top-level tags are used as a fallback.
3. JSON Sidecar Generation: Writes the extracted metadata to a companion `.json` sidecar matching the output image filename (e.g., converting `spect.dcm` to `spect.nii.gz` creates a `spect.json` sidecar).

---

## 2. Python API

You can programmatically utilize the DICOM module directly in your scripts:

### `list_dicom_series` & `print_dicom_series`
- `list_dicom_series(folder_path)`: Gathers and returns a list of dictionaries containing metadata for all discovered DICOM series, sorted chronologically and grouped by study/frame of reference.
- `print_dicom_series(series_list)`: Takes the list returned by `list_dicom_series` and prints a beautifully aligned table to the console.

```python
from pathlib import Path
from radot import list_dicom_series, print_dicom_series

folder = Path("./data/patient_01")

# Gather series metadata silently
series = list_dicom_series(folder)

# Format and print the table to the console
print_dicom_series(series)
```

---

### `convert_dicom`
Converts DICOM files and writes their matching sidecar JSON files.

```python
from pathlib import Path
from radot.dicom import convert_dicom

# Convert a CT series directory to NIfTI
convert_dicom(
    dicom_path=Path("./data/patient_01/ct"),
    output_path=Path("./data/patient_01/ct.nii.gz")
)

# Convert a single SPECT DICOM file
convert_dicom(
    dicom_path=Path("./data/patient_01/spect.dcm"),
    output_path=Path("./data/patient_01/spect.nii.gz")
)
```

#### Sidecar JSON Metadata Fields:
When converted, a matching `.json` file is written containing fields such as:
```json
{
    "modality": "spect",
    "unit": "bqml",
    "weight": 70.0,
    "activity": 185.0,
    "injection_date": "2026-06-25 08:00:00",
    "acquisition_date": "2026-06-25 09:00:00",
    "height": 175.0,
    "sex": "male",
    "radionuclide": "Lu177",
    "decay_correction": "start"
}
```
These fields are automatically recognized by other `radot` modules (like `resample` and `dose`) to perform dosimetry calculation without requiring manual command options.
