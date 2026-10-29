# Metadata Module Documentation

The radot metadata module manages the reading, writing, and calculations of therapy-related attributes for medical images (CT, PET, SPECT). It relies on companion JSON sidecar files to store custom metadata, keeping the primary image headers clean.

---

## 1. CLI Commands

### Reading Metadata: `radot meta` or `radot info`
Prints detailed image grid attributes along with the therapy metadata.

```bash
radot meta <image_path>
radot info <image_path>
```

#### Example Output:
```text
Image Size: [128, 128, 64]
Spacing: [2.0, 2.0, 2.0]
Origin: [0.0, 0.0, 0.0]
Orientation Matrix: [1.0000, 0.0000, 0.0000]
                    [0.0000, 1.0000, 0.0000]
                    [0.0000, 0.0000, 1.0000]
Pixel Type: 32-bit float
Number of Voxels: 1048576
Memory Size: 4.00 MB
Modality: SPECT
Unit: bqml
Weight: 75.0 kg
Activity: 250.0 MBq
Injection Date: 2026-06-25 08:00:00
Acquisition Date: 2026-06-25 09:00:00
Time from injection: 1.00 hours
```

---

### Listing Supported Keys: `radot meta-list`
Lists all supported metadata keys, choices, formats, and descriptions.

```bash
radot meta-list
```

#### Example Output:
```text
Supported image metadata keys:
  - modality: Image modality type.
              Valid choices: ct, pet, spect, dose
  - unit: Physical unit of voxel values.
          Valid choices: hu, suv, bqml, bq, gy, gy/h, sul, bqs, bqsml
  - weight: Patient weight (in kg).
            Format: floating point number
  - activity: Injected activity (in MBq).
              Format: floating point number
  - injection_date: Injection date and time.
                    Format: YYYY-MM-DD HH:MM:SS
  - acquisition_date: Acquisition date and time.
                      Format: YYYY-MM-DD HH:MM:SS
  - time_from_injection: Time elapsed from injection to acquisition (computed helper option, not stored directly).
                         Format: day:hours:min (e.g. 1:12:30)
```

---

### Writing Metadata: `radot meta` with options
Applies metadata edits to the image JSON sidecar file. Only the specified options are modified; unspecified fields are merged and preserved.

```bash
radot meta <image_path> [options]
```

#### Options:
* `--modality <ct|pet|spect|dose>`: Modality of the image.
* `--unit <hu|suv|bqml|bq|gy|gy/h|sul|bqs|bqsml>`: Physical unit of voxel values (gy for dose map, gy/h for dose rate map, bqs / bqsml for a time-integrated activity map in Bq·s per voxel / Bq·s/mL).
* `--weight <float>`: Patient weight in kg.
* `--activity <float>`: Injected activity in MBq.
* `--injection-date <YYYY-MM-DD HH:MM:SS>`: Injection timestamp.
* `--acquisition-date <YYYY-MM-DD HH:MM:SS>`: Acquisition timestamp.
* `--time-from-injection <day:hours:min>`: Time elapsed between injection and acquisition (e.g. 1:12:00 for 1 day and 12 hours).

#### Default Unit Mapping:
When saving metadata, if modality is specified but the unit is left empty, it defaults automatically:
* modality ct defaults unit to hu
* modality dose defaults unit to gy

#### Automatic Time Calculations:
The --time-from-injection option simplifies time calculations:
- If only --injection-date is present: Acquisition date is computed and saved.
- If only --acquisition-date is present: Injection date is computed and saved.
- If neither date is present: Defaults to 2000-01-01 12:00:00 as injection date, and computes acquisition date.

---

## 2. Python API

You can programmatically manage metadata directly in your scripts:

### `read_metadata`, `write_metadata` & `print_metadata`
- `read_metadata(image_path)`: Reads and returns a dictionary containing metadata and image grid information from the JSON sidecar and file headers.
- `write_metadata(image_path, ...)`: Writes or merges new metadata values into the companion JSON sidecar file.
- `print_metadata(info_dict)`: Takes the dictionary returned by `read_metadata` and prints a formatted information table to the console.

```python
from pathlib import Path
from radot import read_metadata, write_metadata, print_metadata

image = Path("spect.nii.gz")

# Write metadata using the explicit write_metadata function
write_metadata(
    image_path=image,
    modality="spect",
    unit="bq",
    weight=75.0,
    activity=250.0,
    injection_date="2026-06-25 08:00:00",
    acquisition_date="2026-06-25 09:00:00"
)

# Read metadata dictionary silently
info = read_metadata(image)

# Format and print the table to the console
print_metadata(info)
```

---

### `get_metadata_case_insensitive`
Helper to fetch a metadata field case-insensitively, automatically checking the sidecar first, and falling back to image headers if not found.

```python
import SimpleITK as sitk
from pathlib import Path
from radot.meta import get_metadata_case_insensitive

image_path = Path("spect.nii.gz")
image = sitk.ReadImage(str(image_path))

# Fetch weight
weight = get_metadata_case_insensitive(
    image=image,
    key="weight",
    image_path=image_path
)
print(f"Loaded patient weight: {weight} kg")
```

---

### `get_sidecar_path`
Returns the expected companion JSON sidecar path for a given image.

```python
from pathlib import Path
from radot.meta import get_sidecar_path

image = Path("spect.nii.gz")
sidecar = get_sidecar_path(image)
print(sidecar) # Prints "spect.json"
```
