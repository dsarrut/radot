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
                      The only date format of radot (sidecars, options, configs; an optional
                      fraction of second is accepted). Any other format (day first, slashes,
                      date without time) is an error.
  - time_from_injection: Time elapsed from injection to acquisition (computed helper option, not stored directly).
                         Format: day:hours:min (e.g. 1:12:30)
                         The tools use the dates. A `time_from_injection` key in hours (written by other
                         software, never by radot) is used only without dates; with dates and a
                         difference above one minute, it is an error.
  - height: Patient height (in cm).
            Format: floating point number
  - sex: Patient sex.
         Valid choices: male, female
  - radionuclide: Radionuclide (option --radionuclide or --rad).
                  Format: name known to ICRP 107 (e.g. Lu177, Y90, Tc99m)
  - decay_correction: Time reference of the activity values (DICOM 0054,1102).
                      Valid choices: start, none, admin (decay-corrected to the injection)
                      Required for activity images (bq, bqml, suv, sul) used at a time point.
```

---

### Writing Metadata: `radot meta` with options
Applies metadata edits to the image JSON sidecar file. Only the specified options are modified; unspecified fields are merged and preserved. A sidecar that is not valid JSON (or not a JSON object) is an error, for writing and for reading: it is never rewritten from scratch nor ignored; fix it, or delete it to start again.

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
* `--height <float>`: Patient height in cm.
* `--sex <male|female>`: Patient sex.
* `--radionuclide <name>` (or `--rad`): Radionuclide, e.g. `Lu177`; must be known to ICRP 107.
* `--decay-correction <start|none|admin>`: Time reference of the activity values (DICOM Decay Correction 0054,1102). See [decay correction](#decay-correction-of-the-activity-images) below.

#### Decay correction of the activity images

The reconstruction may correct the activity for the physical decay to a reference time. The sidecar key `decay_correction` says which:

| `decay_correction` | The image holds | radot uses |
|--------------------|-----------------|------------|
| `start` | the activity at the start of the acquisition | the values as they are |
| `none` | the activity during the acquisition, not corrected | the values as they are (the decay during the acquisition is neglected) |
| `admin` | the activity decay-corrected to the injection time, $A(t)\,e^{\lambda t}$ | the values $\times\, e^{-\lambda t}$: the activity at the acquisition time |

There is no default: an `admin` image read as not corrected gives an activity, and a dose, too high by $e^{\lambda t}$ (for ¹⁷⁷Lu about +10 % at 24 h and +50 % at 96 h). The key is required wherever the activity at a time point is used: `radot dose` on an activity image (`bq`, `bqml`, `suv`, `sul`), multiple time points, and the conversions to and from SUV/SUL. A TIA image (`bqs`, `bqsml`), a dose, a CT or a label image does not need it. `radot dicom` writes it from the DICOM tag when present; it is usually present for PET and may be absent for SPECT (a WARNING says so): set it with `radot meta IMAGE --decay-correction start`. A unit conversion keeps it (the values are scaled, their time reference does not change). For an in-memory image with the Python API: `image.SetMetaData("decay_correction", "start")`.

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
