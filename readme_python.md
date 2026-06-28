# Python API Reference Guide

This document details the programmatic Python interface for the radot package, structured around the main dosimetry and image processing workflow.

---

## 1. DICOM Conversion

Functions to list and convert DICOM datasets to dosimetry-ready formats:

```python
from pathlib import Path
from radot import list_dicom_series, print_dicom_series, convert_dicom

folder_path = Path("dicom_directory")

# Scan directory for DICOM series
series = list_dicom_series(folder_path)

# Display detected series as a formatted table
print_dicom_series(series)

# Convert a DICOM CT directory to NIfTI format
convert_dicom(
    dicom_path=folder_path / "ct_series",
    output_path=Path("ct.nii.gz")
)

# Convert a single SPECT DICOM file
convert_dicom(
    dicom_path=folder_path / "spect.dcm",
    output_path=Path("spect.nii.gz")
)
```

---

## 2. Image Metadata Management

Functions to read and write dosimetry metadata in JSON sidecars:

```python
from pathlib import Path
from radot import read_metadata, write_metadata

image_path = Path("spect.nii.gz")

# Populate companion JSON sidecar with dosimetry parameters
write_metadata(
    image_path=image_path,
    modality="spect",
    unit="bq",
    weight=98.0,
    activity=7400.0,
    injection_date="2021-10-20 09:00:00",
    acquisition_date="2021-10-21 12:10:09"
)

# Load metadata dictionary
meta = read_metadata(image_path)
print(f"Modality: {meta['metadata']['modality']}")
print(f"Time from injection: {meta['metadata']['time_from_injection']} hours")
```

---

## 3. Segmentation and Remapping

Segments CT scans and maps label nomenclature to AAPM TG-263 guidelines:

```python
from pathlib import Path
from radot.segmentation import segment_ct

# Segment CT scan, align grid, and write TG-263 JSON mapping sidecar
segment_ct(
    ct_path=Path("ct.nii.gz"),
    output_path=Path("labels_all.nii.gz"),
    mode="total",
    fast=True
)
```

---

## 4. Image Resampling

Aligns mismatched grids or resamples to isotropic voxel spacings:

```python
from pathlib import Path
from radot.resampling import resample_image_file, resample_image_file_like

# Resample CT to 4mm isotropic grid
resample_image_file(
    image_path=Path("ct.nii.gz"),
    spacing=(4.0, 4.0, 4.0),
    image_type="ct",
    output_path=Path("ct_4mm.nii.gz")
)

# Align SPECT grid to match the label map layout (preserving total activity)
resample_image_file_like(
    image_path=Path("spect.nii.gz"),
    like_image_path=Path("labels_all.nii.gz"),
    type="spect",
    unit="bq",
    preserve_total=True,
    output_path=Path("spect_aligned.nii.gz")
)
```

---

## 5. Regional Statistics

Computes physical volumes, tissue masses, and activity metrics:

```python
from pathlib import Path
from radot.stats import compute_stats_file

# Compute regional statistics and print output table
compute_stats_file(
    labels_path=Path("labels_all.nii.gz"),
    image_path=Path("spect.nii.gz"),
    image_type="spect",
    unit="bq"
)
```

---

## 6. Dosimetry Calculation

Calculates organ-absorbed doses from activity distribution and decay constraints:

```python
from radot.dose import compute_roi_dose_file

# Compute dosimetry using Hanscheid 2017 formulation
results = compute_roi_dose_file(
    image_path="spect.nii.gz",
    labels_path="labels_all.nii.gz",
    method="hanscheid2017",
    rois=["Liver", "Kidney_L"],
    isotope="Lu177",
    ct_path="ct.nii.gz",
    resample_target="like_activity",
    method_params={
        "Liver": {"teff": 45.0},
        "Kidney_L": {"teff": 32.0}
    },
    verbose=True
)

for roi, data in results.items():
    print(f"{roi}: {data['dose_gy']:.4f} Gy (Mass: {data['mass_g']:.2f} g)")
```

---

## 7. Radionuclide Decay Lookup

Accesses physical half-life properties from the ICRP-107 database:

```python
from radot.db import get_physical_half_life_h

# Retrieve physical half-life in hours
tp_hours = get_physical_half_life_h("Lu177")
print(f"Lu-177 Half-Life: {tp_hours:.2f} hours")
```

---

## 8. Logger Level Control

Configures package-level verbosity to silence progress output or show debugging logs:

```python
import radot

# Silence all progress and warning logs (only show errors)
radot.set_log_level("QUIET")

# Restore standard level
radot.set_log_level("INFO")
```
