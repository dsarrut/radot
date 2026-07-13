# Dosimetry Module Documentation

The radot dosimetry module computes absorbed doses (in Gy) or dose rates (in Gy/h) using organ-level S-value methods or a voxel-level Local Deposition Model (LDM).

---

## 1. CLI Commands

### Listing Dosimetry Models: `radot dose-list`
Prints equations and descriptions of all available dosimetry methods.

```bash
radot dose-list
```

---

### Computing Dose: `radot dose`
Calculates absorbed doses (or dose rates) and either writes a 3D dose NIfTI map or prints a formatted summary table to console.

```bash
radot dose <image_path> [options]
```

#### Options:
* -l, --labels <labels_path>: Optional anatomical segmentation label map. Required for S-value methods. Optional for LDM.
* --ct <ct_path>: Companion CT image used for physical density estimation and mass calibration.
* -m, --method <hanscheid2017|hanscheid2018|ldm>: Dosimetry model to apply.
* --isotope <name>: Radionuclide identifier (e.g. Lu177, Y90) to retrieve physical yield constraints.
* --rois <name>: ROI name to include. Can be specified multiple times or comma-separated.
* --teff <ROI:value>: Effective half-life override per ROI in hours (e.g. --teff Liver:45.0 --teff Kidney_L:32.0).
* --total-activity <float>: Global calibration activity override in MBq.
* --hu-threshold <float>: Hounsfield Units cutoff (defaults to -800.0). Voxel doses corresponding to CT values below this limit are zeroed out (air masking).
* --dose-rate: If set, computes instantaneous dose rate (in Gy/h) at the time of acquisition rather than integrated absorbed dose (in Gy). Only available for LDM.
* -o, --output <output_path>: Path to output NIfTI dose map or JSON summary results. Required for LDM if running without labels.
* --resample <like_activity|like_labels|like_ct>: Specifies the grid alignment strategy if input image spaces differ.
* --phantom <male|female>: Reference phantom model used for S-value calculations. Defaults to male.

#### Examples:

1. Organ-level S-value dosimetry (Hänscheid 2017):
```bash
radot dose spect.nii.gz -l labels.nii.gz -m hanscheid2017 --isotope Lu177 --ct ct.nii.gz --teff Liver:45.0 --teff Kidney_L:32.0 --rois Liver --rois Kidney_L
```

2. Voxel-level LDM absorbed dose map:
```bash
radot dose spect.nii.gz -m ldm --isotope Lu177 --ct ct.nii.gz -o dose.nii.gz
```

3. Voxel-level LDM dose rate map:
```bash
radot dose spect.nii.gz -m ldm --isotope Lu177 --ct ct.nii.gz --dose-rate -o dose_rate.nii.gz
```

---

## 2. Python API

Manage dosimetry programmatically:

### `compute_roi_dose_file`
Executes the full dosimetry pipeline on file paths.

```python
from pathlib import Path
from radot.dose import compute_roi_dose_file

compute_roi_dose_file(
    image_path=Path("spect.nii.gz"),
    labels_path=Path("labels.nii.gz"),
    method="ldm",
    isotope="Lu177",
    ct_path=Path("ct.nii.gz"),
    output_path=Path("dose.nii.gz"),
    hu_threshold=-800.0,
    dose_rate=False
)
```

---

### `compute_roi_dose`
Executes dosimetry on in-memory SimpleITK image objects.

```python
import SimpleITK as sitk
from radot.dose import compute_roi_dose

spect = sitk.ReadImage("spect.nii.gz")
ct = sitk.ReadImage("ct.nii.gz")
labels = sitk.ReadImage("labels.nii.gz")

# Returns a dictionary of ROI statistics and the 3D SimpleITK dose image
results, dose_image = compute_roi_dose(
    image=spect,
    method="ldm",
    isotope="Lu177",
    labels=labels,
    ct_image=ct,
    hu_threshold=-800.0,
    dose_rate=True,
    labels_map={1: "Liver", 2: "Kidney_L"}
)
```

---

### `dose_ldm`
Low-level physical dosimetry function representing the core vectorized mathematical model of LDM.

```python
import numpy as np
from radot.dose import dose_ldm

# Inputs are numpy arrays
dose_array = dose_ldm(
    image_arr=np.ones((64, 64, 64)),      # Activity array (in Bq)
    voxel_volume_ml=0.008,                # Voxel volume (in mL)
    time_from_injection_h=24.0,           # Timing parameter
    isotope="Lu177",
    ct_arr=np.zeros((64, 64, 64)),        # Companion CT HUs
    hu_threshold=-800.0,
    dose_rate=False
)
```
