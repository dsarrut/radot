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

---

### `roi_dose_uncertainty_bounds` & `roi_dose_uncertainty_estimate`
Tools for estimating the Monte Carlo statistical uncertainty of an organ / ROI mean dose from 3D dose and voxel relative uncertainty maps.

```python
import SimpleITK as sitk
from radot.dose import roi_dose_uncertainty_bounds, roi_dose_uncertainty_estimate

dose_img = sitk.ReadImage("output_dose.mhd")
unc_img = sitk.ReadImage("output_dose_uncertainty.mhd")
mask_img = sitk.ReadImage("liver_mask.nii.gz")

# Theoretical lower bound (independent voxels) and upper bound (100% correlation)
unc_min, unc_max = roi_dose_uncertainty_bounds(dose_img, unc_img, mask_img)

# Realistic physical estimate accounting for particle range in soft tissue
unc_est = roi_dose_uncertainty_estimate(dose_img, unc_img, mask_img, isotope="Lu177")

print(f"ROI Uncertainty: {unc_est * 100:.2f}% (bounds: [{unc_min * 100:.3f}%, {unc_max * 100:.2f}%])")
```

---

## 3. Dosimetry Methods Summary

### Hanscheid 2017 (`hanscheid2017`)
The Hänscheid 2017 method is an organ-level dosimetry approach specific to Lutetium-177 ($\text{}^{177}\text{Lu}$) that scales the mean activity concentration by a scaling factor of $0.125$ and the effective half-life.

**Equation:**
$$
D = 0.125 \cdot C_{\text{mean}} \cdot 2^{\frac{t}{T_{\text{eff}}}} \cdot T_{\text{eff}}
$$

Where:
*   $D$ is the mean absorbed dose in the ROI (in Gy).
*   $C_{\text{mean}}$ is the mean activity concentration in the ROI at the time of acquisition (in $\text{MBq}/\text{mL}$).
*   $t$ is the time from injection to acquisition (in hours).
*   $T_{\text{eff}}$ is the effective half-life in the ROI (in hours). If not specified, falls back to the physical half-life of the isotope.

---

### Hanscheid 2018 (`hanscheid2018`)
The Hänscheid 2018 method is an organ-level dosimetry method that uses the OpenDose S-value database for phantom-specific regional dosimetry. It assumes that the effective half-life is $2 \cdot t$ ($T_{\text{eff}} = 2 \cdot t$), meaning that the time integrated activity (TIA) is computed from a single SPECT scan. S-values are scaled by the ratio of reference organ mass to patient organ mass.

**Equation:**
$$
D = \frac{\tilde{A} \cdot S}{\text{mass-scaling}} \cdot 10^{-3}
$$

Where:
*   $D$ is the mean absorbed dose in the ROI (in Gy).
*   $\tilde{A}$ is the Time Integrated Activity (TIA) in the ROI (in $\text{MBq}\cdot\text{s}$), calculated as:
    $$
    \tilde{A} = A(t) \cdot \frac{2 \cdot t \cdot 3600}{\ln(2)}
    $$
    with $A(t)$ being the total activity in the ROI (in MBq) at the time of acquisition, and $t$ being the time from injection (in hours).
*   $S$ is the self S-value of the organ from the OpenDose database (in $\text{mGy}/(\text{MBq}\cdot\text{s})$).
*   $\text{mass-scaling}$ is the mass scaling factor:
    $$
    \text{mass-scaling} = \frac{m_{\text{patient}}}{m_{\text{reference}}}
    $$
    where $m_{\text{patient}}$ is the patient's organ mass and $m_{\text{reference}}$ is the reference phantom's organ mass.
*   $10^{-3}$ converts the resulting dose from mGy to Gy.

#### S-Value Database Lookup & Missing Organs (e.g. Lesions)
1. **S-Value Resolution & Cache:**
   * S-value matrices (organized by phantom and isotope) are lazy-loaded from local JSON files under `radot/data/opendose/<phantom>/<isotope>.json`. If not present locally, they are automatically downloaded on-demand from the online Creatis registry.
   * Input ROI names are mapped to standard OpenDose reference organ names using `resolve_organ_name()` (which employs a multi-stage cascade matching against TG-263 standards).
   * Key matching in the database is flexible (case-insensitive and ignores minor suffixes like "contents", "wall", "lobe", etc.).

2. **Handling Missing Organs (e.g. Lesions):**
   * If an organ name cannot be mapped to the S-value database (e.g., custom regions like `lesion_1` or unmapped structures), the system raises a `ValueError`.
   * **In standard (generic) runs** (when looping over all label IDs in the volume): The system catches this error, logs a warning, and skips dosimetry calculations for that specific ROI.
   * **In explicit requests** (when specific `--rois` are requested via the CLI or parameters): The system raises a `ValueError` immediately, halting execution to inform you that a requested ROI does not exist in the reference database (suggesting adjusting sidecar mappings or using custom S-value overrides).

---

### Madsen 2018 (`madsen2018`)
The Madsen 2018 method is an organ-level dosimetry method that uses the OpenDose S-value database for phantom-specific regional dosimetry, similar to `hanscheid2018`, but allows/requires a user-specified effective half-life $T_{\text{eff}}$ per ROI instead of assuming $T_{\text{eff}} = 2 \cdot t$. S-values are scaled by the ratio of reference organ mass to patient organ mass.

**Equation:**
$$
D = \frac{\tilde{A} \cdot S}{\text{mass-scaling}} \cdot 10^{-3}
$$

Where:
*   $D$ is the mean absorbed dose in the ROI (in Gy).
*   $\tilde{A}$ is the Time Integrated Activity (TIA) in the ROI (in $\text{MBq}\cdot\text{s}$), calculated assuming mono-exponential clearance with the specified effective half-life $T_{\text{eff}}$:
    $$
    \tilde{A} = A(t) \cdot \frac{e^{k \cdot t \cdot 3600}}{k}
    $$
    with $A(t)$ being the total activity in the ROI (in MBq) at the time of acquisition, $t$ being the time from injection (in hours), and $k$ being the effective clearance rate:
    $$
    k = \frac{\ln(2)}{T_{\text{eff}} \cdot 3600}
    $$
    where $T_{\text{eff}}$ is the effective half-life in the ROI (in hours).
*   $S$ is the self S-value of the organ from the OpenDose database (in $\text{mGy}/(\text{MBq}\cdot\text{s})$).
*   $\text{mass-scaling}$ is the mass scaling factor:
    $$
    \text{mass-scaling} = \frac{m_{\text{patient}}}{m_{\text{reference}}}
    $$
    where $m_{\text{patient}}$ is the patient's organ mass and $m_{\text{reference}}$ is the reference phantom's organ mass.
*   $10^{-3}$ converts the resulting dose from mGy to Gy.


---

### Local Deposition Model (`ldm`)
The Local Deposition Model (LDM) is a voxel-level dosimetry method based on the physical deposition of non-penetrating radiations. It calculates the energy absorbed locally in each voxel.

#### Integrated Absorbed Dose (Gy)
**Equation:**
$$
D(x,y,z) = \frac{\tilde{A}(x,y,z) \cdot \Delta_{\text{local}} \cdot 1.60218 \times 10^{-10}}{m(x,y,z)}
$$

Where:
*   $D(x,y,z)$ is the absorbed dose in voxel $(x,y,z)$ (in Gy).
*   $\tilde{A}(x,y,z)$ is the Time Integrated Activity (TIA) in the voxel (in $\text{Bq}\cdot\text{s}$), calculated assuming physical decay:
    $$
    \tilde{A}(x,y,z) = A(t, x,y,z) \cdot \frac{e^{\lambda \cdot t}}{\lambda_{\text{seconds}}
    }
    $$
    where $A(t, x,y,z)$ is the activity in the voxel at acquisition time (in Bq), $t$ is the time from injection (in hours), and $\lambda$ is the physical decay constant ($\lambda = \ln(2)/T_p$ in $\text{hours}^{-1}$, $\lambda_{\text{seconds}} = \lambda / 3600$ in $\text{s}^{-1}$).
*   $\Delta_{\text{local}}$ is the non-penetrating energy constant for the isotope (in $\text{MeV}/\text{decay}$).
*   $m(x,y,z)$ is the voxel mass (in grams), computed from voxel volume ($V_{\text{voxel}}$ in mL) and density ($\rho$ in $\text{g}/\text{mL}$):
    $$
    m(x,y,z) = V_{\text{voxel}} \cdot \rho(x,y,z)
    $$
    where density $\rho(x,y,z)$ is estimated from the companion CT image in Hounsfield Units (HU):
    $$
    \rho(x,y,z) = 1.0 + \frac{\text{HU}(x,y,z)}{1000}
    $$
    (capped at a minimum of $0.001\text{ g/mL}$).
*   $1.60218 \times 10^{-10}$ is the conversion factor from $\text{MeV}/\text{g}$ to $\text{Gy}$ ($\text{J}/\text{kg}$):
    $$1\text{ MeV} = 1.60218 \times 10^{-13}\text{ J} \implies 1\text{ MeV/g} = 1.60218 \times 10^{-10}\text{ J/kg}$$

#### Instantaneous Dose Rate (Gy/h)
When `--dose-rate` is active, it computes the instantaneous dose rate at the acquisition time.

**Equation:**
$$
\dot{D}(x,y,z) = \frac{A(t, x,y,z) \cdot \Delta_{\text{local}} \cdot 1.60218 \times 10^{-10} \cdot 3600}{m(x,y,z)}
$$


Where:
*   $\dot{D}(x,y,z)$ is the dose rate in voxel $(x,y,z)$ (in $\text{Gy}/\text{h}$).
*   $A(t, x,y,z)$ is the activity in voxel $(x,y,z)$ at acquisition time (in Bq).
*   $3600$ converts seconds to hours.

## 4. Tests

### Validation with OpenGate Monte Carlo: Self-Dose (`test_madsen2018_self_dose`)

Located in [`tests/test_dose.py`](../tests/test_dose.py), this test compares the organ-level S-value calculation of `madsen2018` against a ground-truth Monte Carlo dose rate simulation performed with OpenGate (`opengate_dose_rate`) on the `p142` patient dataset:

1. **Simulation Input (`data/p142/`):**
   * Companion cropped CT: `ct_4mm_cropped.nii.gz`
   * Cropped liver activity: `liver_4mm_cropped.nii.gz` with uniform $^{177}\text{Lu}$ activity defined in `p142_dose_rate_test1.json` ($A_0 = 10^6\text{ Bq}$).
   * OpenGate simulated output: Instantaneous dose rate map `output_vrt/output_dose.mhd` (in $\text{Gy/s}$) and voxel uncertainty map `output_vrt/output_dose_uncertainty.mhd`.

2. **Integration & Comparison:**
   * The mean dose rate in the liver from OpenGate ($\dot{D}_{\text{Gate}}$) is time-integrated assuming mono-exponential clearance with effective half-life $T_{\text{eff}} = 45.0\text{ h}$:
     $$D_{\text{Gate}} = \dot{D}_{\text{Gate}} \cdot \frac{T_{\text{eff}} \cdot 3600}{\ln(2)}$$
   * This is compared against the absorbed dose computed by `compute_roi_dose(..., method="madsen2018")` using OpenDose S-values and CT-based patient organ mass scaling.
   * Self-dose is predominantly driven by $\beta^-$ particles ($\sim 95\%$), leading to very local deposition; the test asserts relative difference $< 5\%$ (actual agreement is $\approx 3.6\%$).
   * Also computes and prints the Monte Carlo ROI uncertainty estimate and bounds via `roi_dose_uncertainty_bounds` and `roi_dose_uncertainty_estimate`.
   * Falls back automatically to `tests/data/p142_opengate_benchmark.json` if raw large 3D images under `data/p142/` are absent.

### Validation with OpenGate Monte Carlo: Cross-Dose (`test_madsen2018_cross_dose`)

Located in [`tests/test_dose.py`](../tests/test_dose.py), this test validates the cross-irradiation from source organ `Liver` to target organ `Kidney_R` ($S(\text{Kidney\_R} \leftarrow \text{Liver})$):

1. **Physical Context:**
   * Cross-dose is $\sim 100\times$ lower than self-dose ($\sim 36\ \mu\text{Gy}$ vs $\sim 3.9\text{ mGy}$) and is $100\%$ deposited by penetrating photon emissions (gammas and bremsstrahlung).
   * Cross-dose S-values depend heavily on inter-organ distance and anatomy geometry ($1/r^2$). Due to morphological differences between the stylized reference phantom and the actual patient CT anatomy, differences around $\approx 40\%$ are expected.
   * The test verifies that the relative difference is $< 50\%$ (actual agreement $\approx 43\%$) and reports the kidney ROI uncertainty (lower bound, upper bound, and particle-range cluster estimate).

**Create / update the benchmark:**
```bash
# in folder data/p142 (not in git)
opengate_dose_rate p142_dose_rate_test1.json
./create_benchmark_json.py
```

**Run the validation tests:**
```bash
# Run both OpenGate validation tests
pytest -s tests/test_dose.py -k "madsen2018_self_dose or madsen2018_cross_dose"

# Run only the self-dose test
pytest -s tests/test_dose.py -k "madsen2018_self_dose"

# Run only the cross-dose test
pytest -s tests/test_dose.py -k "madsen2018_cross_dose"
```

