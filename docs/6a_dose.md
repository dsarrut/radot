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

Options can be supplied directly via command-line arguments or through a JSON configuration file:

```bash
radot dose <image_path_or_config_json> [options]
# or
radot dose -c <config_json> [options]
```

#### Options:
* -c, --config <config_path>: Path to a JSON configuration file containing dosimetry settings.
* -l, --labels <labels_path>: Optional anatomical segmentation label map. Required for S-value methods. Optional for LDM.
* --ct <ct_path>: Companion CT image used for physical density estimation and mass calibration.
* -m, --method <hanscheid2017|hanscheid2018|madsen2018|jackson2013|ldm>: Dosimetry model to apply.
* --isotope <name>: Radionuclide identifier (e.g. Lu177, Y90) to retrieve physical yield constraints.
* --rois <name>: ROI name to include. Can be specified multiple times or comma-separated.
* --teff <ROI:value>: Effective half-life override per ROI in hours (e.g. --teff Liver:45.0 --teff Kidney_L:32.0).
* --total-activity <float>: Global calibration activity override in MBq.
* --hu-threshold <float>: Hounsfield Units cutoff (defaults to -800.0). Voxel doses corresponding to CT values below this limit are zeroed out (air masking).
* --dose-rate: If set, computes instantaneous dose rate (in Gy/h) at the time of acquisition rather than integrated absorbed dose (in Gy). Only available for LDM.
* -o, --output <output_path>: Path to output NIfTI dose map or JSON summary results. Required for LDM if running without labels.
* --resample <like_activity|like_labels|like_ct>: Specifies the grid alignment strategy if input image spaces differ.
* --phantom <male|female>: Reference phantom model used for S-value calculations. Defaults to male.
* --cascade / --no-cascade: Account for radioactive decay cascade in secular equilibrium (default: false).

#### Using a JSON Configuration File:

Dosimetry studies can be defined in a reproducible JSON file:
```json
{
    "image": "spect.nii.gz",
    "labels": "labels.nii.gz",
    "ct": "ct.nii.gz",
    "method": "hanscheid2017",
    "radionuclide": "Lu177",
    "total_activity": 250.0,
    "rois": {
        "Liver": 45.0,
        "Kidney_L": 32.0
    },
    "resample": "like_activity",
    "output": "dose_results.json"
}
```
*(Alternatively, you can separate `"rois": ["Liver", "Kidney_L"]` and `"teff": {"Liver": 45.0, "Kidney_L": 32.0}`)*
Relative paths inside the JSON file are automatically resolved relative to the directory containing the JSON file.

Run with:
```bash
# Direct positional syntax
radot dose dose_config.json

# Explicit flag
radot dose -c dose_config.json

# Override any setting from CLI
radot dose dose_config.json -o custom_output.json -v
```

#### CLI Examples:

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
Executes the full dosimetry pipeline on file paths, or directly from a configuration file/dictionary.

```python
from pathlib import Path
from radot.dose import compute_roi_dose_file

# From arguments:
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

# Or directly from a JSON configuration file:
results = compute_roi_dose_file("dose_config.json")
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

### `compute_roi_dose_mtp`
Executes multiple timepoint (MTP) dosimetry (e.g. Jackson 2013 tri-exponential model) using a list of SimpleITK images or file paths:

```python
from radot.dose import compute_roi_dose_mtp

# Supports a list of 3 SPECT images, with per-timepoint labels or a registered reference label
results = compute_roi_dose_mtp(
    images=["spect_tp1.nii.gz", "spect_tp2.nii.gz", "spect_tp3.nii.gz"],
    labels=["labels_tp1.nii.gz", "labels_tp2.nii.gz", "labels_tp3.nii.gz"],
    ct_image=["ct_tp1.nii.gz", "ct_tp2.nii.gz", "ct_tp3.nii.gz"],
    method="jackson2013",
    isotope="Lu177",
    rois=["Kidney_L", "Kidney_R", "Liver", "Spleen"]
)

# Returns dictionary with mean_dose_gy, tia_mbq_s, tia_mbq_h, roi_teff_h, and fit_params per ROI
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
    \tilde{A} = A(t) \cdot \frac{T_{\text{eff}} \cdot 3600}{\ln(2)}
    $$
    where $T_{\text{eff}} = 2 \cdot t$ by default (as established in Hänscheid et al. 2018), or an optional user-provided effective half-life per ROI. $A(t)$ is the total activity in the ROI (in MBq) at acquisition time, and $t$ is the time from injection (in hours).
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

### Jackson 2013 (`jackson2013`)
The Jackson 2013 method (*Jackson et al., Med. Phys. 40(11), 112503 (2013)*) is a **multiple-timepoint (MTP)** organ-level dosimetry method. It requires **exactly 3 SPECT acquisitions** ($t_1 < t_2 < t_3$ post-injection) and fits an analytical tri-exponential biological clearance curve to the decay-corrected time-activity points:

$$
A_{\text{bio}}(t) = A_1 e^{k_1 t} + A_2 e^{k_2 t} + A_3 e^{k_3 t}
$$

with physiological constraints:
1. Physical uptake at injection is zero: $A_{\text{bio}}(0) = 0 \implies A_1 = -(A_2 + A_3)$.
2. The initial uptake rate is fixed based on typical small molecule/peptide kinetics: $k_1 = -1.3\text{ h}^{-1}$ ($T_{1/2} \approx 32\text{ min}$).
3. Stepwise log-linear stripping derives the intermediate $(A_2, k_2)$ and terminal $(A_3, k_3)$ clearance components, with non-divergence boundary enforcement ($k_3 \le -0.001\text{ h}^{-1}$).

#### Time-Integrated Activity (TIA)
The physical activity curve $A(t) = A_{\text{bio}}(t) e^{-\lambda_p t} = \sum_{i=1}^3 A_i e^{(k_i - \lambda_p) t}$ is analytically integrated from $t = 0 \to \infty$:

$$
\tilde{A} = \sum_{i=1}^{3} \frac{A_i}{\lambda_p - k_i} \times 3600 \quad [\text{MBq}\cdot\text{s}]
$$

where $\lambda_p = \ln(2) / T_p$ is the physical decay constant (in $\text{h}^{-1}$). This formulation decouples physical decay from biological kinetics, making the implementation applicable to any radionuclide.

#### Absorbed Dose
Organ absorbed dose is calculated using OpenDose reference $S$-values and patient-specific organ mass scaling:

$$
D = \frac{\tilde{A} \cdot S}{\text{mass-scaling}} \cdot 10^{-3} \quad [\text{Gy}]
$$

where $\text{mass-scaling} = m_{\text{patient}} / m_{\text{reference}}$, and $S$ is in $\text{mGy}/(\text{MBq}\cdot\text{s})$.

#### Multi-Timepoint Configuration Example
MTP studies are configured in a JSON file specifying the sequence of SPECT acquisitions under `"timepoints"`.

**1. Independent segmentations per timepoint (recommended when images are not co-registered):**
Each timepoint can define its own companion `labels` and `ct`. If both a global `labels` and per-timepoint `labels` are provided, the per-timepoint segmentation takes precedence:

```json
{
    "method": "jackson2013",
    "isotope": "Lu177",
    "rois": ["Liver", "Kidney_L", "Kidney_R", "Spleen"],
    "timepoints": [
        {
            "image": "spect_tp1.nii.gz",
            "labels": "labels_tp1.nii.gz",
            "ct": "ct_tp1.nii.gz"
        },
        {
            "image": "spect_tp2.nii.gz",
            "labels": "labels_tp2.nii.gz",
            "ct": "ct_tp2.nii.gz"
        },
        {
            "image": "spect_tp3.nii.gz",
            "labels": "labels_tp3.nii.gz",
            "ct": "ct_tp3.nii.gz"
        }
    ],
    "output": "dose_jackson.json"
}
```

**2. Single reference segmentation (when images are pre-registered to TP1):**
```json
{
    "method": "jackson2013",
    "isotope": "Lu177",
    "labels": "labels_tp1.nii.gz",
    "ct": "ct_tp1.nii.gz",
    "rois": ["Liver", "Kidney_L", "Kidney_R", "Spleen"],
    "timepoints": [
        "spect_tp1.nii.gz",
        "spect_tp2.nii.gz",
        "spect_tp3.nii.gz"
    ],
    "output": "dose_jackson.json"
}
```
*Note: Acquisition delays can be inferred automatically from each SPECT image sidecar (`injection_date` and `acquisition_date`), or explicitly specified in each timepoint dictionary (e.g. `{"image": "spect_tp1.nii.gz", "time_h": 4.99}`).*

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
    $$
    1\text{ MeV} = 1.60218 \times 10^{-13}\text{ J} \implies 1\text{ MeV/g} = 1.60218 \times 10^{-10}\text{ J/kg}
    $$

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

#### Decay Chain Cascade (`--cascade`)
For cascade alpha emitters (such as $^{225}\text{Ac}$, $^{223}\text{Ra}$, $^{212}\text{Pb}$), the parent radionuclide decays into radioactive progeny. When `--cascade` is active, the effective non-penetrating energy constant $\Delta_{\text{local, cascade}}$ is computed in secular equilibrium over all daughters in the decay chain:
$$
\Delta_{\text{local, cascade}} = \sum_{i \in \text{chain}} Y_i \cdot \Delta_{\text{local}}(i)
$$
where $Y_i$ is the cumulative branching yield for daughter $i$ sourced from the official **IAEA LiveChart of Nuclides REST API** (curated in `radot/data/decay_chains.json`, updated via maintainer script `scripts/update_decay_chains.py`).
For example, for $^{225}\text{Ac}$, $\Delta_{\text{local, cascade}} \approx 28.70\text{ MeV}$ (compared to $5.92\text{ MeV}$ for parent alone), matching full-chain Monte Carlo simulations (OpenGate) within $\approx 2.7\%$.

## 4. Tests

Dosimetry models (LDM, Madsen 2018, Hänscheid 2017/2018) are validated against reference OpenGate Monte Carlo simulations (`opengate_dose_rate`) on patient CT data (`p142`) across multiple radionuclides ($^{177}\text{Lu}$, $^{90}\text{Y}$, $^{131}\text{I}$, $^{225}\text{Ac}$), covering self-dose, cross-dose, and decay chain cascades. Tests compare computed doses against precomputed benchmark JSON files (`tests/data/p142_*_opengate_benchmark.json`).

### Create / Update Benchmarks

Run OpenGate simulations and generate the benchmark JSON files (in `data/p142/`):
```bash
opengate_dose_rate p142_dose_rate_lu177.json
opengate_dose_rate p142_dose_rate_y90.json
opengate_dose_rate p142_dose_rate_i131.json
opengate_dose_rate p142_dose_rate_ac225.json

./create_benchmark_json.py
```

### Run Validation Tests

```bash
# Run all dose validation tests
pytest -s tests/test_dose.py -k "dose"

# Run Jackson 2013 MTP tests
pytest -v tests/test_dose_jackson.py

# Run individual tests
pytest -s tests/test_dose.py -k "test_ldm_self_dose"
pytest -s tests/test_dose.py -k "test_madsen2018_self_dose"
pytest -s tests/test_dose.py -k "test_hanscheid2017_self_dose"
pytest -s tests/test_dose.py -k "test_hanscheid2018_self_dose"
pytest -s tests/test_dose.py -k "test_hanscheid2018_custom_teff"
pytest -s tests/test_dose.py -k "test_y90_self_dose"
pytest -s tests/test_dose.py -k "test_i131_self_dose"
pytest -s tests/test_dose.py -k "test_ac225_self_dose"
pytest -s tests/test_dose.py -k "test_madsen2018_cross_dose"
pytest -s tests/test_dose.py -k "test_hanscheid2018_cross_dose"
```
