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

#### Decay Chain Cascade (`--cascade`)
For cascade alpha emitters (such as $^{225}\text{Ac}$, $^{223}\text{Ra}$, $^{212}\text{Pb}$), the parent radionuclide decays into radioactive progeny. When `--cascade` is active, the effective non-penetrating energy constant $\Delta_{\text{local, cascade}}$ is computed in secular equilibrium over all daughters in the decay chain:
$$
\Delta_{\text{local, cascade}} = \sum_{i \in \text{chain}} Y_i \cdot \Delta_{\text{local}}(i)
$$
where $Y_i$ is the cumulative branching yield for daughter $i$ sourced from the official **IAEA LiveChart of Nuclides REST API** (curated in `radot/data/decay_chains.json`, updated via maintainer script `scripts/update_decay_chains.py`).
For example, for $^{225}\text{Ac}$, $\Delta_{\text{local, cascade}} \approx 28.70\text{ MeV}$ (compared to $5.92\text{ MeV}$ for parent alone), matching full-chain Monte Carlo simulations (OpenGate) within $\approx 2.7\%$.

## 4. Tests

### Validation with OpenGate Monte Carlo: Self-Dose (`test_madsen2018_self_dose`, `test_hanscheid2017_self_dose`, `test_hanscheid2018_self_dose`, `test_ldm_self_dose`)

Located in [`tests/test_dose.py`](../tests/test_dose.py), these tests compare organ-level absorbed dose and dose rate calculations against a ground-truth Monte Carlo dose rate simulation performed with OpenGate (`opengate_dose_rate`) on the `p142` patient dataset ($A_0 = 10^6\text{ Bq}$, $T_{\text{eff}} = 45.0\text{ h}$):

1. **Simulation Input (`data/p142/`):**
   * Companion cropped CT: `ct_4mm_cropped.nii.gz`
   * Cropped liver activity: `liver_4mm_cropped.nii.gz` with uniform $^{177}\text{Lu}$ activity defined in `p142_dose_rate_test1.json` ($A_0 = 10^6\text{ Bq}$).
   * OpenGate simulated output: Instantaneous dose rate map `output_vrt/output_dose.mhd` (in $\text{Gy/s}$) and voxel uncertainty map `output_vrt/output_dose_uncertainty.mhd`.

2. **Integration & Comparison:**
   * The mean dose rate in the liver from OpenGate ($\dot{D}_{\text{Gate}}$) is time-integrated assuming mono-exponential clearance with effective half-life $T_{\text{eff}} = 45.0\text{ h}$:
     $$D_{\text{Gate}} = \dot{D}_{\text{Gate}} \cdot \frac{T_{\text{eff}} \cdot 3600}{\ln(2)}$$
   * **Madsen 2018 (`test_madsen2018_self_dose`)**: Uses OpenDose S-values and CT-based patient organ mass scaling. Self-dose is predominantly driven by $\beta^-$ particles ($\sim 95\%$), leading to very local deposition; relative difference is $< 5\%$ (actual agreement is $\approx 3.8\%$).
   * **Hänscheid 2017 (`test_hanscheid2017_self_dose`)**: Uses the empirical formula $D = 0.125 \cdot C_{\text{mean}} \cdot 2^{t/T_{\text{eff}}} \cdot T_{\text{eff}}$. Assumes unit tissue density ($\rho = 1.0\text{ g/mL}$) whereas the patient's liver in CT has $\bar{\rho} \approx 1.08\text{ g/mL}$, leading to an expected relative difference $< 15\%$ (actual agreement is $\approx 9.5\%$).
   * **Hänscheid 2018 (`test_hanscheid2018_self_dose`)**: Uses OpenDose S-values assuming $T_{\text{eff}} = 2 \cdot t$ ($t = 22.5\text{ h}$ for $T_{\text{eff}} = 45.0\text{ h}$). Yields identical results to Madsen 2018; relative difference is $< 5\%$ (actual agreement is $\approx 3.8\%$).
   * **Local Deposition Model (`test_ldm_self_dose`)**: Uses physical non-penetrating energy constant $\Delta_{\text{local}} = 0.14713\text{ MeV}$ with voxel-level CT density. Directly compares instantaneous dose rate $\dot{D}$ (and integrated dose); relative difference is $< 1\%$ (actual agreement is $\approx 0.29\%$).
   * All tests compute and print the Monte Carlo ROI uncertainty estimate and bounds via `roi_dose_uncertainty_bounds` and `roi_dose_uncertainty_estimate`.
   * Fall back automatically to `tests/data/p142_opengate_benchmark.json` if raw 3D images under `data/p142/` are absent.

### Validation with OpenGate Monte Carlo: Yttrium-90 Self-Dose (`test_y90_self_dose`)

Located in [`tests/test_dose.py`](../tests/test_dose.py), this test compares absorbed dose and dose rate calculations against an OpenGate Monte Carlo simulation for $^{90}\text{Y}$ ($A_0 = 10^6\text{ Bq}$, $T_{\text{eff}} = T_p = 64.1\text{ h}$ for radioembolization microspheres):

1. **Simulation Input (`data/p142/`):**
   * Configured in `p142_dose_rate_y90.json` with `"output_folder": "output_y90_vrt"`.
   * Mean $\beta^-$ energy of $^{90}\text{Y}$ is $\approx 0.93\text{ MeV}$ (compared to $\approx 0.15\text{ MeV}$ for $^{177}\text{Lu}$), with maximum beta range $\sim 11\text{ mm}$ in tissue.
2. **Comparison:**
   * **Local Deposition Model (`ldm`)**: Matches OpenGate within $1.8\%$ (slight boundary escape effect due to the $11\text{ mm}$ beta range extending across boundary voxels).
   * **Madsen 2018 (`madsen2018`)**: Matches OpenGate within $6.3\%$ using the OpenDose S-value for $^{90}\text{Y}$ and patient mass scaling.

### Validation with OpenGate Monte Carlo: Iodine-131 Self-Dose (`test_i131_self_dose`)

Located in [`tests/test_dose.py`](../tests/test_dose.py), this test compares absorbed dose calculations for $^{131}\text{I}$ ($A_0 = 5 \times 10^5\text{ Bq}$, $T_{\text{eff}} = T_p = 192.5\text{ h}$):

1. **Simulation Input (`data/p142/`):**
   * Configured in `p142_dose_rate_i131.json` with `"output_folder": "output_i131_vrt"`.
2. **Comparison:**
   * **Madsen 2018 (`madsen2018`)**: Matches OpenGate within **$0.43\%$** ($1.437 \times 10^{-2}\text{ Gy}$ vs $1.431 \times 10^{-2}\text{ Gy}$). OpenDose S-values account for both the $\beta^-$ spectrum and photon self-absorption.
   * **Local Deposition Model (`ldm`)**: Predicts $1.079 \times 10^{-2}\text{ Gy}$ ($\approx 25\%$ lower than OpenGate). This illustrates an important physical limitation of the LDM: $^{131}\text{I}$ emits energetic photons ($364\text{ keV}$, $81.7\%$ yield) which deposit $\approx 25\%$ of the total self-dose in a large organ like the liver ($1.4\text{ kg}$). Because LDM only models non-penetrating electrons ($\Delta_{\text{local}} = 0.192\text{ MeV}$), it excludes this photon self-dose component.

### Validation with OpenGate Monte Carlo: Actinium-225 Self-Dose (`test_ac225_self_dose`)

Located in [`tests/test_dose.py`](../tests/test_dose.py), this test compares alpha emitter $^{225}\text{Ac}$ calculations ($A_0 = 10^5\text{ Bq}$, $T_{\text{eff}} = T_p = 240.0\text{ h}$):

1. **Simulation Input (`data/p142/`):**
   * Configured in `p142_dose_rate_ac225.json` in analog mode (`"mode": "analog"`, `"output_folder": "output_ac225_vrt"`).
2. **Comparison:**
   * **Madsen 2018 vs LDM**: For the single $^{225}\text{Ac}$ parent decay ($\Delta_{\text{local}} = 5.92\text{ MeV}$), Madsen 2018 ($8.308 \times 10^{-2}\text{ Gy}$) and LDM ($8.301 \times 10^{-2}\text{ Gy}$) match each other within **$0.09\%$**.
   * **OpenGate Decay Cascade**: In analog mode, OpenGate tracks the complete decay chain cascade of $^{225}\text{Ac}$ down to stable $^{209}\text{Bi}$ ($4\alpha + 2\beta^-$: $^{225}\text{Ac} \to {}^{221}\text{Fr} \to {}^{217}\text{At} \to {}^{213}\text{Bi} \to {}^{213}\text{Po} \to {}^{209}\text{Pb} \to {}^{209}\text{Bi}$), depositing $\approx 28.0\text{ MeV}$ per $^{225}\text{Ac}$ decay in secular equilibrium. Consequently, the OpenGate integrated dose ($3.921 \times 10^{-1}\text{ Gy}$) is $\approx 4.72\times$ higher than the single-decay parent models ($28.0\text{ MeV} / 5.92\text{ MeV} \approx 4.73$), matching the physical cascade factor.

### Validation with OpenGate Monte Carlo: Cross-Dose (`test_madsen2018_cross_dose` & `test_hanscheid2018_cross_dose`)

Located in [`tests/test_dose.py`](../tests/test_dose.py), these tests validate cross-irradiation from source organ `Liver` to target organ `Kidney_R` ($S(\text{Kidney\_R} \leftarrow \text{Liver})$):

1. **Physical Context:**
   * Cross-dose is $\sim 100\times$ lower than self-dose ($\sim 36\ \mu\text{Gy}$ vs $\sim 3.9\text{ mGy}$) and is $100\%$ deposited by penetrating photon emissions (gammas and bremsstrahlung).
   * Cross-dose S-values depend heavily on inter-organ distance and anatomy geometry ($1/r^2$). Due to morphological differences between the stylized reference phantom and the actual patient CT anatomy, differences around $\approx 40\%$ are expected.
   * Both **Madsen 2018** and **Hänscheid 2018** cross-dose calculations verify that relative difference against OpenGate is $< 50\%$ (actual agreement $\approx 44.8\%$) and report the kidney ROI uncertainty.
   * *Note on Hänscheid 2017 & LDM:* Cross-dose is not defined for `hanscheid2017` or `ldm` because both methods assume non-penetrating local energy deposition within the activity-bearing volume ($C_{\text{mean}} = 0$ / zero activity in target organ yields $0\text{ Gy}$).

**Create / update the benchmarks:**
```bash
# In data/p142/ (not committed to git)
# Run simulations for desired radionuclides:
opengate_dose_rate p142_dose_rate_lu177.json
opengate_dose_rate p142_dose_rate_y90.json
opengate_dose_rate p142_dose_rate_i131.json
opengate_dose_rate p142_dose_rate_ac225.json

# Automatically discover all simulations and generate benchmark JSON files:
./create_benchmark_json.py
# Generates:
#   tests/data/p142_lu177_opengate_benchmark.json
#   tests/data/p142_y90_opengate_benchmark.json
#   tests/data/p142_i131_opengate_benchmark.json
#   tests/data/p142_ac225_opengate_benchmark.json
#   tests/data/p142_opengate_benchmark.json (Lu-177 default)
```

**Run the validation tests:**
```bash
# Run all OpenGate validation tests (9 tests)
pytest -s tests/test_dose.py -k "dose"

# Run individual tests
pytest -s tests/test_dose.py -k "test_madsen2018_self_dose"
pytest -s tests/test_dose.py -k "test_hanscheid2017_self_dose"
pytest -s tests/test_dose.py -k "test_hanscheid2018_self_dose"
pytest -s tests/test_dose.py -k "test_ldm_self_dose"
pytest -s tests/test_dose.py -k "test_y90_self_dose"
pytest -s tests/test_dose.py -k "test_i131_self_dose"
pytest -s tests/test_dose.py -k "test_ac225_self_dose"
pytest -s tests/test_dose.py -k "test_madsen2018_cross_dose"
pytest -s tests/test_dose.py -k "test_hanscheid2018_cross_dose"
```


