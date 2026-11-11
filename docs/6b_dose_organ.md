# Dosimetry: Organ Level

The organ-level methods give one mean absorbed dose per ROI (organ or lesion), from the total activity in the ROI. They need a label image (see the [overview](6a_dose.md) for the inputs, the outputs and the CLI).

---

## 1. Principles

* **Self-dose only**: the dose to an ROI comes from the activity in that ROI (the S-value of the organ to itself). The cross-doses from the other organs, and the photons that leave the organ, are not included; use a [Monte Carlo map](6c_dose_voxel.md#2-external-dose-rate-map-monte-carlo) for them. For Lu-177 the self-dose dominates in the kidneys, the liver and the spleen.
* **OpenDose S-values**: the self S-value $S$ of the organ (mGy/(MBq·s)) for the ICRP 110 adult male or female phantom (`--phantom`, default male), from the [OpenDose](5b_opendose.md) data (shipped for Lu-177, downloaded once from Zenodo for the other radionuclides).
* **Mass scaling**: the phantom S-value is scaled to the patient organ mass:
  $$
  S_{\text{patient}} = \frac{S}{\text{mass-scaling}}, \qquad \text{mass-scaling} = \frac{m_{\text{patient}}}{m_{\text{reference}}}
  $$
  where $m_{\text{patient}}$ is the ROI mass from the CT density ($\rho = 1 + \text{HU}/1000$; without CT, 1 g/mL) and $m_{\text{reference}}$ the OpenDose phantom organ mass.
* **ROI names**: each ROI name is mapped to an OpenDose organ (TotalSegmentator → TG-263 → OpenDose, see [ROI names](2b_roi_names.md)); the OpenDose name must match exactly, ignoring case and separators only: no partial matching, so a lesion named `Liver_met1` gets no liver S-value.
* **ROIs without an OpenDose organ** (lesions, custom structures): no S-value. When looping over all the labels, they are skipped with a warning. When requested with `--rois`, radot stops with an error. `hanscheid2017` needs no S-value and works for any ROI.

Then the dose is the time-integrated activity $\tilde{A}$ in the ROI (MBq·s) times the scaled S-value:

$$
D = \frac{\tilde{A} \cdot S}{\text{mass-scaling}} \cdot 10^{-3} \quad [\text{Gy}]
$$

The methods differ in how $\tilde{A}$ is obtained.

---

## 2. Methods

### Hänscheid 2017 (`hanscheid2017`)
A single-time-point method specific to $^{177}\text{Lu}$. It uses no S-value and no organ mass: the dose is the mean activity concentration times a constant (0.125 Gy·mL/(MBq·h), valid for Lu-177 only) and the effective half-life, so it applies to any ROI, lesions included.

$$
D = 0.125 \cdot C_{\text{mean}} \cdot 2^{\frac{t}{T_{\text{eff}}}} \cdot T_{\text{eff}}
$$

*   $D$: mean absorbed dose in the ROI (Gy).
*   $C_{\text{mean}}$: mean activity concentration in the ROI at the acquisition time (MBq/mL).
*   $t$: time from injection to acquisition (h).
*   $T_{\text{eff}}$: effective half-life in the ROI (h), `--teff`; default: the physical half-life.

---

### Hänscheid 2018 (`hanscheid2018`)
A single-time-point method with OpenDose S-values, which needs no effective half-life: it assumes $T_{\text{eff}} = 2t$ (Hänscheid et al. 2018), so that one SPECT scan gives the TIA.

$$
\tilde{A} = A(t) \cdot \frac{T_{\text{eff}} \cdot 3600}{\ln(2)}, \qquad T_{\text{eff}} = 2t \text{ by default}
$$

*   $A(t)$: total activity in the ROI at the acquisition time (MBq); $t$: time from injection (h).
*   $T_{\text{eff}}$ can be given per ROI (`--teff`).

---

### Madsen 2018 (`madsen2018`)
A single-time-point method with OpenDose S-values and a given effective half-life per ROI, assuming a mono-exponential clearance from the injection:

$$
\tilde{A} = A(t) \cdot \frac{e^{k \cdot t \cdot 3600}}{k}, \qquad k = \frac{\ln(2)}{T_{\text{eff}} \cdot 3600}
$$

*   $A(t)$: total activity in the ROI at the acquisition time (MBq); $t$: time from injection (h); $T_{\text{eff}}$: effective half-life in the ROI (h), `--teff`.

---

### S-value on a TIA image (`svalue`)
When the time-integrated activity is already known (for example a TIA map from another software, or the TIA maps of the SNMMI Dosimetry Challenge), $\tilde{A}$ is the sum of the TIA image in the ROI (MBq·s). No time point and no effective half-life are needed. The results also report the ROI TIA (`tia_mbq_s`, `tia_mbq_h`). See [TIA images](6a_dose.md#time-integrated-activity-tia-images).

---

### Jackson 2013 (`jackson2013`)
A **multiple-time-point** method (*Jackson et al., Med. Phys. 40(11), 112503 (2013)*). It needs **exactly 3 SPECT acquisitions** ($t_1 < t_2 < t_3$ after injection) and fits a tri-exponential biological clearance to the decay-corrected ROI activities:

$$
A_{\text{bio}}(t) = A_1 e^{k_1 t} + A_2 e^{k_2 t} + A_3 e^{k_3 t}
$$

with the constraints:
1. No activity at injection: $A_{\text{bio}}(0) = 0 \implies A_1 = -(A_2 + A_3)$.
2. Fixed uptake rate (typical small molecule / peptide kinetics): $k_1 = -1.3\text{ h}^{-1}$ ($T_{1/2} \approx 32\text{ min}$).
3. Log-linear stripping gives the intermediate $(A_2, k_2)$ and terminal $(A_3, k_3)$ components, with $k_3 \le -0.001\text{ h}^{-1}$ (no divergence).

The physical activity $A(t) = A_{\text{bio}}(t) e^{-\lambda_p t} = \sum_{i=1}^3 A_i e^{(k_i - \lambda_p) t}$ is integrated from $0$ to $\infty$:

$$
\tilde{A} = \sum_{i=1}^{3} \frac{A_i}{\lambda_p - k_i} \times 3600 \quad [\text{MBq}\cdot\text{s}]
$$

where $\lambda_p = \ln(2) / T_p$ is the physical decay constant (h⁻¹). The physical decay is separated from the biological kinetics, so the method applies to any radionuclide. The results report `tia_mbq_s`, `tia_mbq_h`, `roi_teff_h`, the ROI activities and times, and `fit_params`.

#### Configuration
The time points are given in a JSON file under `"timepoints"`.

**1. One segmentation per time point** (images not registered):

```json
{
    "method": "jackson2013",
    "isotope": "Lu177",
    "rois": ["Liver", "Kidney_L", "Kidney_R", "Spleen"],
    "timepoints": [
        {"image": "spect_tp1.nii.gz", "labels": "labels_tp1.nii.gz", "ct": "ct_tp1.nii.gz"},
        {"image": "spect_tp2.nii.gz", "labels": "labels_tp2.nii.gz", "ct": "ct_tp2.nii.gz"},
        {"image": "spect_tp3.nii.gz", "labels": "labels_tp3.nii.gz", "ct": "ct_tp3.nii.gz"}
    ],
    "output": "dose_jackson.json"
}
```

A per-time-point `labels` has priority over a global `labels`.

**2. One reference segmentation** (images registered to the first time point):

```json
{
    "method": "jackson2013",
    "isotope": "Lu177",
    "labels": "labels_tp1.nii.gz",
    "ct": "ct_tp1.nii.gz",
    "rois": ["Liver", "Kidney_L", "Kidney_R", "Spleen"],
    "timepoints": ["spect_tp1.nii.gz", "spect_tp2.nii.gz", "spect_tp3.nii.gz"],
    "output": "dose_jackson.json"
}
```

The times from injection are read from each SPECT sidecar (`injection_date`, `acquisition_date`), or given per time point (`{"image": "spect_tp1.nii.gz", "time_h": 4.99}`). radot does not register the time points: use one segmentation per time point, or register the images beforehand.

#### Python API

```python
from radot.dose import compute_roi_dose_mtp

results = compute_roi_dose_mtp(
    images=["spect_tp1.nii.gz", "spect_tp2.nii.gz", "spect_tp3.nii.gz"],
    labels=["labels_tp1.nii.gz", "labels_tp2.nii.gz", "labels_tp3.nii.gz"],
    ct_image=["ct_tp1.nii.gz", "ct_tp2.nii.gz", "ct_tp3.nii.gz"],
    method="jackson2013",
    isotope="Lu177",
    rois=["Kidney_L", "Kidney_R", "Liver", "Spleen"],
)
# per ROI: mean_dose_gy, tia_mbq_s, tia_mbq_h, roi_teff_h, fit_params
```
