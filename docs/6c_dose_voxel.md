# Dosimetry: Voxel Level

The voxel-level methods compute a dose in every voxel, from the activity in that voxel: a 3D dose image (Gy, or Gy/h for a dose rate), for dose distributions, dose-volume histograms and heterogeneous uptake (lesions). With a label image, radot also gives the mean dose of each ROI. See the [overview](6a_dose.md) for the inputs, the outputs and the CLI.

Two engines:

* **Local Deposition Model** (`ldm`): the energy of the non-penetrating emissions (electrons, α) is deposited in the voxel where the decay happens. No photons, no cross-dose: a good approximation for β and α emitters in large uptake regions.
* **External dose-rate map** (`--dose-rate-map`): a dose-rate image from a Monte Carlo simulation (e.g. OpenGATE), with the photons, the cross-doses and the patient anatomy.

**Output**: the dose image has `unit` `gy` (or `gy/h`) and `modality` `dose` in its sidecar. With labels, both the image and the ROI table are written (`-o dose.nii.gz` → `dose.nii.gz` and `dose_rois.json`). **ROI mean**: the mass-weighted mean of the voxel doses in the ROI, $\bar{D} = \sum_i D_i m_i / \sum_i m_i$ (voxel masses from the CT density).

**Viewing**: `vvv ct.nii.gz,dose.nii.gz + labels.nii.gz` shows the dose map with isodoses and the DVH of each ROI (see [Viewing with vvv](viewing.md)).

---

## 1. Local Deposition Model (`ldm`)

### Absorbed dose (Gy)

$$
D(x,y,z) = \frac{\tilde{A}(x,y,z) \cdot \Delta_{\text{local}} \cdot 1.60218 \times 10^{-10}}{m(x,y,z)}
$$

*   $D(x,y,z)$: absorbed dose in the voxel (Gy).
*   $\tilde{A}(x,y,z)$: time-integrated activity in the voxel (Bq·s). From an activity image, physical decay only after the acquisition (no biological clearance):
    $$
    \tilde{A}(x,y,z) = A(t, x,y,z) \cdot \frac{e^{\lambda t}}{\lambda_{\text{s}}}
    $$
    where $A(t, x,y,z)$ is the voxel activity at the acquisition time (Bq), $t$ the time from injection (h), $\lambda = \ln(2)/T_p$ (h⁻¹) and $\lambda_{\text{s}} = \lambda / 3600$ (s⁻¹). From a TIA image (`bqs`, `bqsml`), $\tilde{A}$ is the image itself.
*   $\Delta_{\text{local}}$: non-penetrating energy per decay of the radionuclide (MeV), from ICRP 107 (see [radionuclides](5a_radionuclide.md)).
*   $m(x,y,z) = V_{\text{voxel}} \cdot \rho(x,y,z)$: voxel mass (g), with the density from the CT, $\rho = 1 + \text{HU}/1000$ (at least 0.001 g/mL); without CT, 1 g/mL. Voxels below `--hu-threshold` (−800 HU) get no dose.
*   $1.60218 \times 10^{-10}$: MeV/g to Gy ($1\text{ MeV} = 1.60218 \times 10^{-13}\text{ J}$).

Since only the physical decay is used, LDM on an activity image overestimates the dose where the activity clears biologically; on a TIA image it uses the measured kinetics.

### Dose rate (Gy/h)

With `--dose-rate`, the dose rate at the acquisition time:

$$
\dot{D}(x,y,z) = \frac{A(t, x,y,z) \cdot \Delta_{\text{local}} \cdot 1.60218 \times 10^{-10} \cdot 3600}{m(x,y,z)}
$$

### Decay chain (`--cascade`)
For α emitters with radioactive progeny ($^{225}\text{Ac}$, $^{223}\text{Ra}$, $^{212}\text{Pb}$), `--cascade` sums the non-penetrating energy of all the daughters in secular equilibrium:
$$
\Delta_{\text{local, cascade}} = \sum_{i \in \text{chain}} Y_i \cdot \Delta_{\text{local}}(i)
$$
where $Y_i$ is the cumulative branching yield of daughter $i$ from the IAEA LiveChart of Nuclides (curated in `radot/data/decay_chains.json`, updated with `scripts/update_decay_chains.py`). For $^{225}\text{Ac}$, $\Delta_{\text{local, cascade}} \approx 28.70\text{ MeV}$ (5.92 MeV for the parent alone), within about 2.7 % of a full-chain OpenGATE simulation. Every unstable daughter must have ICRP 107 data: if one has none, `--cascade` stops with an error naming it (it is never dropped from the sum).

Only the daughters in **secular equilibrium** with the parent are counted: half-life below one tenth of the parent's, and every nuclide between them counted too. A longer-lived daughter and its own daughters are excluded, with a WARNING: e.g. for $^{149}$Tb (4.1 h) the daughters $^{149}$Gd (9.3 d) and $^{145}$Eu (5.9 d) are excluded, the cascade is the parent alone; for $^{211}$At the $^{211}$Po alphas are counted, $^{207}$Bi (31.6 y) is not. See [secular equilibrium](5a_radionuclide.md#which-daughters-are-counted-secular-equilibrium).

### Examples

```bash
radot dose spect.nii.gz --ct ct.nii.gz -m ldm --isotope Lu177 -o dose.nii.gz                       # dose image
radot dose spect.nii.gz -l labels.nii.gz --ct ct.nii.gz -m ldm --isotope Lu177 -o dose.nii.gz      # + ROI table
radot dose spect.nii.gz --ct ct.nii.gz -m ldm --isotope Lu177 --dose-rate -o dose_rate.nii.gz      # Gy/h
radot dose spect.nii.gz --ct ct.nii.gz -m ldm --isotope Ac225 --cascade -o dose.nii.gz             # decay chain
```

### Python API: `dose_ldm`
The core calculation on numpy arrays:

```python
import numpy as np
from radot.dose import dose_ldm

dose_array = dose_ldm(
    image_arr=np.ones((64, 64, 64)),      # activity per voxel (Bq)
    voxel_volume_ml=0.008,                # voxel volume (mL)
    time_from_injection_h=24.0,
    isotope="Lu177",
    ct_arr=np.zeros((64, 64, 64)),        # CT (HU)
    hu_threshold=-800.0,
    dose_rate=False,
)
```

---

## 2. External dose-rate map (Monte Carlo)

Every single-time-point method is a dose rate at the acquisition time multiplied by a time-integration factor. The dose rate comes from the S-value (`hanscheid2018`, `madsen2018`) or the local deposition (`ldm`). With `--dose-rate-map`, it comes from an external dose-rate image instead, for example a Monte Carlo simulation with [OpenGATE](https://github.com/OpenGATE/opengate) (`opengate_dose_rate`), which includes the photons, the cross-doses and the patient anatomy.

*   The map is the dose rate produced by an activity $A_{\text{map}}$ (`--dose-rate-map-activity`, MBq) distributed like the radot image. It is scaled to the image total: $\dot{D}(t) = \text{map} \cdot A_{\text{image}} / A_{\text{map}}$. One simulation (e.g. 1 MBq) is reused for any activity.
*   Unit: sidecar `unit` `gy/s` or `gy/h`; without a sidecar, `gy/s` (OpenGATE output).
*   The map must be on the grid of the activity image. ROI value: mass-weighted mean of the map in the ROI (CT density). No HU threshold.

| Method | Image | Dose | Output |
|---|---|---|---|
| `hanscheid2018` | activity | $\dot{D}(t) \cdot T_{\text{eff}} / \ln 2$, $T_{\text{eff}} = 2t$ by default | ROI table |
| `madsen2018` | activity | $\dot{D}(t) \cdot e^{kt} / k$, $k = \ln 2 / T_{\text{eff}}$ | ROI table |
| `ldm` | activity | $\dot{D}(t) \cdot e^{\lambda t} / \lambda$ (physical decay); with `--dose-rate`, $\dot{D}(t)$ in Gy/h | dose image (+ ROI table) |
| `svalue`, `ldm` | TIA | $\text{map} \cdot \tilde{A}_{\text{image}} / A_{\text{map}}$ (no time model) | ROI table; dose image with `ldm` |

With `hanscheid2018` and `madsen2018`, the time integration uses the effective half-life of each ROI, so the result is per ROI only. `hanscheid2017`, `jackson2013` and multiple time points reject a dose-rate map. The results add `mean_dose_rate_gys` (activity image) and `dose_engine: "map"`.

```bash
opengate_dose_rate config.json      # source: spect.nii.gz, CT: ct.nii.gz, activity_bq: 1e6
radot dose spect.nii.gz -l labels.nii.gz --ct ct.nii.gz -m hanscheid2018 --isotope lu177 \
    --dose-rate-map output/output_dose.mhd --dose-rate-map-activity 1
```

For a pure β emitter such as $^{90}\text{Y}$, the dose outside the source comes only from the bremsstrahlung: the electron run of `opengate_dose_rate` (VRT mode) must create it (option `bremsstrahlung`, on by default since OpenGATE PR #1223; older versions give cross-doses 10 to 100 times too low). For $^{177}\text{Lu}$ and $^{131}\text{I}$ it is less than 1 % of the cross-doses (γ lines dominate) and not required.

---

## 3. Statistical uncertainty of a Monte Carlo ROI dose

`roi_dose_uncertainty_bounds` and `roi_dose_uncertainty_estimate` estimate the statistical uncertainty of the mean dose of an ROI from a Monte Carlo dose image and its voxel relative uncertainty image.

```python
import SimpleITK as sitk
from radot.dose import roi_dose_uncertainty_bounds, roi_dose_uncertainty_estimate

dose_img = sitk.ReadImage("output_dose.mhd")
unc_img = sitk.ReadImage("output_dose_uncertainty.mhd")
mask_img = sitk.ReadImage("liver_mask.nii.gz")

# Theoretical lower bound (independent voxels) and upper bound (100% correlation)
unc_min, unc_max = roi_dose_uncertainty_bounds(dose_img, unc_img, mask_img)

# Realistic estimate: one electron history correlates the voxels within the mean beta range
unc_est = roi_dose_uncertainty_estimate(dose_img, unc_img, mask_img, isotope="Lu177")

print(f"ROI Uncertainty: {unc_est * 100:.2f}% (bounds: [{unc_min * 100:.3f}%, {unc_max * 100:.2f}%])")
```

One history deposits in several voxels, so the voxel uncertainties are correlated. The estimate is `unc_min * sqrt(nu)`, where `nu` is the number of ROI voxels reached by one history:

- electrons (default): a cluster of the mean β range of the isotope, `nu = prod(1 + range / spacing)`;
- photons (`photons=True`, e.g. a track-length estimator dose map): a photon deposits along its path across the ROI, `nu = 1 + chord / spacing`, with the ROI mean chord `4 V / S` (`roi_mean_chord_mm`).

For a dose map with both, estimate the electron and photon parts separately when they are available (e.g. the `e/` and `gamma_tle/` outputs of an OpenGATE VRT simulation) and add the variances. On the merged map, the electron estimate underestimates the uncertainty of the photon dose (for a cross-dose, by a factor of about 2). Checked against repeated simulations in `validation/opengate/` (see the [validation page](7_validation_opengate.md)).

```python
from radot.dose import roi_dose_uncertainty_estimate

parts = [("e", False), ("gamma_tle", True)]
variance, total = 0.0, 0.0
for folder, photons in parts:
    dose = sitk.ReadImage(f"output/{folder}/output_dose.mhd")
    unc = sitk.ReadImage(f"output/{folder}/output_dose_uncertainty.mhd")
    mean = sitk.GetArrayViewFromImage(dose)[sitk.GetArrayViewFromImage(mask_img) > 0].mean()
    u = roi_dose_uncertainty_estimate(dose, unc, mask_img, isotope="Lu177", photons=photons)
    variance += (u * mean) ** 2
    total += mean
unc_organ = variance ** 0.5 / total
```
