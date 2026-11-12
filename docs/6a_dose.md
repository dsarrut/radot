# Dosimetry: Overview

radot computes absorbed doses (Gy), or dose rates (Gy/h), at two levels:

* **Organ level** ([details](6b_dose_organ.md)): one mean dose per ROI (organ or lesion), from the activity in the ROI. Methods `hanscheid2017`, `hanscheid2018`, `madsen2018`, `svalue`, `jackson2013`. Output: a table of ROI doses.
* **Voxel level** ([details](6c_dose_voxel.md)): a dose in every voxel, from the activity in that voxel. Method `ldm` (local deposition), or an external Monte Carlo dose-rate map (`--dose-rate-map`). Output: a 3D dose image, and the ROI mean doses when labels are given.

---

## 1. Two choices behind every method

A dose calculation combines two independent parts:

* **How the energy is deposited** (the dose-rate engine): organ S-values (OpenDose, self-dose), a constant (Hänscheid 2017), the local deposition of the non-penetrating emissions (LDM), or an external Monte Carlo map (photons and cross-doses included).
* **How the activity is integrated over time** (the time model): from one time point with an effective half-life (Hänscheid 2017/2018, Madsen 2018), physical decay only (LDM), a fit of several time points (Jackson 2013), or nothing when the image is already a time-integrated activity (TIA).

| Method | Level | Energy deposition | Time model | Input images | Output |
|---|---|---|---|---|---|
| `hanscheid2017` | organ | constant 0.125 (Lu-177 only) | one time point, $T_{\text{eff}}$ (default: physical half-life) | activity, labels | ROI table |
| `hanscheid2018` | organ | OpenDose self S-value, mass scaled | one time point, $T_{\text{eff}} = 2t$ (or given) | activity, labels (CT for the mass) | ROI table |
| `madsen2018` | organ | OpenDose self S-value, mass scaled | one time point, mono-exponential, given $T_{\text{eff}}$ | activity, labels (CT) | ROI table |
| `svalue` | organ | OpenDose self S-value, mass scaled | none (TIA image) | TIA, labels (CT) | ROI table |
| `jackson2013` | organ | OpenDose self S-value, mass scaled | 3 time points, tri-exponential fit | 3 activity images, labels (CT) | ROI table |
| `ldm` | voxel | local deposition (non-penetrating, optional decay chain) | physical decay, or none (TIA image) | activity or TIA (CT for the density) | dose image (+ ROI table with labels) |
| `--dose-rate-map` | voxel engine | external Monte Carlo dose rate | that of the method (`hanscheid2018`, `madsen2018`, `ldm`, `svalue`) | activity or TIA, the map (labels) | ROI table; dose image with `ldm` |

## 2. Which one?

* **One SPECT/PET time point, organ doses**: `hanscheid2018` (no effective half-life needed) or `madsen2018` (with a known $T_{\text{eff}}$ per ROI). `hanscheid2017` for Lu-177 when only the mean concentration is available.
* **Three time points**: `jackson2013`.
* **A TIA image from another software** (unit `bqs` or `bqsml`): `svalue` (organ) or `ldm` (voxel).
* **A dose map** (dose distribution, DVH, lesion heterogeneity): `ldm`, or a Monte Carlo map; view it and its DVH with [vvv](viewing.md).
* **Cross-doses, photons, patient anatomy**: the organ methods give the self-dose only (the S-value of the organ to itself), and LDM the local dose only. A Monte Carlo map (OpenGATE `opengate_dose_rate`) includes all of them.

The validation against Monte Carlo and the SNMMI challenge is in [OpenGATE validation](7_validation_opengate.md) and [SNMMI challenge](8_validation_snmmi.md).

## 3. Inputs and outputs

* **Activity image**: `bq` (per voxel), `bqml`, or `suv`/`sul` (converted to activity with the sidecar, see [units](3a_images_units.md)); time from injection from the sidecar dates. **TIA image**: unit `bqs` (Bq·s per voxel) or `bqsml` (Bq·s/mL); see [TIA images](#time-integrated-activity-tia-images).
* **Labels**: a label image and its sidecar (`{"1": "Liver", ...}`); names resolved to OpenDose organs (see [ROI names](2b_roi_names.md)). Required for the organ methods, optional for `ldm`. A label image without a readable sidecar is an error (the ROI names are unknown).
* **CT**: voxel density $\rho = 1 + \text{HU}/1000$, used for the ROI mass (S-value mass scaling, mass-weighted means) and the LDM voxel mass; with LDM, voxels below `--hu-threshold` (−800 HU, air) get no dose.
* **ROI table** (JSON): per ROI `mean_dose_gy`, `roi_mass_g`, `roi_volume_ml`, plus the method values (`roi_teff_h`, `tia_mbq_s`, `fit_params`, ...).
* **Dose image**: `unit` `gy` (or `gy/h` with `--dose-rate`) and `modality` `dose` in its sidecar.

`-o` sets the output: a `.json` file for the ROI table, an image file (`.nii.gz`, `.mhd`, ...) for the dose image. For `ldm` with labels, both are written: `-o dose.nii.gz` writes `dose.nii.gz` (+ its sidecar `dose.json`) and the table `dose_rois.json`; `-o results.json` writes `results.json` and the image `results_dose.nii.gz`. Without `-o`, the ROI table is printed; `ldm` without labels writes `dose.nii.gz`.

### Time-integrated activity (TIA) images

An image whose unit (sidecar `unit`) is `bqs` (Bq·s per voxel) or `bqsml` (Bq·s/mL) is a TIA image. It can be used with `svalue` (organ level) and `ldm` (voxel level). The single-time-point methods (`hanscheid2017`, `hanscheid2018`, `madsen2018`) reject a TIA image, `svalue` rejects an activity image, and `--dose-rate` and `--total-activity` cannot be used with a TIA image.

```bash
radot meta tia.nii.gz --modality spect --unit bqsml
radot dose tia.nii.gz -l labels.nii.gz --ct ct.nii.gz -m svalue --isotope lu177 --rois Liver --rois Kidney_L
radot dose tia.nii.gz -l labels.nii.gz --ct ct.nii.gz -m ldm --isotope lu177 -o dose.nii.gz
```

---

## 4. CLI

### Listing the methods: `radot dose-list`

Prints the equations of all the methods.

```bash
radot dose-list
```

### Computing doses: `radot dose`

```bash
radot dose <image_path_or_config_json> [options]
# or
radot dose -c <config_json> [options]
```

#### Options:
* -c, --config <config_path>: JSON configuration file (below).
* -l, --labels <labels_path>: Label image. Required for the organ methods, optional for `ldm`.
* --ct <ct_path>: CT image, for the density and the ROI masses.
* -m, --method <hanscheid2017|hanscheid2018|madsen2018|jackson2013|svalue|ldm>: Dosimetry method.
* --isotope <name>: Radionuclide (e.g. Lu177, Y90).
* --rois <name>: ROI to include. Can be repeated or comma-separated. Default: all the labels.
* --teff <ROI:value>: Effective half-life per ROI in hours (e.g. --teff Liver:45.0 --teff Kidney_L:32.0).
* --total-activity <float>: Total activity of the image in MBq (rescales the image; activity images only).
* --hu-threshold <float>: `ldm`: CT values below this (default −800 HU) get no dose (air).
* --dose-rate: `ldm` only: dose rate at the acquisition time (Gy/h) instead of the absorbed dose.
* --cascade / --no-cascade: `ldm`: include the decay chain in secular equilibrium (default: no).
* --dose-rate-map <path>: External dose-rate image (e.g. OpenGATE), see [voxel level](6c_dose_voxel.md#2-external-dose-rate-map-monte-carlo).
* --dose-rate-map-activity <float>: Activity (MBq) that produced the dose-rate map. Default: the map sidecar key `source_activity_mbq`, else the map is an absolute dose rate.
* --phantom <male|female>: OpenDose phantom for the S-value methods. Default: the patient sex in the image sidecar (`sex`); if neither is given, radot stops with an error (no default phantom). Not needed by `hanscheid2017`, `ldm` and a dose-rate map.
* --resample <like_activity|like_labels|like_ct>: Grid to resample the images to when they differ.
* -o, --output <output_path>: ROI table (`.json`) or dose image (see [inputs and outputs](#3-inputs-and-outputs)).

#### JSON configuration file

A dosimetry study can be described in a JSON file, for reproducibility:

```json
{
    "image": "spect.nii.gz",
    "labels": "labels.nii.gz",
    "ct": "ct.nii.gz",
    "method": "madsen2018",
    "radionuclide": "Lu177",
    "rois": {
        "Liver": 45.0,
        "Kidney_L": 32.0
    },
    "resample": "like_activity",
    "output": "dose_results.json"
}
```

*(The ROIs and their effective half-lives can also be given separately: `"rois": ["Liver", "Kidney_L"]` and `"teff": {"Liver": 45.0, "Kidney_L": 32.0}`.)* A Teff that is not a number of hours, or a `teff` item not in `ROI:hours` form, is an error. Relative paths are resolved from the folder of the JSON file. Multiple time points use `"timepoints"` (see [Jackson 2013](6b_dose_organ.md#jackson-2013-jackson2013)).

```bash
radot dose dose_config.json                         # positional
radot dose -c dose_config.json                      # explicit flag
radot dose dose_config.json -o custom_output.json   # any setting can be overridden
```

An option given on the command line (or an argument given to `compute_roi_dose_file` in Python) always wins over the config file, including `--no-cascade`, `--hu-threshold -800` and `--phantom` (also over a `phantom` in `method_params`).

#### Examples

```bash
# Organ level, one time point, effective half-lives per ROI
radot dose spect.nii.gz -l labels.nii.gz --ct ct.nii.gz -m madsen2018 --isotope Lu177 --teff Liver:45.0 --teff Kidney_L:32.0

# Voxel level: LDM dose image and ROI table
radot dose spect.nii.gz -l labels.nii.gz --ct ct.nii.gz -m ldm --isotope Lu177 -o dose.nii.gz

# Voxel level: LDM dose rate image (Gy/h)
radot dose spect.nii.gz --ct ct.nii.gz -m ldm --isotope Lu177 --dose-rate -o dose_rate.nii.gz
```

---

## 5. Python API

The public API is in `radot.dose`. The method-specific functions (`compute_roi_dose_mtp`, `dose_ldm`, the uncertainty estimates) are described in the [organ](6b_dose_organ.md) and [voxel](6c_dose_voxel.md) pages.

### `compute_roi_dose_file`
The full pipeline on files (reading, resampling, output), or from a JSON configuration file or dictionary.

```python
from pathlib import Path
from radot.dose import compute_roi_dose_file

results = compute_roi_dose_file(
    image_path=Path("spect.nii.gz"),
    labels_path=Path("labels.nii.gz"),
    method="madsen2018",
    isotope="Lu177",
    ct_path=Path("ct.nii.gz"),
    method_params={"Liver": {"teff": 45.0}},
    output_path=Path("dose.json"),
)

# Or from a JSON configuration file:
results = compute_roi_dose_file("dose_config.json")
```

### `compute_roi_dose`
The same on SimpleITK images in memory.

```python
import SimpleITK as sitk
from radot.dose import compute_roi_dose

spect = sitk.ReadImage("spect.nii.gz")
ct = sitk.ReadImage("ct.nii.gz")
labels = sitk.ReadImage("labels.nii.gz")

# Organ method: ROI table
results = compute_roi_dose(image=spect, method="hanscheid2018", isotope="Lu177", labels=labels, ct_image=ct,
                           labels_map={1: "Liver", 2: "Kidney_L"})

# Voxel method with labels: ROI table and dose image
results, dose_image = compute_roi_dose(image=spect, method="ldm", isotope="Lu177", labels=labels, ct_image=ct,
                                       labels_map={1: "Liver", 2: "Kidney_L"})
```

---

## 6. Tests and validation

Unit tests of the dosimetry models (no patient data):
```bash
pytest -q tests/test_dose.py tests/test_dose_jackson.py tests/test_dose_config.py \
    tests/test_dose_tia.py tests/test_dose_rate_map.py tests/test_ldm.py
```

The methods (LDM, Madsen 2018, Hänscheid 2017/2018, dose-rate map) are validated against OpenGATE Monte Carlo simulations (`opengate_dose_rate`) on a patient CT (p142), for $^{177}\text{Lu}$, $^{90}\text{Y}$, $^{131}\text{I}$ and $^{225}\text{Ac}$: self-dose, cross-doses and decay chain. Reference values: `validation/opengate/p142/p142_*_opengate_benchmark.json`; results: [OpenGATE validation](7_validation_opengate.md); how the references were produced: `validation/opengate/README.md`. The S-value, LDM and OpenGATE engines are also compared with the participants of the SNMMI 177Lu Dosimetry Challenge: [SNMMI challenge](8_validation_snmmi.md).

```bash
pytest -q validation                       # stored reference values, no data needed (CI)
pytest -q validation --opengate=images     # radot on the real images (data from Zenodo)
pytest -q validation --opengate=simulate   # reruns the OpenGATE simulations (hours of CPU)
```
