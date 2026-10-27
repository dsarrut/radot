# OpenGATE validation

Comparison of radot organ doses with OpenGATE Monte Carlo simulations on a patient CT (p142).
The results are published in [docs/7_validation.md](../../docs/7_validation.md).

| File | Role |
|---|---|
| `benchmarks.py` | Cases, acceptance ranges and computations (single source of truth) |
| `test_opengate.py` | The same comparisons as tests (collected by `pytest`) |
| `report.py` | Writes `docs/7_validation.md` |
| `p142/prepare/` | Preparation of the simulation inputs |
| `p142/p142_dose_rate_*.json` | OpenGATE configuration, one per radionuclide |
| `p142/create_benchmark_json.py` | OpenGATE outputs → `p142/p142_*_opengate_benchmark.json` |
| `p142/p142_*_opengate_benchmark.json` | Reference values used by the tests and the report |
| `batch_uncertainty.py` | Measures the uncertainty of the organ means with repeated simulations (step 6 below) |
| `p142/p142_lu177_batch_uncertainty.json` | Result of `batch_uncertainty.py` for Lu-177 |
| `download_data.py` | Downloads the images and OpenGATE outputs from Zenodo into `p142/data/` |
| `p142/data/` | Images and OpenGATE outputs (not in git, downloaded on demand) |

## Run

```bash
pytest validation                            # radot vs the stored OpenGATE values (no download)
pytest validation --opengate=images          # radot on the real images vs the OpenGATE images
pytest validation --opengate=simulate        # rerun the 4 OpenGATE simulations, then compare
python validation/opengate/report.py         # regenerate docs/7_validation.md
```

- Default level: the organ is modelled as a single voxel with the volume, mass and activity of the
  benchmark JSON. It needs no image and runs in CI.
- `images` level: radot runs on the real CT, liver and labels, and the OpenGATE dose rate is averaged
  over the same masks. The data (34 MB) are downloaded on first use (`python validation/opengate/download_data.py`
  does the same) and checked against their SHA-256. It also checks that the stored JSON values are those
  of the OpenGATE images. On p142, both levels agree within 0.02 percentage points.
- `simulate` level: needs opengate (`opengate_dose_rate`). Runs the 4 configurations on the downloaded
  inputs, with outputs in `p142/data/simulated/` (the published outputs are not modified), then checks
  that the new OpenGATE values agree with the stored ones within 3 combined standard uncertainties, and
  runs the radot comparisons against the new outputs. Each output folder also contains the OpenGATE log
  and `opengate_version.json` (version and git commit), which `create_benchmark_json.py` copies into the
  benchmark JSON.

`pytest` fails if `docs/7_validation.md` is not up to date.

## How the reference was produced

All commands run in `p142/data/`.

1. Inputs, resampled to 4 mm (`labels_all.nii.gz` is the organ segmentation of the CT):
   ```bash
   radot resample -t ct ct.nii.gz -o ct_4mm.nii.gz -s 4 4 4
   radot resample -t labels labels_all.nii.gz -o labels_all_4mm.nii.gz -s 4 4 4
   radot meta ct_4mm.nii.gz --modality ct
   ```
2. Activity source: uniform activity in the liver (`liver_4mm.nii.gz`):
   ```bash
   ../prepare/extract_liver.py
   ```
3. Crop along Z, centred on the liver, to twice the liver height
   (`ct_4mm_cropped.nii.gz`, `liver_4mm_cropped.nii.gz`):
   ```bash
   ../prepare/crop_ct_liver.py
   ```
4. One OpenGATE simulation per radionuclide (1 MBq in the liver, dose rate at t = 0):
   ```bash
   opengate_dose_rate ../p142_dose_rate_lu177.json
   opengate_dose_rate ../p142_dose_rate_y90.json
   opengate_dose_rate ../p142_dose_rate_i131.json
   opengate_dose_rate ../p142_dose_rate_ac225.json
   ```
   Lu-177, Y-90 and I-131 use variance reduction (`"mode": "vrt"`); Ac-225 is simulated in analog mode.
5. Mean dose rate and uncertainty in the liver (self-dose) and the right kidney (cross-dose, label 6):
   ```bash
   python ../create_benchmark_json.py
   ```
   It also crops `labels_all_4mm.nii.gz` to `labels_all_4mm_cropped.nii.gz` (slices 282–366) if needed.
6. Uncertainty of the organ means. OpenGATE gives a history-by-history uncertainty per voxel, not for
   an organ, and one history deposits in several voxels, which are therefore correlated.
   `benchmarks.gate_roi_uncertainty` estimates the organ uncertainty from the voxel maps, separately for
   the two VRT simulations: electrons (`e/`, a cluster of the mean β range) and photons (`gamma_tle/`,
   the mean chord 4V/S of the organ, since a photon deposits along its path). This estimate was checked
   with the batch method: the Lu-177 simulation is repeated 20 times with different random seeds, each
   identical to the reference (1 MBq), and the relative standard deviation of the 20 organ means is the
   uncertainty of one simulation (about 1 min per batch on 8 threads):
   ```bash
   python ../../batch_uncertainty.py --run   # batches in batches/run_NN/, result in p142_lu177_batch_uncertainty.json
   ```
   | Lu-177 (20 batches) | Measured        | Estimated | Former estimate (β range only, merged map) |
   |---------------------|-----------------|-----------|---------------------------------------------|
   | Liver (self-dose)   | 0.127 % ± 0.021 | 0.126 %   | 0.20 % |
   | Right kidney (cross-dose, photons) | 1.06 % ± 0.17 | 0.94 % | 0.58 % |

   The mean of the 20 batches agrees with the reference simulation: +0.23 % ± 0.13 (liver) and
   −0.52 % ± 1.09 (kidney), including the uncertainty of the reference. `pytest` checks that the
   estimates in the benchmark JSONs agree with the measured values.

## Data

The images and OpenGATE outputs are not in git. They are published on Zenodo (CC-BY-4.0):
[doi:10.5281/zenodo.23696023](https://doi.org/10.5281/zenodo.23696023) (all versions). The tests use
version 1.1, [10.5281/zenodo.23702111](https://doi.org/10.5281/zenodo.23702111); version 1.0
([10.5281/zenodo.23696024](https://doi.org/10.5281/zenodo.23696024)) lacks the electron and photon images.

The zip contains the cropped trunk CT (4 mm), the liver and organ label images, the material tables and,
for each radionuclide, the OpenGATE dose rate, its uncertainty and the simulation statistics; for the
VRT simulations (Lu-177, Y-90, I-131), also the electron (`e/`) and photon (`gamma_tle/`) parts, used to
estimate the uncertainty of the organ means (step 6).
Unzip the content of its top folder into `p142/data/`. The full-body CT and the acquisition metadata are not published.
