# OpenGATE validation

Comparison of radot organ doses with OpenGATE Monte Carlo simulations on a patient CT (p142).
The results are published in [docs/7_validation_opengate.md](../../docs/7_validation_opengate.md).

| File | Role |
|---|---|
| `benchmarks.py` | Cases, acceptance ranges and computations (single source of truth) |
| `test_opengate.py` | The same comparisons as tests (collected by `pytest`) |
| `report.py` | Writes `docs/7_validation_opengate.md` |
| `p142/prepare/` | Preparation of the simulation inputs |
| `p142/p142_dose_rate_*.json` | OpenGATE configuration, one per radionuclide |
| `p142/run_references.zsh`, `p142/create_campaigns.py` | Run the reference simulations with `opengate_jobs` (bremsstrahlung included) |
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
python validation/opengate/report.py         # regenerate docs/7_validation_opengate.md
```

- Default level: the organ is modelled as a single voxel with the volume, mass and activity of the
  benchmark JSON. It needs no image and runs in CI.
- `images` level: radot runs on the real CT, liver and labels, and the OpenGATE dose rate is averaged
  over the same masks. The data (74 MB) are downloaded on first use (`python validation/opengate/download_data.py`
  does the same) and checked against their SHA-256. It also checks that the stored JSON values are those
  of the OpenGATE images. On p142, both levels agree within 0.02 percentage points.
- `simulate` level: needs opengate. Reruns the 4 simulations on the downloaded inputs, at the activity
  of the configurations (1 MBq; 0.5 MBq for I-131, 0.1 MBq for Ac-225): Lu-177, Y-90 and I-131 as the
  references (`create_campaigns.py`, parallel jobs on all cores, bremsstrahlung included), Ac-225 with
  `opengate_dose_rate`. Outputs in `p142/data/simulated/` (the reference outputs are not modified). Then
  checks that the new OpenGATE values agree with the stored ones within 3 combined standard uncertainties
  (per MBq), and runs the radot comparisons against the new outputs. Each output folder also contains the
  OpenGATE log, `opengate_version.json` (version and git commit) and, for the VRT runs, `reference.json`
  (simulated activity), which `create_benchmark_json.py` copies into the benchmark JSON.

`pytest` fails if `docs/7_validation_opengate.md` is not up to date.

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
   The crop (34 cm along the body axis, the full 50 cm in the transverse plane) removes the head. A
   Lu-177 simulation on the full-body CT changes the organ dose rates by less than 10 % (left kidney
   +9 %, heart +3 %, liver +0.2 %, the other organs within the statistical uncertainty).
4. One OpenGATE simulation per radionuclide (uniform activity in the liver, dose rate at t = 0; the
   results are given per MBq of the simulated activity). Lu-177, Y-90 and I-131 use variance reduction
   (`"mode": "vrt"`: an electron run and a photon run with the track-length estimator, TLE); they are
   run with `opengate_jobs` in parallel jobs, with more activity than the configurations (Y-90 30 MBq,
   Lu-177 10 MBq, I-131 5 MBq, about 35 min with 12 jobs):
   ```bash
   ../run_references.zsh 12        # outputs in new/output_<radionuclide>_vrt/, then moved to output_<radionuclide>_vrt/
   opengate_dose_rate ../p142_dose_rate_ac225.json        # Ac-225: analog mode, 0.1 MBq
   ```
   `run_references.zsh` calls `create_campaigns.py`, which builds the campaigns from the configurations
   and writes `reference.json` (simulated activity, opengate version) next to the output.
   Bremsstrahlung: `opengate_dose_rate` sets a 1 m production cut for all particles in the electron run,
   so no bremsstrahlung photon is created and its energy stays in the liver; for Y-90 (no γ line) the
   cross-doses were then 10 to 140 times too low. `create_campaigns.py` uses a 2 mm photon cut in the
   electron run (electrons still 1 m) and a TLE dose actor (electrons analog, bremsstrahlung photons with
   TLE). The Ac-225 analog run already has a 2 mm cut. To remove once `opengate_dose_rate` does it.
5. Mean dose rate and uncertainty in the liver (self-dose) and in the organs of `CROSS_TARGETS` in
   `benchmarks.py` (cross-doses: kidneys, spleen, pancreas, adrenals, stomach, heart; whole organs entirely
   inside the cropped CT), with their mass and closest distance to the liver:
   ```bash
   python ../create_benchmark_json.py
   ```
   It also crops `labels_all_4mm.nii.gz` to `labels_all_4mm_cropped.nii.gz` (slices 282–366) if needed.
6. Uncertainty of the organ means. OpenGATE gives a history-by-history uncertainty per voxel, not for
   an organ, and one history deposits in several voxels, which are therefore correlated.
   `benchmarks.gate_roi_uncertainty` estimates the organ uncertainty from the voxel maps, separately for
   the two VRT simulations: electrons (`e/`, a cluster of the mean β range) and photons (`gamma_tle/`,
   the mean chord 4V/S of the organ, since a photon deposits along its path). The electron run also
   contains the bremsstrahlung photons: in organs beyond the electron reach (more than 10 mm from the
   liver), its dose is bremsstrahlung only and the photon model is used for it too. This estimate was
   checked with the batch method: the Lu-177 simulation is repeated 20 times with different random
   seeds, each at 1 MBq, and the relative standard deviation of the 20 organ means is the
   uncertainty of one simulation (about 1 min per batch on 8 threads):
   ```bash
   python ../../batch_uncertainty.py --run   # batches in batches/run_NN/, result in p142_lu177_batch_uncertainty.json
   ```
   | Lu-177 (20 batches) | Measured | Estimated |
   |-----------|-----------------|-----------|
   | Liver     | 0.13 % ± 0.02 | 0.13 % |
   | Kidney_R  | 1.06 % ± 0.17 | 0.94 % |
   | Kidney_L  | 3.06 % ± 0.50 | 2.52 % |
   | Spleen    | 2.45 % ± 0.40 | 2.07 % |
   | Pancreas  | 1.82 % ± 0.30 | 1.31 % |
   | Adrenal_R | 2.83 % ± 0.46 | 2.99 % |
   | Adrenal_L | 4.60 % ± 0.75 | 4.11 % |
   | Stomach   | 0.90 % ± 0.15 | 1.15 % |
   | Heart     | 0.97 % ± 0.16 | 0.67 % |

   (These batches predate the bremsstrahlung fix, negligible for Lu-177; the test scales the reference
   estimate, now at 10 MBq, to 1 MBq.) The estimate is within 2 standard deviations of the measurement for all organs (0.7 to 1.3 times the
   measured value); `pytest` checks it. The former estimate (cluster of the maximum β range on the merged
   map) gave 0.20 % for the liver and 0.58 % for the right kidney. The mean of the 20 batches agrees with
   the reference simulation within 0.7 % per MBq.

   Y-90 (bremsstrahlung only outside the liver), 10 electron runs of 1 MBq with the bremsstrahlung of
   `create_campaigns.py` (run by hand, not with `batch_uncertainty.py`): the estimate is 1.1 times the
   measured uncertainty for the right kidney and 1.0 for the spleen, 0.65 for the pancreas and the left
   kidney, and too high (1.4 to 2.2) for the liver and the organs touching it, with a ±24 % precision
   on the measurement. The mean of the 10 runs agrees with the reference (right kidney +3.2 % ± 1.1,
   pancreas +0.1 % ± 2.6), so the references do not depend on the machine (laptop vs cluster).
   The 3-sigma check of the `simulate` level can fail by chance on an organ where the estimate is low
   (a first rerun gave −10 % for the right kidney, 3.0 sigma; the next one −2 %).

## Data

The images and OpenGATE outputs are not in git. They are published on Zenodo (CC-BY-4.0):
[doi:10.5281/zenodo.23696023](https://doi.org/10.5281/zenodo.23696023) (all versions). The download uses
version 1.1, [10.5281/zenodo.23702111](https://doi.org/10.5281/zenodo.23702111); version 1.0
([10.5281/zenodo.23696024](https://doi.org/10.5281/zenodo.23696024)) lacks the electron and photon images.

Not yet published: the Lu-177, Y-90 and I-131 references with bremsstrahlung (step 4), on which the
benchmark JSONs are now based. Until they are on Zenodo (version 1.2), the `images` level needs them in
`p142/data/output_<radionuclide>_vrt/`; with the version 1.1 outputs, the check that the JSON values are
those of the images fails for these three.

The zip contains the cropped trunk CT (4 mm), the liver and organ label images, the material tables and,
for each radionuclide, the OpenGATE dose rate, its uncertainty and the simulation statistics; for the
VRT simulations (Lu-177, Y-90, I-131), also the electron (`e/`) and photon (`gamma_tle/`) parts, used to
estimate the uncertainty of the organ means (step 6).
Unzip the content of its top folder into `p142/data/`. The full-body CT and the acquisition metadata are not published.
