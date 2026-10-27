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
| `download_data.py` | Downloads the images and OpenGATE outputs from Zenodo into `p142/data/` |
| `p142/data/` | Images and OpenGATE outputs (not in git, downloaded on demand) |

## Run

```bash
pytest validation                            # radot vs the stored OpenGATE values (no download)
pytest validation --opengate=images          # radot on the real images vs the OpenGATE images
python validation/opengate/report.py         # regenerate docs/7_validation.md
```

- Default level: the organ is modelled as a single voxel with the volume, mass and activity of the
  benchmark JSON. It needs no image and runs in CI.
- `images` level: radot runs on the real CT, liver and labels, and the OpenGATE dose rate is averaged
  over the same masks. The data (34 MB) are downloaded on first use (`python validation/opengate/download_data.py`
  does the same) and checked against their SHA-256. It also checks that the stored JSON values are those
  of the OpenGATE images. On p142, both levels agree within 0.02 percentage points.

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

## Data

The images and OpenGATE outputs are not in git. They are published on Zenodo (CC-BY-4.0):
[doi:10.5281/zenodo.23696023](https://doi.org/10.5281/zenodo.23696023) (all versions; version 1.0 is
[10.5281/zenodo.23696024](https://doi.org/10.5281/zenodo.23696024)).

The zip contains the cropped trunk CT (4 mm), the liver and organ label images, the material tables and,
for each radionuclide, the OpenGATE dose rate, its uncertainty and the simulation statistics.
Unzip the content of its top folder into `p142/data/`. The full-body CT and the acquisition metadata are not published.
