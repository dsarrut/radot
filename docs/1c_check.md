# Input Checks (`radot check`)

`radot check` runs quick quality checks on the inputs of a dosimetry before computing doses: metadata, units, activity balance, ROIs, and later geometry. Wrong doses come more often from the data (units, injection time, weight, misregistration) than from the dose models.

The files are only read, never modified. A missing or suspicious value is reported; fix the sidecar with `radot meta` (or by editing the JSON) and run the check again.

## Usage

```bash
radot check spect.nii.gz                                   # metadata of one activity image
radot check spect.nii.gz --ct ct.nii.gz -m madsen2018      # fields required by the planned method
radot check spect_d1.nii.gz spect_d4.nii.gz spect_d7.nii.gz -m jackson2013   # multi time points
radot check spect.nii.gz -a                                # show all findings, not only the problems
radot check spect.nii.gz -o qc.json --strict               # JSON report, exit code 1 on warnings too
radot check --list                                         # all checks
```

Options:

- `ACTIVITY_PATHS`: activity images (several = multi time points); `--ct`, `-l/--labels`, `--dose`: the other images. Only the checks whose inputs are given run.
- `-m/--method`: the planned dosimetry method; sets which sidecar fields are required (e.g. the dates for `madsen2018`, not for a TIA image with `svalue`).
- `--only`, `--skip`: checks or categories (`meta`, `meta.patient`), comma-separated or repeated.
- `--thresholds file.json`: thresholds merged with the default ones (see below).
- `-a/--all`: also print the OK and INFO findings (measured values, skipped checks).
- `-o report.json`: all the findings, whatever their level.
- Exit code: 1 if a check FAILs (`--strict`: also on WARN), else 0.

Example:

```
$ radot check spect.nii.gz -m madsen2018
  meta      FAIL  meta.required              spect.nii.gz: injection_date missing (time from injection, for madsen2018)
  meta      WARN  meta.activity              spect.nii.gz: injected activity 7.4 MBq outside [500, 12000] for lu177 (MBq expected, not GBq or Bq)

  5 findings: 1 warnings, 1 failures

$ radot meta spect.nii.gz --injection-date "2026-01-10 10:00:00" --activity 7400
$ radot meta spect.nii.gz --rad Lu177       # a missing radionuclide
```

## Levels

| Level | Meaning |
|-------|---------|
| FAIL | The dose would be wrong or cannot be computed (missing unit or dates, acquisition before injection, unknown radionuclide). |
| WARN | Suspicious value, to look at (weight, injected activity, time from injection). |
| INFO | Value reported, or optional field missing. |
| OK | Check passed, with the measured value. |

## Checks

| Code | Check |
|------|-------|
| `meta.required` | Sidecar present; fields needed by the method (unit, radionuclide, dates; weight and activity for `suv`/`sul`; height and sex for `sul`); TIA images only with `svalue`/`ldm`. |
| `meta.modality` | Modality and unit consistent with the role of each image (activity, CT in `hu`, dose in `gy`). |
| `meta.dates_order` | Valid dates (`YYYY-MM-DD HH:MM:SS`), injection before acquisition. |
| `meta.time_from_injection` | Time from injection below 10 half-lives; time points distinct. |
| `meta.patient` | Weight, height, BMI and sex in plausible ranges. |
| `meta.activity` | Injected activity in the usual range of the radionuclide, in MBq (catches GBq or Bq typos). |
| `meta.radionuclide` | Radionuclide known to radot (ICRP 107). |
| `meta.consistency` | Multi time points: same radionuclide, injection date, activity and patient data in all the images. |
| `units.ct_hu` | CT in HU: air near −1000, soft tissue (or fat) peak within ±150 HU, max below 4000. Catches attenuation maps and the +1024 offset; padding below −1100 is ignored. |
| `units.activity_values` | Activity image not empty, < 1 % negative voxels; integer values with a max ≤ 65535 in `bq`/`bqml` are probably counts. |
| `units.suv_range` | SUV/SUL images: max in [0.5, 150]. |
| `units.dose_range` | Dose image in Gy: positive, max below 1000 Gy, no negative values. |
| `activity.fov_fraction` | Total activity in the image / injected activity decayed to the acquisition (not decayed if the image is decay-corrected to the injection). FAIL above 110 % (units, calibration) or below 2 % (counts); WARN below 10 %. SUV images: Σ SUV·mL / weight. TIA images: compared with the injected activity × physical mean life. |
| `activity.mtp_trend` | Multi time points: the total activity decreases at least as fast as the physical decay (20 % tolerance); else different calibrations or decay-corrected images. Reports the effective half-life. |
| `roi.names` | Labels sidecar present (label id → ROI name); ROIs without an OpenDose organ listed (no S-value). |
| `roi.empty` | Labels image not empty; ROIs of the sidecar empty, smaller than 1 mL or touching the image border (truncated) listed. On a cropped image these are normal, so they are reported as INFO. |
| `roi.laterality` | Left/right pairs (`Kidney_L`/`Kidney_R`, `kidney_left`/`kidney_right`) on the correct side in patient coordinates (LPS: patient left = +x); liver and gallbladder on the patient right, spleen and stomach on the left (catches mirrored images). |
| `roi.mass` | With a CT: mass (density from HU, as in `radot stats`) of whole organs (liver, spleen, kidneys, pancreas, thyroid, prostate, brain) vs the OpenDose reference mass for the patient sex (male if unknown). Lobes, vessels and muscles are parts of an OpenDose organ and are not compared. Truncated organs are not compared. |
| `roi.density` | With a CT: mean HU of lungs, bones and soft-tissue organs in a plausible range; many failures mean labels misaligned with the CT. |

## Thresholds

The default thresholds are in `radot/data/check_thresholds.json`, one section per topic (`patient`, `time_from_injection`, `injected_activity_mbq` per radionuclide, `ct`, `activity_image`, `suv`, `dose`, `activity_balance`, `roi`). They are heuristics, generous on purpose. A user file with the same sections overrides the values it contains, e.g. for a dosimetric study with low activities:

```json
{"injected_activity_mbq": {"lu177": [100, 12000]}}
```

```bash
radot check spect.nii.gz --thresholds my_thresholds.json
```

## Python API

```python
from radot.check import run_checks, worst_level, Level

findings = run_checks(activity=["spect_d1.nii.gz", "spect_d4.nii.gz"], ct="ct.nii.gz", method="jackson2013")
for f in findings:
    if f.level >= Level.WARN:
        print(f.level.name, f.code, f.subject, f.message, f.value, f.expected)
print(worst_level(findings).name)
```
