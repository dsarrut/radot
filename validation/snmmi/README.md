# SNMMI 177Lu Dosimetry Challenge

Results: [docs/8_validation_snmmi.md](../../docs/8_validation_snmmi.md).

Evaluation of radot with the data of the SNMMI 2021 177Lu Dosimetry Challenge (2 patients treated with 177Lu-DOTATATE, SPECT/CT at 4 time points). There is no ground truth: radot is compared with the distribution of the challenge participants' results (Brosch-Lenz et al., J Nucl Med 2023;64:1109, Supplemental Table 1).

| File | Role |
|---|---|
| `cases.py` | Patients, VOIs, data layout, access to the participants' statistics (single source of truth) |
| `challenge_stats.json` | Participants' absorbed doses per task, patient and VOI (mean, quartiles, QCD), parsed from Supplemental Table 1 |
| `prepare_data.py` | Extracts the Deep Blue zips and converts the images and masks to NIfTI with radot |
| `run_gate.py` | Optional: OpenGATE Monte Carlo dose-rate map with the TIA map as source (needs opengate) |
| `t5.py` | Task 5: doses from the provided TIA map and VOIs (`svalue`, `ldm`, and `gate` when the map exists), compared with the participants |
| `test_snmmi.py` | Tests (the T5 test is skipped without the data) |
| `results/t5_radot.json` | radot T5 results |
| `report.py` | Writes `docs/8_validation_snmmi.md` from the stored results |

## Data

From the University of Michigan Deep Blue collection [hm50ts030](https://deepblue.lib.umich.edu/data/collections/hm50ts030), license CC BY-NC 4.0 (research only, no commercial use). Deep Blue does not allow scripted downloads: download by hand into `data/snmmi_challenge/` (or set `RADOT_SNMMI_DATA`):

| Zip | Content | Needed for |
|---|---|---|
| `gb19f607x.zip` | SPECT TIA maps, Bq/mL·s, registered to the SPECT of time point 1 | T5 |
| `4t64gn408.zip` | Organ and lesion VOIs (RTSTRUCT and masks) at each time point | T4, T5 |
| `dn39x173n.zip` | SPECT (Bq/mL) and CT at 4 time points | T4, T5 (CT for the masses) |
| `1j92g763g.zip`, `xp68kg44x.zip`, `js956g027.zip` | planar, projections, pre-therapy images | not used |

Then:

```bash
python validation/snmmi/prepare_data.py    # extract + convert (data/snmmi_challenge/nifti/)
python validation/snmmi/run_gate.py        # optional: OpenGATE maps in data/.../gate/ (~10 min per patient)
python validation/snmmi/t5.py              # T5 table, writes results/t5_radot.json
pytest validation/snmmi
python validation/snmmi/report.py          # regenerate docs/8_validation_snmmi.md
```

Acknowledgement requested by the data owners for publications and presentations:

> We acknowledge Yuni Dewaraja at the University of Michigan for providing access to Lu-177 patient imaging data and segmentations obtained with support from R01EB022075 awarded by NIBIB and R01CA240706 awarded by NCI, NIH. Data were shared via the University of Michigan Deep Blue Data sharing repository.

## Notes on the data

- Patient 4 = challenge patient A (male, 7.21 GBq, days 0/1/4/5, spleen, 2 lesions); patient 6 = B (female, 7.31 GBq, days 0/1/4/8, splenectomised, 4 lesions). The VOI readme sentence "Patient 4 does not have a spleen" is wrong.
- The TIA map, the SPECT and the masks share the same grid (256×256×199, 1.95 mm); the CT (0.98×0.98×3 mm) is resampled to it for the masses (density 1 + HU/1000).
- MIM exports the SPECT and TIA as PT / BQML: `prepare_data.py` sets the sidecars to `spect` / `bqml` and `spect` / `bqsml`.
- The liver VOI contains the liver lesions (A: 14845 of the 15550 voxels of lesion 1). `t5.py` reports the liver as given and without the lesions.

## T5 results

radot `svalue` (TIA × OpenDose self S-value, mass-scaled; male phantom for A, female for B) and `ldm` (local deposition), against the participants' median:

- All organs and lesions are within 5 % of the participants' median (participants' QCD 3–7 %), except the liver of A.
- Liver of A: the participants' spread (QCD 32 %) comes from the liver definition. The liver VOI as given gives 2.39 Gy (median 2.25 Gy); without the lesions, 1.32 Gy (Q1 1.27 Gy).
- Patient A doses are 2–5 % below the median, B within 3 %. LDM is 1–2 % below `svalue` (no photon contribution).
