# AGENTS.md

Instructions for AI coding agents working on `radot`: setup, architecture, module responsibilities, and the design rules to follow.

## 0. Quick Start

```bash
pip install -e ".[dev]"          # core + pytest/scipy; ".[seg]" adds TotalSegmentator
pytest -q                         # unit tests + validation levels that need no data (CI: Python 3.11–3.14)
```

*   Never commit patient data or anything under `data/`. The SNMMI challenge data (CC BY-NC) and the outputs derived from it are never committed or redistributed.
*   Commits use conventional prefixes (`feat:`, `fix:`, `docs:`, `test:`, ...).
*   Planning and status notes are in `dev_docs/` (gitignored); start with `dev_docs/dev_readme_070_consolidate.md`.

---

## 1. Project Principles & Architecture

`radot` (Radiation Dosimetry for Therapy) provides tools to manage dosimetry calculations, image resampling, unit normalization, and metadata workflows in radionuclide therapy. The codebase adheres to the standard Python Packaging Authority (PyPA) **src-layout** (`src/radot/`).

### Directory Structure & Module Responsibilities
*   **`src/radot/meta.py`**: Reads/writes metadata from/to companion JSON sidecars. Ensures standard unit defaults and provides time-from-injection calculation helper methods.
*   **`src/radot/dicom.py`**: Scans and converts patient DICOM series (CT, SPECT, PET) to continuous formats (NIfTI, NRRD) and generates initial companion JSON sidecars.
*   **`src/radot/segmentation.py`**: Runs TotalSegmentator on CT scans and maps label IDs/naming conventions to standard AAPM TG-263 nomenclature.
*   **`src/radot/resampling.py`**: Performs isotropic resampling or grid-alignment resampling. Handles continuous values vs. discrete labels, and isotope activity preservation.
*   **`src/radot/stats.py`**: Computes volumes, masses (density-calibrated from CT HU), and intensity/peak metrics (e.g. SUVpeak via Gaussian smoothing).
*   **`src/radot/db.py`**: Coordinates radionuclide physical data (using `icrp107-database`), decay chain network trees and cumulative secular equilibrium yields from curated IAEA LiveChart data (`decay_chains.json`), and local caching/lazy-loading of OpenDose reference masses and S-values.
*   **`src/radot/nomenclature.py`**: Encapsulates the 3-stage cascading mapping dictionary logic between segmentations (TotalSegmentator), clinical naming standards (AAPM TG-263), and reference phantoms (OpenDose).
*   **`src/radot/convert.py`**: Handles conversion of image voxel values and files between radioactivity units (`bq`, `bqml`, `suv`, `sul`), using patient weight, injected activity, height, sex, and Lean Body Mass (LBM) calculation methods.
*   **`src/radot/log.py`**: Configures package-level logging levels (`QUIET`, `INFO`, `DEBUG`) and formats standard outputs.
*   **`src/radot/dose/`**: Dosimetry sub-package. The public API is re-exported by `radot/dose/__init__.py`: always import from `radot.dose`, and add any new public name there.
    *   `models.py`: dose formulas, voxel-level Local Deposition Model (LDM, with optional secular-equilibrium daughter cascade) and organ-level Hänscheid 2017/2018, Madsen 2018, S-value (TIA).
    *   `strategies.py`: per-ROI strategies (S-value lookup, mass scaling, log), the `DOSIMETRY_STRATEGIES` registry, `TIA_UNITS` / `TIA_METHODS`, and the loop over labels.
    *   `rate_map.py`: `dose_from_rate_map`, external Monte Carlo dose-rate map.
    *   `mtp.py`: multi-time-point Jackson 2013 (exactly 3 time points), fit and `compute_roi_dose_mtp`.
    *   `uncertainty.py`: ROI dose uncertainty (`roi_dose_uncertainty_estimate`, also used by the validation).
    *   `compute.py`: input checks and `compute_roi_dose` (SimpleITK images, dispatch to the above).
    *   `config.py`: JSON dose configs (`load_dose_config`); `files.py`: `compute_roi_dose_file` (files, resampling, output) and the results table.
*   **`src/radot/data/`**: Packaged JSON data (mapping dictionaries, decay chains) and the OpenDose cache `data/opendose/icrp_110_{male,female}/` (declared in `[tool.setuptools.package-data]`).

### Repository Layout (outside the package)
*   **`tests/`**: Unit tests. They must run without patient data and without TotalSegmentator (tests needing it are skipped when it is not installed). `tests/conftest.py` mocks the OpenDose S-values with a module-scoped fixture (it must not leak into `validation/`).
*   **`validation/opengate/`**: Comparison of radot organ doses with OpenGATE Monte Carlo on one patient CT (p142). Not part of the `radot` package. See section 3.
*   **`docs/`**: mkdocs user documentation (`mkdocs.yml`). `docs/7_validation_opengate.md` and `docs/8_validation_snmmi.md` are generated, never edit them by hand.
*   **`dev_docs/`**: Developer notes, numbered `dev_readme_NNN_<topic>.md` (gitignored, local only). Obsolete notes go to `dev_docs/old/`.
*   **`scripts/`**: Maintenance scripts (`update_decay_chains.py`, `compare_dosimetry.py`).
*   **`.github/workflows/tests.yml`**: CI, `pip install -e ".[dev]"` then `pytest -q` on Python 3.11–3.14.

### Packaging
*   License: LGPL-3.0-or-later (`LICENSE.md`, same as opengate).
*   Core dependencies are light: SimpleITK, pydicom, click, numpy, tqdm, icrp107-database. Do not add heavy runtime dependencies.
*   Extras: `[seg]` = TotalSegmentator (imported lazily in `segmentation.py`, with a clear error if missing), `[dev]` = pytest + scipy, `[docs]` = mkdocs-material.
*   `opengate` is **not** a dependency of radot; only the `--opengate=simulate` validation level needs it.

---

## 2. Critical Development Rules

### Metadata & Unit Safety
1.  **Read from Header/Sidecar First**: Tools and functions requiring modality, units, patient weight, injected activity, height, or sex must first look up these parameters from the image's companion JSON sidecar.
2.  **No Static Defaults in CLI**: Do not hardcode static defaults (e.g., `--type` or `--unit` defaulting to `ct`/`hu`) in click/argparse options, as this disables dynamic metadata resolution. Let them default to `None` and resolve dynamically.
3.  **Explicit Fallbacks**: If required parameters are missing from both options and sidecar/metadata, raise a clear `ValueError` specifying the missing attribute.
4.  **Metadata Inheritance**: Any tool that creates or resamples an image (e.g., `resample`, `convert-units`) must copy, update, or generate a matching sidecar JSON file to accompany the new image.

### Resampling & Unit Conversions
*   **Preserve Total Activity**: When resampling images in absolute activity units (e.g., `bq`), you must scale voxel values by the volume change ratio (new voxel volume / old voxel volume) to keep the total activity constant.
*   **Interpolation Strategy**: Use Linear interpolation for continuous modalities (CT, PET, SPECT, Dose) and Nearest Neighbor interpolation for discrete labels.
*   **Anti-Aliasing**: Apply a Gaussian anti-aliasing filter during downsampling for continuous images (enabled by default); disable it for discrete label maps.
*   **Unit Conversions**: Convert image values between `bq`, `bqml`, `suv`, and `sul` in `convert.py`.
    *   `suv`: Requires patient weight ($W$ in kg) and activity ($A_{\text{inj}}$ in MBq).
    *   `sul`: Normalizes by Lean Body Mass (LBM). Supports analytical methods (`james`, `janmahasatian` requiring height and sex) and patient CT fat segmentation subtraction (`segmentation`).

### Radionuclide Decay Chains & Dosimetry
*   **Data Prioritization**: ICRP 107 is authoritative for single-nuclide emission spectra and non-penetrating energy constants ($\Delta_{\text{local}}$); IAEA LiveChart is authoritative for decay tree topology and cumulative branch yields ($Y_i$).
*   **Cascade Dosimetry**: LDM calculations support multi-generation cascades under secular equilibrium (`cascade=True` / `--cascade`). Total non-penetrating energy is summed over all progeny weighted by cumulative yield $Y_i$.
*   **Dosimetry Model Registration**: Use the strategy registry pattern to route dosimetry methods (`hanscheid2017`, `hanscheid2018`, `madsen2018`, `ldm`; `jackson2013` for multi-time-point; `svalue` for a time-integrated activity image). New models register dynamically.
*   **TIA Images**: Units `bqs` (Bq·s per voxel) and `bqsml` (Bq·s/mL) mark a time-integrated activity image. Only the methods in `TIA_METHODS` (`svalue`, `ldm`) accept it; no time point, Teff, dose rate or total-activity override is used. Results then report `tia_mbq_s` / `tia_mbq_h` instead of `roi_teff_h`.
*   **Dose-rate Map**: `--dose-rate-map` (`dose_from_rate_map` in `dose/rate_map.py`) takes an external dose-rate image (e.g. OpenGATE, Gy/s) computed with the activity image as source, scaled to the image total with the map activity (`--dose-rate-map-activity` or sidecar `source_activity_mbq`). It replaces the S-value / local deposition; the method only selects the time integration (`DOSE_RATE_MAP_METHODS`). radot never runs GATE. Planned engine × time-model refactor: `dev_docs/dev_readme_090_dose_engine_time_model.md`.
*   **Validation**: Any change to a dosimetry model, S-values, masses or decay data must keep `pytest validation` passing, then regenerate `docs/7_validation_opengate.md` and `docs/8_validation_snmmi.md` (see sections 3 and 4).

### ROI Nomenclature & OpenDose Mapping
*   **Master Anchor**: AAPM TG-263 is the primary hub standard.
*   **3-Stage Resolution Cascade**:
    1.  *Stage 1*: Segmented Name $\rightarrow$ TG-263 standard code $\rightarrow$ OpenDose organ name.
    2.  *Stage 2*: Direct TS-to-OD bypass fallback when no standard TG-263 name exists.
    3.  *Stage 3*: Raw retain fallback for custom structures.
*   **Skeletal Aggregation Prohibition**: Do NOT automatically group individual bones (`rib_*`, `vertebra_*`) into a single `"Skeleton"` compartment. OpenDose defines detailed sub-components and lacks a unified `"Skeleton"` S-value. If an unmapped bone is requested, halt with an exception.
*   **Failsafe Rule**:
    *   *Generic Loop*: If looping over all labels in a volume, log warnings and skip unmapped targets.
    *   *Explicit Request*: If the user explicitly requests an unmapped ROI via `--rois`, raise a `click.ClickException` immediately and ask for custom S-value overrides.
*   **Data Decoupling**: Store mapping translation dictionaries in flat JSON files within the `src/radot/data/` folder and load them dynamically. Avoid hardcoded structures.

### Logging & Console Outputs
*   Use Python's built-in `logging` module. Avoid raw `print()` statements in programmatic libraries to prevent standard output pollution.

---

## 3. OpenGATE Validation (`validation/opengate/`)

*   **Single source of truth**: `benchmarks.py` holds the cases, acceptance ranges and computations. `test_opengate.py` and `report.py` both import it. Do not duplicate benchmark values or tolerances in `tests/`.
*   **Case**: patient p142, CT cropped and resampled to 4 mm, uniform activity in the liver (1 MBq at t = 0; I131 0.5 MBq, Ac225 0.1 MBq in the configs; the Lu177, Y90 and I131 references were run with more activity, `reference.json`, and all results are per MBq). Radionuclides Lu177, Y90, I131, Ac225. Self-dose liver ← liver, and cross-doses to 8 organs (`CROSS_TARGETS`; only Kidney_R is asserted). Male OpenDose phantom, S-values scaled to the patient organ mass (mass scaling of cross-doses kept on purpose).
*   **Reference values**: `p142/p142_*_opengate_benchmark.json`, produced from the OpenGATE outputs by `p142/create_benchmark_json.py` (with the opengate version/commit). OpenGATE configs: `p142/p142_dose_rate_*.json`; input preparation: `p142/prepare/`.
*   **Test levels** (option defined in `validation/conftest.py`; `testpaths = ["tests", "validation"]`):
    *   `pytest` (default, CI): single-voxel proxy of each organ (benchmark volume, mass, activity) vs the stored GATE values. Needs no images.
    *   `pytest --opengate=images`: radot on the real images vs the GATE dose images. Data downloaded from Zenodo (record 23809698, v1.2) into `p142/data/` (gitignored) by `download_data.py`, checked by SHA-256.
    *   `pytest --opengate=simulate`: reruns the OpenGATE simulations (needs opengate, hours of CPU) into `p142/data/simulated/`, checks them against the stored values within 3 combined standard uncertainties, then compares radot.
*   **Uncertainty**: GATE organ-mean uncertainty uses a per-part model (electrons: mean β range; photons: mean organ chord 4V/S). Checked with 20 Lu177 batches (`batch_uncertainty.py`, `p142/p142_lu177_batch_uncertainty.json`).
*   **Report**: `python validation/opengate/report.py` writes `docs/7_validation_opengate.md`; a test fails when the page is out of date.
*   **Trap**: under pytest, `import opengate` resolves to the `validation/opengate/` folder (namespace package). Never `import opengate` directly in validation code; the opengate version is read from the install record.
*   **Bremsstrahlung**: the references need opengate with the `bremsstrahlung` option of `opengate_dose_rate` (OpenGATE PR #1223, on by default); without it the Y90 cross-doses are 10 to 140 times too low (for Lu177 and I131 it is below 1 %). `p142/create_campaigns.py` refuses an opengate without it.
*   **Cross-doses**: large differences except heart/stomach come from the anatomy (phantom vs patient geometry, organ distances, mass scaling); shown, only Kidney_R asserted. Open point: OpenGATE issue #1215 (default material tables), then drop the tables from the configs and the zip. Details in `dev_docs/dev_readme_080_opengate_validation.md` (done, kept until PR #5 is merged).

## 4. SNMMI Challenge Evaluation (`validation/snmmi/`)

*   **Data**: SNMMI 177Lu Dosimetry Challenge (Deep Blue, CC BY-NC 4.0, research only). Downloaded by hand into `data/snmmi_challenge/` (gitignored, or `$RADOT_SNMMI_DATA`), converted by `prepare_data.py`. Never commit or redistribute the images.
*   **Reference**: no ground truth; radot is compared with the participants' statistics (`challenge_stats.json`, parsed from Brosch-Lenz 2023 Supplemental Table 1). Patient 4 = A, patient 6 = B.
*   **Single source of truth**: `cases.py` (patients, VOIs, paths, statistics access). `t5.py` computes task 5 (TIA map + VOIs → doses; adds `gate` when `run_gate.py` has produced the OpenGATE map, used through `--dose-rate-map`). `test_snmmi.py`: the statistics and page tests run in CI, the T5 tests are skipped without the data.
*   **Report**: `python validation/snmmi/report.py` writes `docs/8_validation_snmmi.md` from the stored `results/t5_radot.json` (no data needed); a test fails when the page is out of date. After rerunning `t5.py`, regenerate the page.
*   **Status and findings**: `dev_docs/dev_readme_082_snmmi_challenge.md`.

---

## 5. Reference Documentation Map

Refer to these Markdown documents in the `docs/` directory for detailed specifications and user guides:
*   [quickstart.md](docs/quickstart.md): End-to-end CLI usage and Python API workflows.
*   [6a_dose.md](docs/6a_dose.md): Mathematical formulations for LDM, Hänscheid, and Madsen dosimetry models.
*   [3a_images_units.md](docs/3a_images_units.md): Unit conversions (`bq`, `bqml`, `suv`, `sul`) and LBM strategies.
*   [3b_sul.md](docs/3b_sul.md): Stoichiometric CT density integration for SUL without weight dependence.
*   [5a_radionuclide.md](docs/5a_radionuclide.md): ICRP-107 and IAEA LiveChart decay cascade data integration and prioritization.
*   [1a_dicom.md](docs/1a_dicom.md): DICOM conversion, series sorting, and sidecar metadata extraction.
*   [1b_meta.md](docs/1b_meta.md): Companion JSON sidecar specifications and schema.
*   [2a_segmentation.md](docs/2a_segmentation.md): TotalSegmentator automation and TG-263 label mapping.
*   [4_images_resample.md](docs/4_images_resample.md): Grid resampling and total activity preservation.
*   [2b_roi_names.md](docs/2b_roi_names.md): Comprehensive ROI mapping policy and cascade stages.
*   [6b_stats.md](docs/6b_stats.md): Volume, mass, and statistical dosimetry reporting.
*   [5b_opendose.md](docs/5b_opendose.md): OpenDose S-values and phantom caching.
*   [3c_peak.md](docs/3c_peak.md): SUVpeak equations, spherical volume integration, and Gaussian filtering.
*   [logging.md](docs/logging.md): Configuration and levels of the package logging system.
*   [python.md](docs/python.md): Python API reference and integration recipes.
*   [7_validation_opengate.md](docs/7_validation_opengate.md): Generated comparison with OpenGATE (do not edit).
*   [8_validation_snmmi.md](docs/8_validation_snmmi.md): Generated comparison with the SNMMI challenge participants (do not edit).
*   [validation/opengate/README.md](validation/opengate/README.md): How the OpenGATE reference was produced and how to run each validation level.
*   [validation/snmmi/README.md](validation/snmmi/README.md): SNMMI challenge data, how to prepare them and run T5.
*   [dev_docs/dev_readme_010_checklist.md](dev_docs/dev_readme_010_checklist.md): Developer parameter and metadata validation checklist.
*   [dev_docs/dev_readme_020_db_api.md](dev_docs/dev_readme_020_db_api.md): Developer reference for database and phantom layers.
*   [dev_docs/dev_readme_070_consolidate.md](dev_docs/dev_readme_070_consolidate.md): Status of the work towards v0.1 (start here for planning).
*   [dev_docs/dev_readme_080_opengate_validation.md](dev_docs/dev_readme_080_opengate_validation.md): OpenGATE validation status, decisions and open points.
*   [dev_docs/dev_readme_090_dose_engine_time_model.md](dev_docs/dev_readme_090_dose_engine_time_model.md): `--dose-rate-map` (external Monte Carlo dose rate) and the planned dose-rate engine × time model refactor.
*   [dev_docs/dev_readme_082_snmmi_challenge.md](dev_docs/dev_readme_082_snmmi_challenge.md): Plan and status (T5 done) to evaluate radot with the SNMMI 177Lu Dosimetry Challenge data.
*   [dev_docs/dev_readme_085_lesion_validation.md](dev_docs/dev_readme_085_lesion_validation.md): Plan (not started) for synthetic lesions in the OpenGATE validation and a sphere S-value model.

