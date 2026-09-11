# AI Developer Guide (AI.md)

This document provides a high-level overview of the `radot` library architecture, core design rules, metadata guidelines, and module responsibilities to help future AI coding assistants develop and maintain this repository.

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
*   **`src/radot/dose.py`**: Implements voxel-level (Local Deposition Model / LDM with optional secular equilibrium daughter cascade support) and organ-level (Hänscheid 2017/2018, Madsen 2018) dosimetry models using a Strategy registry pattern.

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
*   **Dosimetry Model Registration**: Use the strategy registry pattern to route dosimetry methods (`hanscheid2017`, `hanscheid2018`, `madsen2018`, `ldm`). New models register dynamically.

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

## 3. Reference Documentation Map

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
*   [dev_docs/dev_readme_checklist.md](dev_docs/dev_readme_checklist.md): Developer parameter and metadata validation checklist.
*   [dev_docs/dev_readme_db_api.md](dev_docs/dev_readme_db_api.md): Developer reference for database and phantom layers.

