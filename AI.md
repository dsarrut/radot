# AI Developer Guide (AI.md)

This document provides a high-level overview of the `radot` library architecture, core design rules, metadata guidelines, and module responsibilities to help future AI coding assistants develop and maintain this repository.

---

## 1. Project Principles & Architecture

`radot` (Radiation Dosimetry for Therapy) provides tools to manage dosimetry calculations, image resampling, and metadata workflows in radionuclide therapy.

### Directory Structure & Module Responsibilities
*   **`radot/meta.py`**: Reads/writes metadata from/to companion JSON sidecars. Ensures standard unit defaults and provides time-from-injection calculation helper methods.
*   **`radot/dicom.py`**: Scans and converts patient DICOM series (CT, SPECT, PET) to continuous formats (NIfTI, NRRD) and generates initial companion JSON sidecars.
*   **`radot/segmentation.py`**: Runs TotalSegmentator on CT scans and maps label IDs/naming conventions to standard AAPM TG-263 nomenclature.
*   **`radot/resampling.py`**: Performs isotropic resampling or grid-alignment resampling. Handles continuous values vs. discrete labels, and isotope activity preservation.
*   **`radot/stats.py`**: Computes volumes, masses (density-calibrated from CT HU), and intensity/peak metrics (e.g. SUVpeak via Gaussian smoothing).
*   **`radot/db.py`**: Coordinates radionuclide physical data (using `icrp107-database`) and local caching/lazy-loading of OpenDose reference masses and S-values.
*   **`radot/nomenclature.py`**: Encapsulates the 3-stage cascading mapping dictionary logic between segmentations (TotalSegmentator), clinical naming standards (AAPM TG-263), and reference phantoms (OpenDose).
*   **`radot/dose.py`**: Implements voxel-level (Local Deposition Model / LDM) and organ-level (Hänscheid 2017/2018) dosimetry models using a Strategy registry pattern.

---

## 2. Critical Development Rules

### Metadata & Unit Safety
1.  **Read from Header/Sidecar First**: Tools and functions requiring modality, units, patient weight, or injected activity must first look up these parameters from the image's companion JSON sidecar.
2.  **No Static Defaults in CLI**: Do not hardcode static defaults (e.g., `--type` or `--unit` defaulting to `ct`/`hu`) in click/argparse options, as this disables dynamic metadata resolution. Let them default to `None` and resolve dynamically.
3.  **Explicit Fallbacks**: If required parameters are missing from both options and sidecar/metadata, raise a clear `ValueError` specifying the missing attribute.
4.  **Metadata Inheritance**: Any tool that creates or resamples an image (e.g., `resample`) must copy, update, or generate a matching sidecar JSON file to accompany the new image.

### Resampling & Unit Conversions
*   **Preserve Total Activity**: When resampling images in absolute activity units (e.g., `bq`), you must scale voxel values by the volume change ratio (new voxel volume / old voxel volume) to keep the total activity constant.
*   **Interpolation Strategy**: Use Linear interpolation for continuous modalities (CT, PET, SPECT, Dose) and Nearest Neighbor interpolation for discrete labels.
*   **Anti-Aliasing**: Apply a Gaussian anti-aliasing filter during downsampling for continuous images (enabled by default); disable it for discrete label maps.

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
*   **Data Decoupling**: Store mapping translation dictionaries in flat JSON files within the `radot/data/` folder and load them dynamically. Avoid hardcoded structures.

### Dosimetry Model Registration
*   **Strategy Pattern**: Use the strategy registry pattern to route dosimetry methods (`hanscheid2017`, `hanscheid2018`, `ldm`). New models should register themselves dynamically.

### Logging & console outputs
*   Use Python's built-in `logging` module. Avoid raw `print()` statements in programmatic libraries to prevent standard output pollution.

---

## 3. Reference Documentation Map

Refer to these Markdown documents in the `docs/` directory for detailed requirements:
*   [quickstart.md](file:///Users/dsarrut/src/py/radot/docs/quickstart.md): End-to-end CLI usage and python API workflows.
*   [dev_readme_checklist.md](file:///Users/dsarrut/src/py/radot/docs/dev_readme_checklist.md): Core checklist for parameter and metadata handling.
*   [readme_roi_names.md](file:///Users/dsarrut/src/py/radot/docs/readme_roi_names.md): Comprehensive ROI mapping policy and cascade stages.
*   [dev_readme_db_api.md](file:///Users/dsarrut/src/py/radot/docs/dev_readme_db_api.md): Guidelines for database, isotope layers, and phantoms.
*   [readme_peak.md](file:///Users/dsarrut/src/py/radot/docs/readme_peak.md): Rationale and equations for SUVpeak and Gaussian sigma.
*   [readme_logging.md](file:///Users/dsarrut/src/py/radot/docs/readme_logging.md): Configuration and levels of the package logging system.
