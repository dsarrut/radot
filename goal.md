# radot (Radiation Dosimetry for Therapy)

A Python library and command-line tools to manage dosimetry in radionuclide therapy.

## Core Stack & Design
- **Packaging**: Standard `pyproject.toml` with `setuptools`.
- **Target OS**: Linux, macOS, and Windows.
- **Key Libraries**: `SimpleITK`, `pydicom`, `click`.
- **Interface**: Command-line interface (CLI) only, no GUI.

## Proposed Modules (Ideas/Drafts)

The following modules represent the current planned architecture and may evolve over time:

1. **DICOM Conversion**:
   - Convert DICOM images (CT, SPECT, PET, RTSTRUCT, etc.) to NIfTI or similar image formats.

2. **CT Segmentation & Label Grouping**:
   - Run automatic segmentations using `TotalSegmentator`.
   - Organize and rename/group labels according to **AAPM TG-263** nomenclature standards.

3. **Image Resampling**:
   - Resample CT, labels, and functional images (SPECT, PET).
   - Ensure preservation of activity concentration/physical units during resampling.

4. **Dose Computations (Main Module)**:
   - Implement dosimetry methods from literature (e.g., Madsen, Jackson, S-values, point kernels).
   - Compute dose rate maps and time-integration over time-activity curves, including support for GATE outputs.

5. **TTV (Total Tumor Volume) Computation**:
   - Compute Total Tumor Volume (TTV) / Metabolic Tumor Volume (MTV) based on segmentation labels and SUV/activity concentration thresholds.
