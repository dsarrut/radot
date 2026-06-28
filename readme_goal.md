# radot (Radiation Dosimetry for Therapy)

The radot library provides a unified Python API and command-line interface to manage dosimetry calculations, image resampling, and metadata workflows in radionuclide therapy.

## Core Modules

1. DICOM Conversion:
   - Scans and converts patient DICOM series (CT, SPECT, PET) to dosimetry-ready formats (NIfTI, NRRD).
   - Generates companion JSON sidecars to store key acquisition and therapy metadata.

2. CT Segmentation:
   - Automates anatomical segmentation using TotalSegmentator.
   - Remaps label indices and names dynamically to standardized AAPM TG-263 nomenclature.

**TODO LATER -> PET/SPECT Segmentation**

3. Image Resampling:
   - Resamples multi-modality grids to target isotropic spacings or reference image grids.
   - Preserves total activity concentration metrics during spatial scaling.

4. Dose Computation:
   - Implements voxel-mean absorbed dose computation models (e.g. Madsen 2018, Hanscheid 2017/2018).
   - Integrates lazy on-demand downloads of OpenDose reference masses and S-values.

5. Radionuclide Database:
   - Queries physical half-lives and decay emission yields/energies from the reference ICRP-107 database.

## Installation

The package requires Python 3.10 or higher

To install the package in editable mode for local development:
```bash
pip install -e .
```

To enable command autocompletion for the radot CLI in zsh, append the following line to your shell configuration profile (e.g., ~/.zshrc) or at the end of your "activate" script from the venv:
```bash
eval "$(_RADOT_COMPLETE=zsh_source radot)"
```
