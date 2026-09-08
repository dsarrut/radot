# RADOT (Radiation Dosimetry for Therapy)

The radot library provides a unified Python API and command-line interface to manage dosimetry calculations, image resampling, and metadata workflows in radionuclide therapy.

## Installation

The package requires Python 3.11 or higher

To install the package in editable mode for local development:
```bash
pip install -e .
```

To enable command autocompletion for the radot CLI in zsh, append the following line to your shell configuration profile (e.g., ~/.zshrc) or at the end of your "activate" script from the venv:
```bash
eval "$(_RADOT_COMPLETE=zsh_source radot)"
```

## Quick start

Read the [quickstart](docs/quickstart.md) and start from there.

## Core Modules

1. DICOM Conversion
   - Scans and converts patient DICOM series (CT, SPECT, PET) to dosimetry-ready formats (NIfTI, NRRD).
   - Generates companion JSON sidecars to store key acquisition and therapy metadata such as injected activity or date of injection/acquisition. 

2. CT Segmentation
   - Automates anatomical segmentation using TotalSegmentator.
   - Remaps label indices and names dynamically to standardized AAPM TG-263 nomenclature.

3. PET (basic) Segmentation **TODO LATER**
   - PET normalisation with Gafita's method (SUL / Liver)
   - TTV ? etc

4. Image Resampling
   - Resamples multi-modality grids to target isotropic spacings or reference image grids.
   - Preserves total activity concentration metrics during spatial scaling.

5. Radionuclide Database
   - Queries physical half-lives and decay emission yields/energies from the reference ICRP-107 database.
   - Query/store Opendose S-values, mass data etc. for each nuclide

6. Dose Computation
   - Implements voxel-mean absorbed dose computation models (e.g. Madsen 2018, Hanscheid 2017/2018).
   - Local Deposition Model (LDM) voxel-level modeling

