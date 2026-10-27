# RADOT (Radiation Dosimetry for Therapy)

`radot` is a Python library and command-line tool suite for radiation dosimetry, image processing, and metadata workflows in targeted radionuclide therapy.

## Installation

```bash
pip install -e .
```

CT segmentation (`radot segment`) needs TotalSegmentator (and PyTorch), which is an optional extra:
```bash
pip install -e ".[seg]"
```

For development (tests): `pip install -e ".[dev]"`.

To enable shell autocompletion for the CLI (zsh):
```bash
eval "$(_RADOT_COMPLETE=zsh_source radot)"
```

## Quick Start

`radot` can be used as a command-line tool (CLI) or as a Python library.

### 1. Command-Line Interface (CLI)

Run `radot` commands directly from your terminal:

```bash
radot dose spect.nii.gz -m ldm --isotope Lu177 --ct ct.nii.gz -o dose.nii.gz
```

### 2. Python Library

Import `radot` and call dosimetry functions directly in your Python code:

```python
from pathlib import Path
from radot.dose import compute_roi_dose_file

compute_roi_dose_file(
    image_path=Path("spect.nii.gz"),
    labels_path=Path("labels.nii.gz"),
    method="madsen2018",
    isotope="Lu177",
    ct_path=Path("ct.nii.gz"),
    output_path=Path("dose.nii.gz"),
    dose_rate=True,
)
```

For more examples and workflows, see the [Quickstart Guide](quickstart.md).

## Core Modules

1. **DICOM Conversion**: Converts CT, SPECT, and PET series to NIfTI/NRRD with structured JSON metadata sidecars ([Documentation](1a_dicom.md)).
2. **Segmentation**: Automates CT segmentation via TotalSegmentator and remaps labels to standard AAPM TG-263 nomenclature ([Documentation](2a_segmentation.md)).
3. **Unit Conversion**: Converts image voxel values between `bq`, `bqml`, `suv`, and `sul` ([Documentation](3a_images_units.md)).
4. **Image Resampling**: Isotropic and grid-matching resampling while strictly preserving total radioactivity ([Documentation](4_images_resample.md)).
5. **Radionuclide Physics Data**: Queries ICRP Publication 107 emissions spectra, IAEA LiveChart decay cascades, and OpenDose S-values ([Documentation](5a_radionuclide.md)).
6. **Dosimetry Models**: Implements Local Deposition Model (LDM, voxel-level with daughter cascade support), Hänscheid (2017/2018), Madsen (2018), and Jackson (2013) tri-exponential MTP models ([Documentation](6a_dose.md)).

