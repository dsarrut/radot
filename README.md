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

For more examples and workflows, see the [Quickstart Guide](docs/quickstart.md).

## Core Features

1. **DICOM Conversion**: Converts CT, SPECT, and PET series to NIfTI/NRRD with structured JSON metadata sidecars.

2. **Segmentation**: Automates CT segmentation via TotalSegmentator and remaps labels to standard AAPM TG-263 nomenclature.

3. **Unit Conversion**: Converts image voxel values between `bq`, `bqml`, `suv`, and `sul`

4. **Image Resampling**: Isotropic and grid-matching resampling while strictly preserving total radioactivity.

5. **Radionuclide Physics Data**: Queries ICRP Publication 107 emissions spectra, IAEA LiveChart decay cascades, and OpenDose S-values.

6. **Dosimetry Models**: Implements Local Deposition Model (LDM, voxel-level with daughter cascade support), Hänscheid (2017/2018), and Madsen (2018) models.

## Documentation

- [Quickstart](docs/quickstart.md)

- 1. [DICOM Conversion](docs/1a_dicom.md) & [Metadata Sidecars](docs/1b_meta.md)

- 2. [CT Segmentation & TG-263](docs/2a_segmentation.md) | [ROI Names Policy](docs/2b_roi_names.md)

- 3. [Units & SUL Normalization](docs/3a_images_units.md) | [SUL CT Integration](docs/3b_sul.md) | [SUVpeak](docs/3c_peak.md)

- 4. [Image Resampling](docs/4_images_resample.md)

- 5. [Radionuclides & Decay Chains](docs/5a_radionuclide.md) | [OpenDose Database](docs/5b_opendose.md)

- 6. [Dosimetry Models](docs/6a_dose.md) | [Statistics & Reporting](docs/6b_stats.md) 

- [Python API Reference](docs/python.md) | [Logging](docs/logging.md)
