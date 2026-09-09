# RADOT (Radiation Dosimetry for Therapy)

`radot` is a Python library and command-line tool suite for radiation dosimetry, image processing, and metadata workflows in targeted radionuclide therapy.

## Installation

```bash
pip install -e .
```

To enable shell autocompletion for the CLI (zsh):
```bash
eval "$(_RADOT_COMPLETE=zsh_source radot)"
```

## Quick Start

See the [Quickstart Guide](docs/quickstart.md) to get started with CLI commands and Python workflows.

## Core Features

1. **DICOM Conversion**: Converts CT, SPECT, and PET series to NIfTI/NRRD with structured JSON metadata sidecars.
2. **Segmentation**: Automates CT segmentation via TotalSegmentator and remaps labels to standard AAPM TG-263 nomenclature.
3. **Unit Conversion**: Converts image voxel values between `bq`, `bqml`, `suv`, and `sul` (using James, Janmahasatian, or CT fat segmentation).
4. **Image Resampling**: Isotropic and grid-matching resampling while strictly preserving total radioactivity.
5. **Radionuclide Physics & Decay Chains**: Queries ICRP Publication 107 emissions spectra, IAEA LiveChart decay cascades (for alpha emitters like $^{225}\text{Ac}$), and OpenDose S-values.
6. **Dosimetry Models**: Implements Local Deposition Model (LDM, voxel-level with daughter cascade support), Hänscheid (2017/2018), and Madsen (2018) models.

## Documentation

- [Quickstart](docs/quickstart.md)
- [Dosimetry Models](docs/readme_dose.md)
- [Units & SUL Normalization](docs/readme_images_units.md) | [SUL CT Integration](docs/readme_sul.md)
- [Radionuclides & Decay Chains](docs/readme_radionuclide.md)
- [DICOM Conversion](docs/readme_dicom.md) & [Metadata Sidecars](docs/readme_meta.md)
- [CT Segmentation & TG-263](docs/readme_segmentation.md) | [ROI Names Policy](docs/readme_roi_names.md)
- [Image Resampling](docs/readme_images_resample.md)
- [OpenDose Database](docs/readme_opendose.md)
- [Statistics & Reporting](docs/readme_stats.md) | [SUVpeak](docs/readme_peak.md)
- [Python API Reference](docs/readme_python.md)


