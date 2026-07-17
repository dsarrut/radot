# Resampling Module Documentation

The radot resampling module provides tools to resample medical images (CT, PET, SPECT, labels) to a target voxel spacing or to align them to match the grid of another reference image. It handles continuous values, discrete labels, and isotope activity preservation.

---

## 1. CLI Commands

### Resampling Images: radot resample

Resamples the input image. You must specify either a target spacing or a reference image to align with.

```bash
radot resample <image_path> [options]
```

#### Mutually Exclusive Target Options:
* -s, --spacing <x> <y> <z>: Resample to the specified voxel spacing in millimeters (e.g., -s 4 4 4).
* -l, --like <reference_image_path>: Resample the input image to match the spatial grid (size, spacing, origin, and direction) of the reference image.

#### Other Options:
* -t, --type <ct|labels|pet|spect|dose>: Modality type of the input image. If omitted, it is automatically resolved from the image's companion JSON sidecar.
* -o, --output <output_path>: Path for the resampled output file. If omitted, a name is automatically generated (e.g., image_4mm.nii.gz or image_like_reference.nii.gz).
* --unit <hu|suv|bqml|bq|gy>: Units of the image voxels. If omitted, it is resolved from the companion JSON sidecar. Specifying bq automatically enables total activity preservation.
* --preserve-total: Scales voxel values by the volume change ratio (new voxel volume / old voxel volume). This is critical when working with absolute activity (e.g., Bq) to ensure the total activity in an organ remains constant.
* --gaussian / --no-gaussian: Applies a Gaussian anti-aliasing filter along axes being downsampled. Enabled by default for continuous images to prevent aliasing artifacts. Disabled automatically for discrete label maps.

#### Behavior:
1. Interpolation Choice: Automatically uses linear interpolation for continuous modalities (CT, PET, SPECT) and nearest neighbor interpolation for discrete label maps.
2. Background Value: Sets background voxels outside the original image boundary to -1000 for CT scans (standard air value in HU) and 0 for all other image types.
3. Companion JSON Copying: If the input image has a companion JSON sidecar (e.g. image.json), the sidecar is automatically copied to match the output filename (e.g. output.json), preserving the metadata.

---

## 2. Python API

You can programmatically perform resampling in your scripts.

### resample_image_file_like

Resamples an image file to match the grid layout of a reference image file, and writes the output.

```python
from pathlib import Path
from radot import resample_image_file_like

# Resample a SPECT image to align perfectly with CT labels
resample_image_file_like(
    image_path=Path("spect.nii.gz"),
    like_image_path=Path("labels_all.nii.gz"),
    output_path=Path("spect_aligned.nii.gz"),
    type="spect",
    unit="bqml"
)
```

---

### resample_image

Performs in-memory resampling on SimpleITK images.

```python
import SimpleITK as sitk
from radot.resampling import resample_image

img = sitk.ReadImage("spect.nii.gz")

# Resample to 4mm isotropic spacing, preserving total activity
resampled_img = resample_image(
    image=img,
    spacing=(4.0, 4.0, 4.0),
    image_type="spect",
    preserve_total=True,
    gaussian=True
)
```

---

### resample_image_file

Resamples an image file to a target spacing and writes the output.

```python
from pathlib import Path
from radot.resampling import resample_image_file

resample_image_file(
    image_path=Path("ct.nii.gz"),
    spacing=(2.0, 2.0, 2.0),
    image_type="ct",
    output_path=Path("ct_2mm.nii.gz")
)
```
