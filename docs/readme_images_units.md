
## Unit Conversion Command

### Converting Units: `radot convert-units`

Converts voxel values of radioactive images (SPECT/PET) between different units: `bq`, `bqml`, `suv`, and `sul`.

```bash
radot convert-units <image_path> -t <to_unit> [options]
```

#### Options:
* `-f, --from-unit <bq|bqml|suv|sul>`: Source unit of the image. If omitted, resolved from the image companion JSON sidecar.
* `-t, --to-unit <bq|bqml|suv|sul>`: Target unit to convert to (Required).
* `-o, --output <output_path>`: Optional path for the converted output file. If omitted, defaults to a modified filename like `<base>_<to-unit>.<ext>`.
* `--weight <float>`: Patient weight in kg (Required for SUV/SUL conversions, if not present in the sidecar).
* `--activity <float>`: Injected activity in MBq (Required for SUV/SUL conversions, if not present in the sidecar).
* `--height <float>`: Patient height in cm (Required for analytical SUL methods, if not present in the sidecar).
* `--sex <male|female>`: Patient sex (Required for analytical SUL methods, if not present in the sidecar).
* `--lbm-method <james|janmahasatian|segmentation>`: Method used to compute Lean Body Mass (LBM). Defaults to `james`.
* `-l, --labels <labels_path>`: Path to a label map image (Required if `--lbm-method segmentation` is used to subtract fat volume).

#### Fat Labels Resolution under `segmentation` method:
When using `--lbm-method segmentation`, the labels are resolved from the companion sidecar JSON of the `--labels` image. A label is identified as fat/adipose tissue (to be subtracted from total body weight) if its lowercase name satisfies either of:
1. Contains the substring `"fat"`, `"adipose"`, or `"lipoma"`.
2. Contains any of the exact clinical abbreviation tokens: `"sat"` (Subcutaneous Adipose Tissue), `"vat"` (Visceral Adipose Tissue), or `"tat"` (Torso Adipose Tissue).

This includes standard TotalSegmentator and TG-263 classes such as:
*   `subcutaneous_fat` / `Fat_Subq`
*   `visceral_fat` / `Fat_Visc`
*   `torso_fat`
*   `VAT_tissue` / `SAT_region`


#### Examples:

1. Convert a SPECT image from Bq/mL to Bq:
```bash
radot convert-units spect.nii.gz -t bq
```

2. Convert a PET image from Bq to SUV, specifying patient weight and injected activity explicitly:
```bash
radot convert-units pet.nii.gz -t suv --weight 75.0 --activity 250.0
```

3. Convert a PET image from SUV to SUL using the James formula:
```bash
radot convert-units pet.nii.gz -t sul --weight 75.0 --height 175.0 --sex male --lbm-method james
```

4. Convert a PET image from SUV to SUL using patient-specific segmentation (subtracting fat mass):
```bash
radot convert-units pet.nii.gz -t sul --weight 75.0 --lbm-method segmentation --labels labels_all.nii.gz
```


---

## Python API

You can also perform unit conversions programmatically within Python scripts.

### convert_image_units

Converts voxel values of a SimpleITK image in-memory.

```python
import SimpleITK as sitk
from radot import convert_image_units

# Load SPECT image
img = sitk.ReadImage("spect.nii.gz")

# Convert from Bq/mL (bqml) to SUV
suv_img = convert_image_units(
    image=img,
    from_unit="bqml",
    to_unit="suv",
    weight=75.0,
    activity=250.0
)

# Convert from SUV to SUL using the Janmahasatian formula
sul_img = convert_image_units(
    image=suv_img,
    from_unit="suv",
    to_unit="sul",
    weight=75.0,
    height=175.0,
    sex="male",
    lbm_method="janmahasatian"
)
```

---

### convert_image_units_file

Converts units directly between image files and saves the output (automatically updating/creating companion JSON sidecar metadata).

```python
from pathlib import Path
from radot.convert import convert_image_units_file

# Convert a PET file from Bq to SUL using James formula
convert_image_units_file(
    image_path=Path("pet.nii.gz"),
    to_unit="sul",
    from_unit="bq",
    output_path=Path("pet_sul.nii.gz"),
    weight=70.0,
    activity=220.0,
    height=180.0,
    sex="male",
    lbm_method="james"
)
```
