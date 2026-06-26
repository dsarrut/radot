# Segmentation Module Documentation

The radot segmentation module provides tools to segment CT images using TotalSegmentator, automatically remapping the output label IDs and naming conventions to standard AAPM TG-263 guidelines.

---

## 1. CLI Commands

### Segmenting CT: radot segment

Runs TotalSegmentator on a CT scan, groups the resulting labels into a single multi-label image file, and applies the AAPM TG-263 mapping.

```bash
radot segment <ct_path> [options]
```

#### Options:
* -o, --output <output_path>: The output label image path. Defaults to labels_all.nii.gz, labels_head.nii.gz, or labels_<mode>.nii.gz depending on the chosen mode.
* --mode <mode_name>: The TotalSegmentator task/mode to run (e.g. total, head_glands_cavities, body, etc.). Defaults to total.
* --fast: Runs TotalSegmentator in fast mode (lower resolution, faster execution).

#### Behavior:
1. Multi-Label Execution: Runs TotalSegmentator in multi-label mode (ml=True) under the hood. This prevents creating dozens of individual files and instead outputs a single compressed label map.
2. AAPM TG-263 Remapping: Converts the TotalSegmentator class names and IDs to the standardized AAPM TG-263 nomenclature (e.g., remapping kidney_left to ID 5 and name Kidney_L).
3. Image Grid Alignment: Verifies that the spatial dimensions of the generated label map match the input CT scan. If there is a discrepancy, the label map is automatically resampled to the CT grid using Nearest Neighbor interpolation.
4. JSON Sidecar Generation: Writes a companion JSON file matching the output path (e.g., labels_all.json next to labels_all.nii.gz). This file maps the numeric label IDs to their TG-263 string names.

---

## 2. Python API

You can programmatically invoke the segmentation module and use its mapping helper functions.

### segment_ct

Runs the full segmentation, remapping, size checking, and JSON sidecar writing pipeline.

```python
from pathlib import Path
from radot.segmentation import segment_ct

ct_file = Path("ct.nii.gz")
output_labels = Path("output/labels.nii.gz")

# Run segmentation in fast mode
segment_ct(
    ct_path=ct_file,
    output_path=output_labels,
    mode="total",
    fast=True
)
```

#### Companion JSON Map Output Example:
The generated labels.json contains sorted key-value pairs matching the IDs to the standardized names:
```json
{
    "5": "Kidney_L",
    "6": "Kidney_R",
    "9": "Liver",
    "12": "Spleen",
    "45": "External"
}
```

---

### get_tg263_name & get_tg263_id

Helper functions to retrieve the mapped AAPM TG-263 names and integer IDs from raw TotalSegmentator class names.

```python
from radot.segmentation import get_tg263_name, get_tg263_id

# Retrieve the official TG-263 standardized name
name = get_tg263_name("kidney_left")
print(name)  # Prints: Kidney_L

# Retrieve the official TG-263 standardized ID
label_id = get_tg263_id("kidney_left")
print(label_id)  # Prints: 5

# Fallback behavior
unknown_name = get_tg263_name("unrecognized_class")
print(unknown_name)  # Prints: unrecognized_class
```

---

## 3. Remapping Details

The segmentation module maps major organs, glands, cardiovascular structures, bones, muscles, and cavities to standard TG-263 entries.

### Key Mappings Examples:
* liver -> ID 9, name Liver
* kidney_left -> ID 5, name Kidney_L
* kidney_right -> ID 6, name Kidney_R
* spleen -> ID 12, name Spleen
* brain -> ID 3, name Brain
* body -> ID 45, name External

### Dynamic Ribs and Vertebrae:
* Ribs: rib_left_1 through rib_left_12 and rib_right_1 through rib_right_12 are mapped to IDs 70-93 and named Bone_Rib_L_01 through Bone_Rib_R_12.
* Vertebrae: Cervical (C1-C7), Thoracic (T1-T12), Lumbar (L1-L5), and Sacral (S1) vertebrae are mapped to IDs 94-118 and named Bone_V_C01 through Bone_V_S01.

### Fallback Mechanism:
* Unmapped/Unrecognized classes are assigned a sequential ID starting at 1000 (or their original TotalSegmentator ID if the class name is not found in the mode map) and labeled as unknown_<id> to ensure compatibility and completeness.
