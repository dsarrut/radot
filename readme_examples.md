# Examples

## Segmentation

To segment a CT image using TotalSegmentator with dynamic AAPM TG-263 remapping:

```bash
# By default, runs total mode. Outputs labels_all.nii.gz and labels_all.json
radot segment ct.nii.gz

# Segment using head_glands_cavities mode. Outputs labels_head.nii.gz and labels_head.json
radot segment ct.nii.gz --mode head_glands_cavities

# Specify a custom output path. Outputs my_labels.nii.gz and my_labels.json
radot segment ct.nii.gz -o my_labels.nii.gz
```

For each run, two output files are created:
1. A multi-label `.nii.gz` image containing the segmentation labels remapped to official TG-263 IDs.
2. A matching `.json` file outlining the index-to-name mapping of the generated labels.

Can be visualized with vvv : 

```bash
vvv ct.nii.gz + labels_total.nii.gz 


## Stats

```bash
radot stats labels_all.nii.gz -i ct.nii.gz -t ct 
```

## Convert Units

To convert voxel units between absolute activity (Bq), concentration (Bq/mL), and standardized uptake value (SUV):

```bash
# Convert activity concentration (bqml) to SUV
radot convert-units image.nii.gz -f bqml -t suv --weight 75 --activity 250
```

* **`-f` / `--from-unit`**: Source unit (`bq`, `bqml`, `suv`).
* **`-t` / `--to-unit`**: Target unit (`bq`, `bqml`, `suv`).
* **`-o` / `--output`**: Optional output path. Defaults to `<basename>_<to-unit>.<ext>` (e.g. `image_suv.nii.gz`).
* **`--weight`**: Patient weight in kg (required for SUV conversions).
* **`--activity`**: Injected activity in MBq (required for SUV conversions).

## Image Metadata

To read or write custom dosimetry metadata (modality, unit, patient weight, injected activity) inside NIfTI, MHD, or NRRD image headers:

```bash
# Write metadata to an image
radot meta image.mhd --modality PET --unit bq --weight 70 --activity 180

# Read/print stored metadata of an image
radot meta image.nii.gz
```

* **`--modality`**: Modality type (`ct`, `pet`, `spect`).
* **`--unit`**: Voxel values unit (`hu`, `suv`, `bqml`, `bq`).
* **`--weight`**: Patient weight in kg.
* **`--activity`**: Injected activity in MBq.

