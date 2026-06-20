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

### Outputs
For each run, two output files are created:
1. A multi-label `.nii.gz` image containing the segmentation labels remapped to official TG-263 IDs.
2. A matching `.json` file outlining the index-to-name mapping of the generated labels.

Can be visualized with vvv : 

```bash
vvv ct.nii.gz + labels_total.nii.gz 

