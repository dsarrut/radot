# Statistics Module Documentation

The radot stats module calculates volume, intensity, mass, and peak statistics for regions of interest (ROIs) defined by a label mask on medical images (CT, PET, SPECT).

---

## 1. CLI Commands

### Computing ROI Statistics: radot stats

Calculates statistics for each label present in a segmentation mask mapped onto an intensity image.

```bash
radot stats <labels_path> -i <image_path> [options]
```

#### Options:
* -i, --image <image_path>: Path to the intensity image (CT, PET, or SPECT) to analyze.
* -t, --type <ct|pet|spect>: Specifies the intensity image modality. If omitted, the module reads this from the image companion JSON sidecar.
* -o, --output <output_path>: Path to save the statistics as a JSON file. If omitted, a formatted table is printed directly to the console.
* --labels <id1,id2,...>: Comma-separated list of numeric label IDs to include (filters out other labels).
* --rois <name1,name2,...>: Comma-separated list of ROI names to include (filters out other ROIs).
* --resample-to <labels|image>: Specifies which grid layout to align to if the image grid does not match the labels grid. Defaults to labels.
* --sigma <float>: Specifies the standard deviation in millimeters of the Gaussian smoothing kernel used for peak value computation. Defaults to 5.0 mm for SPECT and 3.0 mm for PET.
* --unit <hu|suv|bqml|bq>: Specifies the image intensity units. Required for PET and SPECT if not present in the sidecar.

#### Behavior:
1. Grid Alignment: If the label mask and intensity image grids differ, the module automatically resamples one to match the other using the --resample-to strategy. Label maps use Nearest Neighbor interpolation; intensity images use Linear interpolation.
2. Modality-Specific Calculations:
   - CT Modality: Estimates physical tissue density (g/cm3) using the formula density = 1.0 + (mean_hu / 1000.0) and computes the total mass of each ROI in grams or kilograms.
   - PET/SPECT Modality: Computes the total activity in Bq or MBq within each ROI if the units are specified as bq or bqml.
3. Peak Value Estimation: Computes the maximum intensity value within each ROI after applying a Gaussian smoothing filter (using the specified --sigma).
4. Sidecar Mapping: If a companion JSON mapping file exists next to the labels file (e.g., labels.json next to labels.nii.gz), the module automatically maps numeric IDs to clinical names in the output table.

---

## 2. Python API

You can programmatically run statistics calculations directly in your Python scripts.

### compute_stats_file

Loads label and intensity files from disk, computes statistics, and prints a formatted table or saves a JSON file.

```python
from pathlib import Path
from radot.stats import compute_stats_file

# Compute CT stats and print the table to the console
compute_stats_file(
    labels_path=Path("labels.nii.gz"),
    image_path=Path("ct.nii.gz"),
    image_type="ct"
)
```

---

### compute_stats

Computes statistics on in-memory SimpleITK image objects and returns a Python dictionary.

```python
import SimpleITK as sitk
from radot.stats import compute_stats

labels_img = sitk.ReadImage("labels.nii.gz")
spect_img = sitk.ReadImage("spect.nii.gz")

# Compute statistics with manual labels_map
stats = compute_stats(
    labels_image=labels_img,
    image=spect_img,
    image_type="spect",
    unit="bqml",
    labels_map={5: "Kidney_L", 9: "Liver"},
    resample_to="labels"
)

# Access stats programmatically
liver_stats = stats["9"]
print(f"Liver Volume: {liver_stats['volume_ml']} mL")
print(f"Liver Peak Activity: {liver_stats['peak']} Bq/mL")
```

#### Output Dictionary Format:
```json
{
    "9": {
        "roi": "Liver",
        "volume_ml": 1820.5,
        "mean": 120.4,
        "std": 15.2,
        "min": 85.0,
        "max": 160.0,
        "sum": 219188.0,
        "peak": 154.1,
        "total_bq": 1753504.0
    },
    "all": {
        "roi": "Full Image",
        "volume_ml": 22000.0,
        "mean": 15.4,
        "std": 8.1,
        "min": -5.0,
        "max": 200.0,
        "sum": 338800.0,
        "peak": 195.0,
        "total_bq": 2710400.0
    }
}
```
