import json
from pathlib import Path
from typing import List, Optional
import numpy as np
import SimpleITK as sitk
from radot.resampling import resample_image

def has_same_grid(img1: sitk.Image, img2: sitk.Image) -> bool:
    """Check if two SimpleITK images share the exact same grid (size, spacing, origin, direction)."""
    if img1.GetSize() != img2.GetSize():
        return False
    if not np.allclose(img1.GetSpacing(), img2.GetSpacing(), atol=1e-5):
        return False
    if not np.allclose(img1.GetOrigin(), img2.GetOrigin(), atol=1e-5):
        return False
    if not np.allclose(img1.GetDirection(), img2.GetDirection(), atol=1e-5):
        return False
    return True

def compute_stats(
    labels_path: Path,
    image_path: Path,
    image_type: str,
    output_path: Optional[Path] = None,
    labels: Optional[List[int]] = None,
    rois: Optional[List[str]] = None
) -> None:
    """
    Compute volume and intensity statistics for each ROI/label from the input image.
    If the image does not share the same grid as the labels, it is resampled to the labels grid.
    """
    if not labels_path.exists():
        raise FileNotFoundError(f"Labels file not found: {labels_path}")
    if not image_path.exists():
        raise FileNotFoundError(f"Image file not found: {image_path}")

    print(f"Reading labels from: {labels_path}")
    labels_image = sitk.ReadImage(str(labels_path))
    
    print(f"Reading image from: {image_path}")
    image = sitk.ReadImage(str(image_path))

    # Check grid alignment and resample if necessary
    if not has_same_grid(labels_image, image):
        print("Image grid does not match labels grid. Resampling image to labels grid...")
        image = resample_image(
            image=image,
            image_type=image_type,
            preserve_total=False,
            gaussian=True,
            reference_image=labels_image
        )

    # Load mapping JSON if it exists in the same folder
    json_path = labels_path.with_suffix("").with_suffix(".json")
    labels_map = {}
    if json_path.exists():
        try:
            with open(json_path, "r") as f:
                raw_map = json.load(f)
                labels_map = {int(k): v for k, v in raw_map.items()}
            print(f"Loaded label name mappings from: {json_path}")
        except Exception as e:
            print(f"Warning: Could not read mapping JSON {json_path}: {e}")

    # Compute intensity statistics per label
    stat_filter = sitk.LabelIntensityStatisticsImageFilter()
    stat_filter.Execute(labels_image, image)
    
    stats_dict = {}
    
    # Sort label IDs for ordered console/JSON printing
    sorted_labels = sorted(stat_filter.GetLabels())
    
    for label_id in sorted_labels:
        roi_name = labels_map.get(label_id, f"unknown_{label_id}")
        
        # Apply filters if specified
        if labels is not None and label_id not in labels:
            continue
        if rois is not None and roi_name not in rois:
            continue
            
        volume_ml = stat_filter.GetPhysicalSize(label_id) / 1000.0  # mm^3 to mL
        mean_val = stat_filter.GetMean(label_id)
        std_val = stat_filter.GetStandardDeviation(label_id)
        min_val = stat_filter.GetMinimum(label_id)
        max_val = stat_filter.GetMaximum(label_id)
        sum_val = stat_filter.GetSum(label_id)
        
        label_stats = {
            "roi": roi_name,
            "volume_ml": float(volume_ml),
            "mean": float(mean_val),
            "std": float(std_val),
            "min": float(min_val),
            "max": float(max_val),
            "sum": float(sum_val)
        }
        
        if image_type == "ct":
            # Estimate density: density_g_cm3 = 1.0 + (mean_hu / 1000.0)
            # Ensure density is non-negative
            density = max(0.0, 1.0 + (mean_val / 1000.0))
            mass_g = volume_ml * density
            label_stats["mass_g"] = float(mass_g)
            
        stats_dict[str(label_id)] = label_stats


    # Output results
    if output_path:
        print(f"Saving statistics to: {output_path}")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(stats_dict, f, indent=4)
    else:
        # Print to console
        print(json.dumps(stats_dict, indent=4))
