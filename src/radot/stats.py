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
    labels_image: sitk.Image,
    image: sitk.Image,
    image_type: Optional[str] = None,
    labels: Optional[List[int]] = None,
    rois: Optional[List[str]] = None,
    resample_to: str = "labels",
    labels_map: Optional[dict] = None,
    sigma: Optional[float] = None,
    unit: Optional[str] = None,
    image_path: Optional[Path] = None
) -> dict:
    """
    Compute volume and intensity statistics for each ROI/label from the input image.
    If the image does not share the same grid as the labels, it is resampled based on `resample_to`.
    Returns a dictionary of statistics per label.
    """
    from radot.meta import get_metadata_case_insensitive

    # 1. Resolve image_type (modality)
    if image_type is None:
        image_type = get_metadata_case_insensitive(image, "modality", image_path)
        if image_type is None:
            raise ValueError("Image modality type is not specified and was not found in the image metadata.")
    image_type = image_type.lower()

    # 2. Resolve unit
    if unit is None:
        unit = get_metadata_case_insensitive(image, "unit", image_path)
        if unit is not None:
            unit = unit.lower()
        elif image_type in ("pet", "spect"):
            raise ValueError("The unit option is required for PET and SPECT modalities (not found in image metadata).")

    # Check grid alignment and resample if necessary
    if not has_same_grid(labels_image, image):
        if resample_to == "labels":
            print("Image grid does not match labels grid. Resampling image to labels grid...")
            preserve_total = (unit is not None and unit.lower() == "bq")
            image = resample_image(
                image=image,
                image_type=image_type,
                preserve_total=preserve_total,
                gaussian=True,
                reference_image=labels_image
            )
        elif resample_to == "image":
            print("Labels grid does not match image grid. Resampling labels to image grid...")
            labels_image = resample_image(
                image=labels_image,
                image_type="labels",
                preserve_total=False,
                gaussian=False,
                reference_image=image
            )
        else:
            raise ValueError(f"Unknown resample_to target: {resample_to}")

    if labels_map is None:
        labels_map = {}

    # For peak value estimation: apply Gaussian smoothing with a physical sigma.
    # PET/CT: sigma = 3.0 mm (FWHM ~7mm)
    # SPECT: sigma = 5.0 mm (FWHM ~12mm)
    if sigma is None:
        sigma = 5.0 if image_type == "spect" else 3.0
    float_image = sitk.Cast(image, sitk.sitkFloat32)
    smoothed_image = sitk.SmoothingRecursiveGaussian(float_image, sigma)

    # Compute intensity statistics per label
    stat_filter = sitk.LabelIntensityStatisticsImageFilter()
    stat_filter.Execute(labels_image, image)
    
    # Compute statistics on smoothed image for peak value calculation
    stat_filter_peak = sitk.LabelIntensityStatisticsImageFilter()
    stat_filter_peak.Execute(labels_image, smoothed_image)
    
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
        peak_val = stat_filter_peak.GetMaximum(label_id)
        
        label_stats = {
            "roi": roi_name,
            "volume_ml": float(volume_ml),
            "mean": float(mean_val),
            "std": float(std_val),
            "min": float(min_val),
            "max": float(max_val),
            "sum": float(sum_val),
            "peak": float(peak_val)
        }
        
        if image_type == "ct":
            # Estimate density: density_g_cm3 = 1.0 + (mean_hu / 1000.0)
            # Ensure density is non-negative
            density = max(0.0, 1.0 + (mean_val / 1000.0))
            mass_g = volume_ml * density
            label_stats["mass_g"] = float(mass_g)
            
        stats_dict[str(label_id)] = label_stats

    return stats_dict


def compute_stats_file(
    labels_path: Path,
    image_path: Path,
    image_type: Optional[str] = None,
    output_path: Optional[Path] = None,
    labels: Optional[List[int]] = None,
    rois: Optional[List[str]] = None,
    resample_to: str = "labels",
    sigma: Optional[float] = None,
    unit: Optional[str] = None
) -> None:
    """
    Load images from paths, compute statistics, and save/print the output.
    """
    if not labels_path.exists():
        raise FileNotFoundError(f"Labels file not found: {labels_path}")
    if not image_path.exists():
        raise FileNotFoundError(f"Image file not found: {image_path}")

    print(f"Reading labels from: {labels_path}")
    labels_image = sitk.ReadImage(str(labels_path))
    
    print(f"Reading image from: {image_path}")
    image = sitk.ReadImage(str(image_path))

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

    stats_dict = compute_stats(
        labels_image=labels_image,
        image=image,
        image_type=image_type,
        labels=labels,
        rois=rois,
        resample_to=resample_to,
        labels_map=labels_map,
        sigma=sigma,
        unit=unit,
        image_path=image_path
    )

    # Output results
    if output_path:
        print(f"Saving statistics to: {output_path}")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(stats_dict, f, indent=4)
    else:
        # Print to console pretty printed
        pretty_print_stats(stats_dict, image_type)


def pretty_print_stats(stats_dict: dict, image_type: str) -> None:
    """Pretty print stats_dict as a table on the screen."""
    has_mass = (image_type == "ct")
    headers = ["Label", "ROI", "Volume (mL)", "Mean", "StdDev", "Min", "Max", "Peak"]
    if has_mass:
        headers.append("Mass (g)")
        
    rows = []
    for label_id, data in stats_dict.items():
        row = [
            str(label_id),
            data["roi"],
            f"{data['volume_ml']:.3f}",
            f"{data['mean']:.3f}",
            f"{data['std']:.3f}",
            f"{data['min']:.3f}",
            f"{data['max']:.3f}",
            f"{data['peak']:.3f}"
        ]
        if has_mass:
            mass_val = data.get("mass_g", 0.0)
            row.append(f"{mass_val:.3f}")
        rows.append(row)
        
    if not rows:
        print("No statistics computed.")
        return

    # Find max width for each column
    col_widths = [len(h) for h in headers]
    for row in rows:
        for i, val in enumerate(row):
            col_widths[i] = max(col_widths[i], len(val))
            
    # Format line template
    def format_row(row_data):
        parts = []
        for i, val in enumerate(row_data):
            width = col_widths[i]
            if i == 1:  # ROI name is left-aligned
                parts.append(val.ljust(width))
            else:  # Everything else is right-aligned
                parts.append(val.rjust(width))
        return " | ".join(parts)

    header_line = format_row(headers)
    separator = "-+-".join("-" * w for w in col_widths)
    
    print("\n" + header_line)
    print(separator)
    for row in rows:
        print(format_row(row))
    print()
