import json
from typing import List, Optional, Dict, Any
import numpy as np
import SimpleITK as sitk
from pathlib import Path
from datetime import datetime

ISOTOPE_HALF_LIVES_H = {
    "lu177": 6.647 * 24.0,  # 159.528 hours
    "y90": 64.08,  # 64.08 hours
    "i131": 8.02 * 24.0,  # 192.48 hours
}


def get_time_from_injection_h(
    image: sitk.Image, image_path: Optional[Path] = None
) -> float:
    from radot.meta import get_metadata_case_insensitive

    # 1. Try reading directly
    t_str = get_metadata_case_insensitive(image, "time_from_injection", image_path)
    if t_str is not None:
        try:
            return float(t_str)
        except ValueError:
            pass

    # 2. Try calculating from acquisition and injection dates
    inj_val = get_metadata_case_insensitive(image, "injection_date", image_path)
    acq_val = get_metadata_case_insensitive(image, "acquisition_date", image_path)
    if inj_val and acq_val:
        try:
            inj_dt = datetime.strptime(inj_val, "%Y-%m-%d %H:%M:%S")
            acq_dt = datetime.strptime(acq_val, "%Y-%m-%d %H:%M:%S")
            time_diff = acq_dt - inj_dt
            return time_diff.total_seconds() / 3600.0
        except Exception:
            pass

    raise ValueError(
        "Time from injection is required but was not found in the image metadata."
    )


def dose_hanscheid2017(
    spect_bq_a: np.ndarray,
    roi_a: np.ndarray,
    time_from_injection_h: float,
    volume_voxel_m_l: float,
    effective_time_h: float,
) -> float:
    """
    Input img and ROI must be numpy arrays
    Input SPECT img must be in Bq
    """
    # compute mean activity concentration in the ROI
    v = spect_bq_a[roi_a == 1] / volume_voxel_m_l
    concentration = np.mean(v) / 1e6

    # compute dose
    dose = (
        0.125
        * concentration
        * np.power(2, time_from_injection_h / effective_time_h)
        * effective_time_h
    )

    return float(dose)


def compute_dose_hanscheid2017(
    image_arr: np.ndarray,
    labels_arr: np.ndarray,
    isotope: str,
    voxel_volume_ml: float,
    time_from_injection_h: float,
    labels_map: Dict[int, str],
    ct_arr: Optional[np.ndarray] = None,
    rois: Optional[List[str]] = None,
    method_params: Optional[Dict[str, Any]] = None,
) -> dict:
    """
    Implementation of the Hanscheid 2017 dosimetry method (Option B).
    """
    if isotope.lower() != "lu177":
        raise ValueError(
            "The Hanscheid 2017 method (0.125 scaling factor) is specific to Lu177 isotope."
        )

    unique_labels = sorted(np.unique(labels_arr))
    results = {}

    for label_id in unique_labels:
        if label_id == 0:
            continue

        roi_name = labels_map.get(label_id, f"unknown_{label_id}")

        # Filter by ROIs if specified
        if rois is not None and roi_name not in rois:
            continue

        mask = labels_arr == label_id
        if not np.any(mask):
            continue

        # Get effective half-life (Teff)
        teff = None
        if method_params and roi_name in method_params:
            teff = method_params[roi_name].get("teff")

        if teff is None:
            # Fall back to physical half-life
            teff = ISOTOPE_HALF_LIVES_H.get(isotope.lower())
            if teff is None:
                raise ValueError(
                    f"No effective half-life (teff) provided for ROI '{roi_name}' and unknown physical half-life for isotope '{isotope}'."
                )
            print(
                f"No Teff specified for '{roi_name}'. Using physical half-life of {isotope} ({teff:.2f} h) as fallback."
            )

        # compute dose using the specific hanscheid2017 method
        dose = dose_hanscheid2017(
            spect_bq_a=image_arr,
            roi_a=mask.astype(np.uint8),
            time_from_injection_h=time_from_injection_h,
            volume_voxel_m_l=voxel_volume_ml,
            effective_time_h=teff,
        )

        # compute volume and mass
        roi_volume_ml = np.sum(mask) * voxel_volume_ml

        if ct_arr is not None:
            # density = 1.0 + (HU / 1000)
            density = np.maximum(0.0, 1.0 + (ct_arr[mask] / 1000.0))
            roi_mass_g = np.sum(density) * voxel_volume_ml
        else:
            # Fall back to water density (1.0 g/mL)
            roi_mass_g = roi_volume_ml

        results[roi_name] = {
            "mean_dose_gy": float(dose),
            "roi_mass_g": float(roi_mass_g),
            "roi_volume_ml": float(roi_volume_ml),
            "roi_teff_h": float(teff),
        }

    return results


DOSIMETRY_METHODS = {"hanscheid2017": compute_dose_hanscheid2017}


def compute_roi_dose(
    image: sitk.Image,
    labels: sitk.Image,
    method: str,
    isotope: str,
    rois: Optional[List[str]] = None,
    ct_image: Optional[sitk.Image] = None,
    total_activity_mbq: Optional[float] = None,
    method_params: Optional[Dict[str, Any]] = None,
    image_path: Optional[Path] = None,
    labels_path: Optional[Path] = None,
) -> dict:
    """
    Computes average/voxel absorbed doses for specified regions of interest.
    """
    from radot.stats import has_same_grid
    from radot.convert import convert_image_units
    from radot.meta import get_metadata_case_insensitive

    # 1. Fast-Fail Boundary Checks
    method_key = method.lower()
    if method_key not in DOSIMETRY_METHODS:
        raise ValueError(
            f"Dosimetry method '{method}' is not implemented yet. Available: {list(DOSIMETRY_METHODS.keys())}"
        )

    if not isotope:
        raise ValueError("Isotope is required.")

    # Load label mapping
    labels_map = {}
    if labels_path is not None:
        json_path = labels_path.with_suffix("").with_suffix(".json")
        if json_path.exists():
            try:
                with open(json_path, "r") as f:
                    raw_map = json.load(f)
                    labels_map = {int(k): v for k, v in raw_map.items()}
            except Exception:
                pass

    if rois:
        if labels_map:
            known_rois = set(labels_map.values())
            invalid_rois = [r for r in rois if r not in known_rois]
            if invalid_rois:
                raise ValueError(
                    f"Requested ROI(s) {invalid_rois} not found in the label mapping sidecar."
                )

    # 2. Strict Grid Validation
    if not has_same_grid(image, labels):
        raise ValueError(
            "Activity image grid does not match labels grid, resample image first."
        )

    if ct_image is not None:
        if not has_same_grid(ct_image, labels):
            raise ValueError("CT image grid does not match labels grid.")

    # 3. Timing Extraction
    time_from_injection_h = get_time_from_injection_h(image, image_path)

    # 4. Calibrated Activity Scaling
    unit = get_metadata_case_insensitive(image, "unit", image_path)
    if unit is None:
        raise ValueError("Image unit is not specified and was not found in metadata.")

    if unit.lower() != "bq":
        print(
            f"Converting image units on-the-fly from {unit} to bq for dose computation..."
        )
        image = convert_image_units(
            image, to_unit="bq", from_unit=unit, image_path=image_path
        )

    # Apply total activity override normalization if requested
    if total_activity_mbq is not None:
        stat_filter = sitk.StatisticsImageFilter()
        stat_filter.Execute(image)
        current_total_bq = stat_filter.GetSum()
        if current_total_bq > 0:
            scale_factor = (total_activity_mbq * 1e6) / current_total_bq
            print(
                f"Normalizing total image activity to {total_activity_mbq} MBq (scale factor: {scale_factor:.6f})"
            )
            image = image * scale_factor

    # 5. Extract 3D data and volumes
    spacing = image.GetSpacing()
    voxel_volume_ml = (spacing[0] * spacing[1] * spacing[2]) / 1000.0

    image_arr = sitk.GetArrayFromImage(image)
    labels_arr = sitk.GetArrayFromImage(labels)
    ct_arr = sitk.GetArrayFromImage(ct_image) if ct_image is not None else None

    # Call the specific dosimetry method
    method_func = DOSIMETRY_METHODS[method_key]
    return method_func(
        image_arr=image_arr,
        labels_arr=labels_arr,
        isotope=isotope,
        voxel_volume_ml=voxel_volume_ml,
        time_from_injection_h=time_from_injection_h,
        labels_map=labels_map,
        ct_arr=ct_arr,
        rois=rois,
        method_params=method_params,
    )


def compute_roi_dose_file(
    image_path: Path,
    labels_path: Path,
    method: str,
    isotope: str,
    ct_path: Optional[Path] = None,
    total_activity_mbq: Optional[float] = None,
    rois: Optional[List[str]] = None,
    teff_params: Optional[Dict[str, float]] = None,
    output_path: Optional[Path] = None,
    resample_target: Optional[str] = None,
) -> None:
    """
    File interface to compute ROI dose and print/save the results.
    """
    if not image_path.exists():
        raise FileNotFoundError(f"Activity image not found: {image_path}")
    if not labels_path.exists():
        raise FileNotFoundError(f"Labels image not found: {labels_path}")

    print(f"Reading activity image: {image_path}")
    image = sitk.ReadImage(str(image_path))

    print(f"Reading labels image: {labels_path}")
    labels = sitk.ReadImage(str(labels_path))

    ct_image = None
    if ct_path:
        if not ct_path.exists():
            raise FileNotFoundError(f"CT image not found: {ct_path}")
        print(f"Reading companion CT image: {ct_path}")
        ct_image = sitk.ReadImage(str(ct_path))

    # Perform resampling if target is specified
    if resample_target:
        print(f"Resampling images to align grids: target is {resample_target}")
        from radot.resampling import resample_image
        from radot.meta import get_metadata_case_insensitive
        
        modality = get_metadata_case_insensitive(image, "modality", image_path) or "spect"
        unit = get_metadata_case_insensitive(image, "unit", image_path)
        preserve_total = (unit is not None and unit.lower() == "bq")

        if resample_target == "like_activity":
            labels = resample_image(labels, image_type="labels", reference_image=image)
            if ct_image is not None:
                ct_image = resample_image(ct_image, image_type="ct", preserve_total=False, gaussian=True, reference_image=image)
        elif resample_target == "like_labels":
            image = resample_image(image, image_type=modality.lower(), preserve_total=preserve_total, gaussian=True, reference_image=labels)
            if ct_image is not None:
                ct_image = resample_image(ct_image, image_type="ct", preserve_total=False, gaussian=True, reference_image=labels)
        elif resample_target == "like_ct":
            if ct_image is None:
                raise ValueError("Cannot resample like CT image because companion CT path was not provided.")
            image = resample_image(image, image_type=modality.lower(), preserve_total=preserve_total, gaussian=True, reference_image=ct_image)
            labels = resample_image(labels, image_type="labels", reference_image=ct_image)
        else:
            raise ValueError(f"Unknown resample target: {resample_target}")

    # Convert teff_params into method_params structure
    method_params = {}
    if teff_params:
        for roi, teff in teff_params.items():
            method_params[roi] = {"teff": teff}

    results = compute_roi_dose(
        image=image,
        labels=labels,
        method=method,
        isotope=isotope,
        rois=rois,
        ct_image=ct_image,
        total_activity_mbq=total_activity_mbq,
        method_params=method_params,
        image_path=image_path,
        labels_path=labels_path,
    )

    # Print results to stdout
    if output_path:
        print(f"Saving dose results to: {output_path}")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(results, f, indent=4)
    else:
        # Pretty print results as a table
        pretty_print_dose(results)


def pretty_print_dose(results: dict) -> None:
    """Pretty print dosimetry results as a table."""
    headers = ["ROI", "Mean Dose (Gy)", "Mass (g)", "Volume (mL)", "Teff (h)"]
    rows = []
    for roi, data in results.items():
        rows.append(
            [
                roi,
                f"{data['mean_dose_gy']:.3f}",
                f"{data['roi_mass_g']:.3f}",
                f"{data['roi_volume_ml']:.3f}",
                f"{data['roi_teff_h']:.1f}",
            ]
        )

    if not rows:
        print("No dose computed.")
        return

    # Compute column widths
    col_widths = [len(h) for h in headers]
    for row in rows:
        for i, val in enumerate(row):
            col_widths[i] = max(col_widths[i], len(val))

    def format_row(row_data):
        parts = []
        for i, val in enumerate(row_data):
            width = col_widths[i]
            if i == 0:  # left aligned
                parts.append(val.ljust(width))
            else:  # right aligned
                parts.append(val.rjust(width))
        return " | ".join(parts)

    header_line = format_row(headers)
    separator = "-+-".join("-" * w for w in col_widths)

    print("\n" + header_line)
    print(separator)
    for row in rows:
        print(format_row(row))
    print()
