"""Input checks and compute_roi_dose, the entry point on SimpleITK images."""

import json
import logging
from typing import List, Optional, Dict, Any
import numpy as np
import SimpleITK as sitk
from pathlib import Path

from radot.meta import get_time_from_injection_h
from radot.dose.models import dose_ldm
from radot.dose.mtp import compute_roi_dose_mtp
from radot.dose.rate_map import DOSE_RATE_MAP_METHODS, dose_from_rate_map
from radot.dose.strategies import (
    DOSIMETRY_STRATEGIES,
    TIA_METHODS,
    TIA_UNITS,
    compute_dose_generic,
)

logger = logging.getLogger("radot")


def check_dosimetry_inputs(
    method: str,
    isotope: str,
    labels_path: Optional[Path] = None,
    rois: Optional[List[str]] = None,
    method_params: Optional[Dict[str, Any]] = None,
) -> Dict[int, str]:
    """
    Perform fast-fail boundary checks (method, isotope, ROI typos) before loading heavy images.
    Returns the loaded labels map.
    """
    method_key = method.lower()
    if (
        method_key not in DOSIMETRY_STRATEGIES
        and method_key != "ldm"
        and method_key != "jackson2013"
    ):
        raise ValueError(
            f"Dosimetry method '{method}' is not implemented yet. Available: {list(DOSIMETRY_STRATEGIES.keys()) + ['ldm', 'jackson2013']}"
        )

    if not isotope:
        raise ValueError("Isotope is required.")

    if method_key == "hanscheid2017" and isotope.lower() != "lu177":
        raise ValueError(
            "The Hanscheid 2017 method (0.125 scaling factor) is specific to Lu177 isotope."
        )

    # Load label mapping
    labels_map = {}
    if labels_path is not None:
        labels_path = Path(labels_path)
        from radot.meta import get_sidecar_path

        json_path = get_sidecar_path(labels_path)
        if json_path.exists():
            try:
                with open(json_path, "r") as f:
                    raw_map = json.load(f)
                    labels_map = {int(k): v for k, v in raw_map.items()}
            except Exception:
                pass

    if labels_map:
        norm_to_std = {v.lower(): v for v in labels_map.values()}
        known_rois = sorted(list(labels_map.values()))

        # 1. Validate & Normalize rois in-place
        if rois:
            normalized_rois = []
            for r in rois:
                r_lower = r.lower()
                if r_lower in norm_to_std:
                    normalized_rois.append(norm_to_std[r_lower])
                else:
                    import difflib

                    close_matches = difflib.get_close_matches(
                        r, known_rois, n=5, cutoff=0.4
                    )
                    if close_matches:
                        suggestions = " or ".join(f"'{m}'" for m in close_matches)
                        raise ValueError(
                            f"Requested ROI '{r}' not found. Did you mean {suggestions}?"
                        )
                    else:
                        raise ValueError(
                            f"Requested ROI '{r}' not found. Available ROIs: {known_rois}"
                        )
            rois.clear()
            rois.extend(normalized_rois)

        # 2. Validate & Normalize method_params in-place
        if method_params:
            normalized_params = {}
            for r, val in method_params.items():
                if r.lower() == "phantom":
                    normalized_params[r] = val
                    continue
                r_lower = r.lower()
                if r_lower in norm_to_std:
                    normalized_params[norm_to_std[r_lower]] = val
                else:
                    import difflib

                    close_matches = difflib.get_close_matches(
                        r, known_rois, n=5, cutoff=0.4
                    )
                    if close_matches:
                        suggestions = " or ".join(f"'{m}'" for m in close_matches)
                        raise ValueError(
                            f"Teff ROI '{r}' not found. Did you mean {suggestions}?"
                        )
                    else:
                        raise ValueError(
                            f"Teff ROI '{r}' not found. Available ROIs: {known_rois}"
                        )
            method_params.clear()
            method_params.update(normalized_params)

    return labels_map


def compute_roi_dose(
    image: Any,
    method: str,
    isotope: str,
    labels: Optional[Any] = None,
    rois: Optional[List[str]] = None,
    ct_image: Optional[Any] = None,
    total_activity_mbq: Optional[float] = None,
    method_params: Optional[Dict[str, Any]] = None,
    image_path: Optional[Path] = None,
    labels_path: Optional[Path] = None,
    labels_map: Optional[Dict[int, str]] = None,
    verbose: bool = False,
    ct_path: Optional[Any] = None,
    teff_params: Optional[Dict[str, float]] = None,
    hu_threshold: Optional[float] = -800.0,
    dose_rate: bool = False,
    cascade: bool = False,
    dose_rate_map: Optional[Any] = None,
    dose_rate_map_path: Optional[Path] = None,
    dose_rate_map_activity_mbq: Optional[float] = None,
) -> Any:
    """
    Computes average/voxel absorbed doses for specified regions of interest.
    With dose_rate_map (an external Monte Carlo dose-rate image), the map replaces the S-value
    or the local deposition and the method selects the time integration (see dose_from_rate_map).
    """
    # Check if multiple images (MTP) were passed
    if isinstance(image, (list, tuple)):
        if dose_rate_map is not None or dose_rate_map_path is not None:
            raise ValueError("A dose-rate map cannot be used with multiple time points.")
        return compute_roi_dose_mtp(
            images=list(image),
            labels=labels,
            method=method,
            isotope=isotope,
            times_h=method_params.get("times_h") if method_params else None,
            ct_image=ct_image,
            rois=rois,
            labels_map=labels_map,
            labels_path=labels_path,
            ct_path=ct_path,
            method_params=method_params,
            verbose=verbose,
        )

    from radot.stats import has_same_grid
    from radot.convert import convert_image_units
    from radot.meta import get_metadata_case_insensitive

    # Convert paths/strings to SimpleITK images on-the-fly
    if isinstance(image, (str, Path)):
        if image_path is None:
            image_path = Path(image)
        image = sitk.ReadImage(str(image))

    if labels is not None and isinstance(labels, (str, Path)):
        if labels_path is None:
            labels_path = Path(labels)
        labels = sitk.ReadImage(str(labels))

    if ct_image is None and ct_path is not None:
        ct_image = ct_path

    if isinstance(ct_image, (str, Path)):
        ct_image = sitk.ReadImage(str(ct_image))

    # Convert/merge teff_params into method_params structure
    if method_params is None:
        method_params = {}
    if teff_params:
        for roi, teff in teff_params.items():
            if roi not in method_params:
                method_params[roi] = {}
            method_params[roi]["teff"] = teff

    if dose_rate_map is None and dose_rate_map_path is not None:
        dose_rate_map = dose_rate_map_path
    if isinstance(dose_rate_map, (str, Path)):
        if dose_rate_map_path is None:
            dose_rate_map_path = Path(dose_rate_map)
        dose_rate_map = sitk.ReadImage(str(dose_rate_map))

    # 1. Fast-Fail Boundary Checks
    method_key = method.lower()
    if dose_rate_map is not None and method_key not in DOSE_RATE_MAP_METHODS:
        raise ValueError(
            f"Method '{method}' cannot be used with a dose-rate map. "
            f"Methods with a dose-rate map: {', '.join(DOSE_RATE_MAP_METHODS)}."
        )
    if labels_map is None:
        if method_key == "ldm":
            if labels_path is not None:
                labels_map = check_dosimetry_inputs(
                    method=method,
                    isotope=isotope,
                    labels_path=labels_path,
                    rois=rois,
                    method_params=method_params,
                )
            else:
                labels_map = {}
        else:
            labels_map = check_dosimetry_inputs(
                method=method,
                isotope=isotope,
                labels_path=labels_path,
                rois=rois,
                method_params=method_params,
            )

    # 2. Strict Grid Validation (skip labels validation for ldm)
    if method_key != "ldm":
        if not has_same_grid(image, labels):
            raise ValueError(
                "Activity image grid does not match labels grid, resample image first (option --resample)."
            )
        if ct_image is not None:
            if not has_same_grid(ct_image, labels):
                raise ValueError("CT image grid does not match labels grid.")
    else:
        if ct_image is not None:
            if not has_same_grid(image, ct_image):
                raise ValueError(
                    "CT image grid does not match activity image grid. Please resample CT first."
                )

    # 3. Unit and timing extraction
    unit = get_metadata_case_insensitive(image, "unit", image_path)
    if unit is None:
        raise ValueError("Image unit is not specified and was not found in metadata.")

    is_tia = unit.lower() in TIA_UNITS
    if is_tia and method_key not in TIA_METHODS:
        raise ValueError(
            f"Method '{method}' cannot be used with a time-integrated activity image "
            f"(unit '{unit}'). Methods for TIA images: {', '.join(TIA_METHODS)}."
        )
    if not is_tia and method_key == "svalue":
        raise ValueError(
            f"Method 'svalue' requires a time-integrated activity image "
            f"(unit {' or '.join(TIA_UNITS)}), got unit '{unit}'."
        )

    if is_tia:
        if total_activity_mbq is not None:
            raise ValueError(
                "--total-activity cannot be used with a time-integrated activity image."
            )
        if dose_rate:
            raise ValueError(
                "Dose rate cannot be computed from a time-integrated activity image."
            )
        # No time point: the image is already integrated over time
        time_from_injection_h = None
    else:
        time_from_injection_h = get_time_from_injection_h(image, image_path)

    # 4. Calibrated Activity Scaling
    if unit.lower() == "bqsml":
        # Bq.s/mL -> Bq.s per voxel
        sp = image.GetSpacing()
        image = image * ((sp[0] * sp[1] * sp[2]) / 1000.0)
    elif is_tia:
        pass
    elif unit.lower() != "bq":
        logger.info(
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
            logger.info(
                f"Normalizing total image activity to {total_activity_mbq} MBq (scale factor: {scale_factor:.6f})"
            )
            image = image * scale_factor

    if dose_rate_map is not None:
        return dose_from_rate_map(
            image=image,
            dose_rate_map=dose_rate_map,
            method=method,
            isotope=isotope,
            is_tia=is_tia,
            time_from_injection_h=time_from_injection_h,
            labels=labels,
            labels_map=labels_map,
            rois=rois,
            ct_image=ct_image,
            method_params=method_params,
            dose_rate_map_path=dose_rate_map_path,
            dose_rate_map_activity_mbq=dose_rate_map_activity_mbq,
            dose_rate=dose_rate,
            verbose=verbose,
        )

    # 5. Extract 3D data and volumes
    spacing = image.GetSpacing()
    voxel_volume_ml = (spacing[0] * spacing[1] * spacing[2]) / 1000.0

    image_arr = sitk.GetArrayFromImage(image)

    if method_key == "ldm":
        from radot.db import get_physical_half_life_h

        t_p_h = get_physical_half_life_h(isotope)
        ct_arr = sitk.GetArrayFromImage(ct_image) if ct_image is not None else None

        # Call the dedicated ldm dosimetry function
        dose_arr = dose_ldm(
            image_arr=image_arr,
            voxel_volume_ml=voxel_volume_ml,
            time_from_injection_h=time_from_injection_h,
            isotope=isotope,
            ct_arr=ct_arr,
            hu_threshold=hu_threshold,
            dose_rate=dose_rate,
            cascade=cascade,
            tia=is_tia,
        )

        # Create SimpleITK image matching activity image space
        dose_image = sitk.GetImageFromArray(dose_arr.astype(np.float32))
        dose_image.CopyInformation(image)
        dose_image.SetMetaData("unit", "gy/h" if dose_rate else "gy")
        dose_image.SetMetaData("modality", "dose")

        if labels is None:
            return dose_image

        labels_arr = sitk.GetArrayFromImage(labels)

        # Compute voxel mass array for region-absorbed dose aggregation
        if ct_image is not None:
            density_arr = np.clip(1.0 + (ct_arr / 1000.0), 1e-3, None)
        else:
            density_arr = np.ones_like(image_arr)
        voxel_mass_g = voxel_volume_ml * density_arr

        results = {}
        for label_val, roi_name in labels_map.items():
            if rois and roi_name not in rois:
                continue
            mask = labels_arr == label_val
            if not np.any(mask):
                continue

            roi_mass_g = np.sum(voxel_mass_g[mask])
            mean_dose_gy = (
                np.sum(dose_arr[mask] * voxel_mass_g[mask]) / roi_mass_g
                if roi_mass_g > 0
                else 0.0
            )
            roi_volume_ml = np.sum(mask) * voxel_volume_ml

            if dose_rate:
                results[roi_name] = {
                    "mean_dose_rate_gyh": float(mean_dose_gy),
                    "roi_mass_g": float(roi_mass_g),
                    "roi_volume_ml": float(roi_volume_ml),
                }
            elif is_tia:
                tia_mbq_s = float(np.sum(image_arr[mask]) / 1e6)
                results[roi_name] = {
                    "mean_dose_gy": float(mean_dose_gy),
                    "roi_mass_g": float(roi_mass_g),
                    "roi_volume_ml": float(roi_volume_ml),
                    "tia_mbq_s": tia_mbq_s,
                    "tia_mbq_h": tia_mbq_s / 3600.0,
                }
            else:
                results[roi_name] = {
                    "mean_dose_gy": float(mean_dose_gy),
                    "roi_mass_g": float(roi_mass_g),
                    "roi_volume_ml": float(roi_volume_ml),
                    "roi_teff_h": float(t_p_h),
                }
        return results, dose_image

    labels_arr = sitk.GetArrayFromImage(labels)
    ct_arr = sitk.GetArrayFromImage(ct_image) if ct_image is not None else None

    # Call the generic dosimetry pipeline
    return compute_dose_generic(
        image_arr=image_arr,
        labels_arr=labels_arr,
        isotope=isotope,
        voxel_volume_ml=voxel_volume_ml,
        time_from_injection_h=time_from_injection_h,
        labels_map=labels_map,
        ct_arr=ct_arr,
        rois=rois,
        method_params=method_params,
        verbose=verbose,
        method=method,
    )
