"""Dose from an external dose-rate map (e.g. OpenGATE Monte Carlo)."""

import logging
from typing import List, Optional, Dict, Any, Tuple
import numpy as np
import SimpleITK as sitk
from pathlib import Path

from radot.db import get_physical_half_life_h

logger = logging.getLogger("radot")


# Methods accepting an external dose-rate map (e.g. OpenGATE): the map replaces the S-value
# or the local deposition, the method only selects the time integration
DOSE_RATE_MAP_METHODS = ("hanscheid2018", "madsen2018", "ldm", "svalue")
# Dose-rate map units, factor to Gy/s
DOSE_RATE_MAP_UNITS = {"gy/s": 1.0, "gy/h": 1.0 / 3600.0}


def dose_from_rate_map(
    image: sitk.Image,
    dose_rate_map: sitk.Image,
    method: str,
    isotope: str,
    is_tia: bool,
    time_from_injection_h: Optional[float],
    labels: Optional[sitk.Image] = None,
    labels_map: Optional[Dict[int, str]] = None,
    rois: Optional[List[str]] = None,
    ct_image: Optional[sitk.Image] = None,
    method_params: Optional[Dict[str, Any]] = None,
    dose_rate_map_path: Optional[Path] = None,
    dose_rate_map_activity_mbq: Optional[float] = None,
    dose_rate: bool = False,
    verbose: bool = False,
) -> Any:
    """
    Absorbed dose from an external dose-rate map (e.g. OpenGATE opengate_dose_rate output).
    - image: activity in Bq per voxel, or TIA in Bq.s per voxel if is_tia
    - dose_rate_map: dose rate (Gy/s or Gy/h, sidecar 'unit', default Gy/s) produced by
      dose_rate_map_activity_mbq distributed like the image. It is scaled to the image total:
      the result is the dose rate at the acquisition time (activity image), or the absorbed
      dose (TIA image). Without the map activity (argument or sidecar 'source_activity_mbq'),
      the map is taken as the absolute dose rate at the acquisition time.
    - method: time integration, hanscheid2018 (Teff / ln2, Teff = 2t by default),
      madsen2018 (exp(kt) / k), ldm (physical decay), or none for a TIA image (svalue, ldm)
    - ROI value: mass-weighted mean of the map in the ROI
    Returns the results dict, (results, dose image) for ldm, or the dose image for ldm
    without labels.
    """
    from radot.meta import get_metadata_case_insensitive
    from radot.stats import has_same_grid

    method_key = method.lower()
    if not has_same_grid(dose_rate_map, image):
        raise ValueError(
            "Dose-rate map grid does not match the activity image grid "
            "(the map must be computed with the activity image as source)."
        )

    # Map unit
    unit = get_metadata_case_insensitive(dose_rate_map, "unit", dose_rate_map_path)
    if unit is None:
        unit = "gy/s"
        logger.info("Dose-rate map unit not specified: assuming gy/s (OpenGATE output).")
    if unit.lower() not in DOSE_RATE_MAP_UNITS:
        raise ValueError(
            f"Dose-rate map unit must be one of {', '.join(DOSE_RATE_MAP_UNITS)}, got '{unit}'."
        )

    # Normalisation to the image total (Bq, or Bq.s for a TIA image)
    if dose_rate_map_activity_mbq is None:
        a = get_metadata_case_insensitive(dose_rate_map, "source_activity_mbq", dose_rate_map_path)
        dose_rate_map_activity_mbq = float(a) if a is not None else None
    image_arr = sitk.GetArrayFromImage(image)
    image_total = float(np.sum(image_arr))
    if dose_rate_map_activity_mbq is not None:
        if dose_rate_map_activity_mbq <= 0:
            raise ValueError("Dose-rate map activity must be positive.")
        scale = image_total / (dose_rate_map_activity_mbq * 1e6)
    elif is_tia:
        raise ValueError(
            "With a time-integrated activity image, the activity of the dose-rate map is "
            "required (--dose-rate-map-activity or sidecar key 'source_activity_mbq')."
        )
    else:
        scale = 1.0
        logger.info("No dose-rate map activity: the map is used as the absolute dose rate.")
    # Gy/s at the acquisition time (activity image), or Gy (TIA image)
    rate_arr = (
        sitk.GetArrayFromImage(dose_rate_map).astype(np.float64)
        * DOSE_RATE_MAP_UNITS[unit.lower()]
        * scale
    )

    ln2 = np.log(2)
    t_s = None if is_tia else time_from_injection_h * 3600.0

    def time_factor_s(roi_name: Optional[str]) -> Tuple[float, Optional[float]]:
        """Time-integration factor (s) and Teff (h) for one ROI."""
        if is_tia:
            return 1.0, None
        if method_key == "ldm":
            teff = get_physical_half_life_h(isotope)
            k = ln2 / (teff * 3600.0)
            return float(np.exp(k * t_s) / k), teff
        teff = None
        if method_params and roi_name in method_params:
            teff = method_params[roi_name].get("teff")
        if method_key == "hanscheid2018":
            if teff is None:
                teff = 2.0 * time_from_injection_h
            return teff * 3600.0 / ln2, teff
        if teff is None:
            raise ValueError(
                f"Effective half-life (teff) must be provided for ROI '{roi_name}' when using the madsen2018 method."
            )
        k = ln2 / (teff * 3600.0)
        return float(np.exp(k * t_s) / k), teff

    spacing = image.GetSpacing()
    voxel_volume_ml = (spacing[0] * spacing[1] * spacing[2]) / 1000.0
    if ct_image is not None:
        density_arr = np.clip(1.0 + (sitk.GetArrayFromImage(ct_image) / 1000.0), 1e-3, None)
    else:
        density_arr = np.ones_like(rate_arr)
    voxel_mass_g = voxel_volume_ml * density_arr

    dose_image = None
    if method_key == "ldm":
        if dose_rate:
            dose_arr = rate_arr * 3600.0
        else:
            dose_arr = rate_arr * time_factor_s(None)[0]
        dose_image = sitk.GetImageFromArray(dose_arr.astype(np.float32))
        dose_image.CopyInformation(image)
        dose_image.SetMetaData("unit", "gy/h" if dose_rate else "gy")
        dose_image.SetMetaData("modality", "dose")
        if labels is None:
            return dose_image

    labels_arr = sitk.GetArrayFromImage(labels)
    results = {}
    for label_val, roi_name in (labels_map or {}).items():
        if rois and roi_name not in rois:
            continue
        mask = labels_arr == label_val
        if not np.any(mask):
            continue
        roi_mass_g = float(np.sum(voxel_mass_g[mask]))
        roi_volume_ml = float(np.sum(mask) * voxel_volume_ml)
        mean_rate = float(np.sum(rate_arr[mask] * voxel_mass_g[mask]) / roi_mass_g)
        res = {"roi_mass_g": roi_mass_g, "roi_volume_ml": roi_volume_ml}
        if dose_rate:
            res["mean_dose_rate_gyh"] = mean_rate * 3600.0
        else:
            factor, teff = time_factor_s(roi_name)
            res["mean_dose_gy"] = mean_rate * factor
            if is_tia:
                tia_mbq_s = float(np.sum(image_arr[mask]) / 1e6)
                res["tia_mbq_s"] = tia_mbq_s
                res["tia_mbq_h"] = tia_mbq_s / 3600.0
            else:
                res["mean_dose_rate_gys"] = mean_rate
                res["roi_teff_h"] = float(teff)
        res["dose_engine"] = "map"
        results[roi_name] = res

        if verbose or logger.isEnabledFor(logging.DEBUG):
            log_msg = (
                f"\nDosimetry calculation log for ROI '{roi_name}' using {method} with a dose-rate map:\n"
                f"  - Radionuclide: {isotope}\n"
                f"  - Map unit: {unit}, map activity: "
                f"{'absolute' if dose_rate_map_activity_mbq is None else f'{dose_rate_map_activity_mbq:.3e} MBq'}\n"
                f"  - Image total: {image_total:.3e} {'Bq.s' if is_tia else 'Bq'} (map scale factor {scale:.4e})\n"
                f"  - ROI Volume: {roi_volume_ml:.3f} mL ({int(np.sum(mask))} voxels)\n"
                f"  - ROI Mass: {roi_mass_g:.3f} g\n"
                + (
                    f"  - Dose: mass-weighted mean of the scaled map\n"
                    if is_tia
                    else f"  - Mean dose rate at {time_from_injection_h:.3f} h: {mean_rate:.4e} Gy/s\n"
                )
                + "\n".join(f"  - {k}: {v}" for k, v in res.items() if k.startswith(("mean_dose", "roi_teff")))
            )
            if verbose:
                logger.info(log_msg)
            else:
                logger.debug(log_msg)

    if method_key == "ldm":
        return results, dose_image
    return results
