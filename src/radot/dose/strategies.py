"""ROI dose strategies (S-value methods) and the generic loop over the labels."""

import logging
from typing import List, Optional, Dict, Any
import numpy as np

from radot.db import get_physical_half_life_h
from radot.dose.models import (
    dose_hanscheid2017,
    dose_hanscheid2018,
    dose_madsen2018,
    dose_svalue,
)

logger = logging.getLogger("radot")


def strategy_hanscheid2017(
    image_arr: np.ndarray,
    mask: np.ndarray,
    voxel_volume_ml: float,
    roi_name: str,
    time_from_injection_h: float,
    isotope: str,
    roi_mass_g: float,
    roi_volume_ml: float,
    method_params: Optional[Dict[str, Any]],
    verbose: bool,
    ct_arr: Optional[np.ndarray] = None,
) -> tuple[float, float]:
    """Strategy implementation for Hanscheid 2017 method."""
    # Teff parameter resolution
    teff = None
    if method_params and roi_name in method_params:
        teff = method_params[roi_name].get("teff")

    if teff is None:
        try:
            teff = get_physical_half_life_h(isotope)
        except ValueError as e:
            raise ValueError(
                f"No effective half-life (teff) provided for ROI '{roi_name}' and {e}"
            )
        logger.warning(
            f"No Teff specified for '{roi_name}'. Using physical half-life of {isotope} ({teff:.2f} h) as fallback."
        )

    dose = dose_hanscheid2017(
        spect_bq_a=image_arr,
        roi_a=mask.astype(np.uint8),
        time_from_injection_h=time_from_injection_h,
        volume_voxel_m_l=voxel_volume_ml,
        effective_time_h=teff,
    )

    if verbose or logger.isEnabledFor(logging.DEBUG):
        voxels_in_roi = np.sum(mask)
        total_activity_bq = np.sum(image_arr[mask])
        total_activity_mbq = total_activity_bq / 1e6
        mean_activity_concentration_bq_ml = np.mean(image_arr[mask] / voxel_volume_ml)
        mean_activity_concentration_mbq_ml = mean_activity_concentration_bq_ml / 1e6
        density_str = (
            f"computed using companion CT, mean density: {roi_mass_g / roi_volume_ml:.3f} g/mL"
            if ct_arr is not None
            else "no companion CT, assuming water density 1.0 g/mL"
        )
        formula_str = f"0.125 * {mean_activity_concentration_mbq_ml:.3e} * 2^({time_from_injection_h:.2f} / {teff:.2f}) * {teff:.2f}"

        log_msg = (
            f"\nDosimetry calculation log for ROI '{roi_name}' using Hanscheid 2017:\n"
            f"  - Radionuclide: {isotope}\n"
            f"  - Scaling factor: 0.125\n"
            f"  - ROI Volume: {roi_volume_ml:.3f} mL ({voxels_in_roi} voxels)\n"
            f"  - ROI Mass: {roi_mass_g:.3f} g\n"
            f"    ({density_str})\n"
            f"  - Total Activity in ROI: {total_activity_mbq:.3e} MBq\n"
            f"  - Mean Activity Concentration: {mean_activity_concentration_mbq_ml:.3e} MBq/mL ({mean_activity_concentration_bq_ml:.3e} Bq/mL)\n"
            f"  - Time from injection: {time_from_injection_h:.3f} h\n"
            f"  - Effective half-life (Teff): {teff:.3f} h\n"
            f"  - Formula: {formula_str}\n"
            f"  - Voxel-mean Absorbed Dose: {dose:.6f} Gy"
        )
        if verbose:
            logger.info(log_msg)
        else:
            logger.debug(log_msg)

    return float(dose), float(teff)


def strategy_hanscheid2018(
    image_arr: np.ndarray,
    mask: np.ndarray,
    voxel_volume_ml: float,
    roi_name: str,
    time_from_injection_h: float,
    isotope: str,
    roi_mass_g: float,
    roi_volume_ml: float,
    method_params: Optional[Dict[str, Any]],
    verbose: bool,
    ct_arr: Optional[np.ndarray] = None,
) -> tuple[float, float]:
    """Strategy implementation for Hanscheid 2018 method."""
    from radot.db import resolve_organ_name, get_s_value, get_reference_mass

    phantom = "male"
    teff = None
    if method_params:
        phantom = method_params.get("phantom", "male")
        if roi_name in method_params:
            teff = method_params[roi_name].get("teff")

    if teff is None:
        teff = 2.0 * time_from_injection_h

    # Standardize phantom/organ/rad
    opendose_organ = resolve_organ_name(roi_name, phantom)
    s_value = get_s_value(phantom, isotope, opendose_organ, opendose_organ)
    s_mass = get_reference_mass(phantom, opendose_organ)
    mass_scaling = roi_mass_g / s_mass

    dose = dose_hanscheid2018(
        spect_bq_a=image_arr,
        roi_a=mask.astype(np.uint8),
        time_from_injection_h=time_from_injection_h,
        s_value=s_value,
        mass_scaling=mass_scaling,
        effective_time_h=teff,
    )

    if verbose or logger.isEnabledFor(logging.DEBUG):
        voxels_in_roi = np.sum(mask)
        total_activity_bq = np.sum(image_arr[mask])
        total_activity_mbq = total_activity_bq / 1e6
        density_str = (
            f"computed using companion CT, mean density: {roi_mass_g / roi_volume_ml:.3f} g/mL"
            if ct_arr is not None
            else "no companion CT, assuming water density 1.0 g/mL"
        )
        teff_is_custom = method_params is not None and roi_name in method_params and "teff" in method_params[roi_name]
        teff_str = f"{teff:.3f} h" if teff_is_custom else f"{teff:.3f} h (2*t)"
        formula_str = (
            "At * (Teff * 3600) / ln(2) * S / mass_scaling / 1000"
            if teff_is_custom
            else "At * (2 * t * 3600) / ln(2) * S / mass_scaling / 1000"
        )

        log_msg = (
            f"\nDosimetry calculation log for ROI '{roi_name}' using Hanscheid 2018:\n"
            f"  - Radionuclide: {isotope}\n"
            f"  - Phantom model: {phantom}\n"
            f"  - OpenDose Organ Name: {opendose_organ}\n"
            f"  - Reference Organ Mass: {s_mass:.3f} g\n"
            f"  - ROI Volume: {roi_volume_ml:.3f} mL ({voxels_in_roi} voxels)\n"
            f"  - ROI Mass: {roi_mass_g:.3f} g\n"
            f"    ({density_str})\n"
            f"  - Mass scaling factor: {mass_scaling:.4f}\n"
            f"  - S-value (self): {s_value:.3e} mGy/MBq.s\n"
            f"  - Total Activity in ROI: {total_activity_mbq:.3e} MBq\n"
            f"  - Time from injection: {time_from_injection_h:.3f} h\n"
            f"  - Effective half-life (Teff): {teff_str}\n"
            f"  - Formula: {formula_str}\n"
            f"  - Voxel-mean Absorbed Dose: {dose:.6f} Gy"
        )
        if verbose:
            logger.info(log_msg)
        else:
            logger.debug(log_msg)

    return float(dose), float(teff)


def strategy_madsen2018(
    image_arr: np.ndarray,
    mask: np.ndarray,
    voxel_volume_ml: float,
    roi_name: str,
    time_from_injection_h: float,
    isotope: str,
    roi_mass_g: float,
    roi_volume_ml: float,
    method_params: Optional[Dict[str, Any]],
    verbose: bool,
    ct_arr: Optional[np.ndarray] = None,
) -> tuple[float, float]:
    """Strategy implementation for Madsen 2018 method."""
    from radot.db import resolve_organ_name, get_s_value, get_reference_mass

    phantom = "male"
    teff = None
    if method_params:
        phantom = method_params.get("phantom", "male")
        if roi_name in method_params:
            teff = method_params[roi_name].get("teff")

    if teff is None:
        raise ValueError(
            f"Effective half-life (teff) must be provided for ROI '{roi_name}' when using the madsen2018 method."
        )

    # Standardize phantom/organ/rad
    opendose_organ = resolve_organ_name(roi_name, phantom)
    s_value = get_s_value(phantom, isotope, opendose_organ, opendose_organ)
    s_mass = get_reference_mass(phantom, opendose_organ)
    mass_scaling = roi_mass_g / s_mass

    dose = dose_madsen2018(
        spect_bq_a=image_arr,
        roi_a=mask.astype(np.uint8),
        time_from_injection_h=time_from_injection_h,
        s_value=s_value,
        mass_scaling=mass_scaling,
        effective_time_h=teff,
    )

    if verbose or logger.isEnabledFor(logging.DEBUG):
        voxels_in_roi = np.sum(mask)
        total_activity_bq = np.sum(image_arr[mask])
        total_activity_mbq = total_activity_bq / 1e6
        density_str = (
            f"computed using companion CT, mean density: {roi_mass_g / roi_volume_ml:.3f} g/mL"
            if ct_arr is not None
            else "no companion CT, assuming water density 1.0 g/mL"
        )
        formula_str = f"At * exp(k * t * 3600) / k * S / mass_scaling / 1000"

        log_msg = (
            f"\nDosimetry calculation log for ROI '{roi_name}' using Madsen 2018:\n"
            f"  - Radionuclide: {isotope}\n"
            f"  - Phantom model: {phantom}\n"
            f"  - OpenDose Organ Name: {opendose_organ}\n"
            f"  - Reference Organ Mass: {s_mass:.3f} g\n"
            f"  - ROI Volume: {roi_volume_ml:.3f} mL ({voxels_in_roi} voxels)\n"
            f"  - ROI Mass: {roi_mass_g:.3f} g\n"
            f"    ({density_str})\n"
            f"  - Mass scaling factor: {mass_scaling:.4f}\n"
            f"  - S-value (self): {s_value:.3e} mGy/MBq.s\n"
            f"  - Total Activity in ROI: {total_activity_mbq:.3e} MBq\n"
            f"  - Time from injection: {time_from_injection_h:.3f} h\n"
            f"  - Effective half-life (Teff): {teff:.3f} h\n"
            f"  - Formula: {formula_str}\n"
            f"  - Voxel-mean Absorbed Dose: {dose:.6f} Gy"
        )
        if verbose:
            logger.info(log_msg)
        else:
            logger.debug(log_msg)

    return float(dose), float(teff)


def strategy_svalue(
    image_arr: np.ndarray,
    mask: np.ndarray,
    voxel_volume_ml: float,
    roi_name: str,
    time_from_injection_h: Optional[float],
    isotope: str,
    roi_mass_g: float,
    roi_volume_ml: float,
    method_params: Optional[Dict[str, Any]],
    verbose: bool,
    ct_arr: Optional[np.ndarray] = None,
) -> tuple[float, None]:
    """Strategy implementation for the S-value method on a TIA image (Bq.s per voxel)."""
    from radot.db import resolve_organ_name, get_s_value, get_reference_mass

    phantom = "male"
    if method_params:
        phantom = method_params.get("phantom", "male")

    # Standardize phantom/organ/rad
    opendose_organ = resolve_organ_name(roi_name, phantom)
    s_value = get_s_value(phantom, isotope, opendose_organ, opendose_organ)
    s_mass = get_reference_mass(phantom, opendose_organ)
    mass_scaling = roi_mass_g / s_mass

    dose = dose_svalue(
        tia_bqs_a=image_arr,
        roi_a=mask.astype(np.uint8),
        s_value=s_value,
        mass_scaling=mass_scaling,
    )

    if verbose or logger.isEnabledFor(logging.DEBUG):
        voxels_in_roi = np.sum(mask)
        tia_mbq_s = np.sum(image_arr[mask]) / 1e6
        density_str = (
            f"computed using companion CT, mean density: {roi_mass_g / roi_volume_ml:.3f} g/mL"
            if ct_arr is not None
            else "no companion CT, assuming water density 1.0 g/mL"
        )

        log_msg = (
            f"\nDosimetry calculation log for ROI '{roi_name}' using S-value (TIA input):\n"
            f"  - Radionuclide: {isotope}\n"
            f"  - Phantom model: {phantom}\n"
            f"  - OpenDose Organ Name: {opendose_organ}\n"
            f"  - Reference Organ Mass: {s_mass:.3f} g\n"
            f"  - ROI Volume: {roi_volume_ml:.3f} mL ({voxels_in_roi} voxels)\n"
            f"  - ROI Mass: {roi_mass_g:.3f} g\n"
            f"    ({density_str})\n"
            f"  - Mass scaling factor: {mass_scaling:.4f}\n"
            f"  - S-value (self): {s_value:.3e} mGy/MBq.s\n"
            f"  - Time-Integrated Activity (TIA): {tia_mbq_s / 3600.0:.2f} MBq.h ({tia_mbq_s:.3e} MBq.s)\n"
            f"  - Formula: TIA * S / mass_scaling / 1000\n"
            f"  - Voxel-mean Absorbed Dose: {dose:.6f} Gy"
        )
        if verbose:
            logger.info(log_msg)
        else:
            logger.debug(log_msg)

    return float(dose), None


DOSIMETRY_STRATEGIES = {
    "hanscheid2017": strategy_hanscheid2017,
    "hanscheid2018": strategy_hanscheid2018,
    "madsen2018": strategy_madsen2018,
    "svalue": strategy_svalue,
}

# Units of time-integrated activity (TIA) images: Bq.s per voxel, Bq.s per mL
TIA_UNITS = ("bqs", "bqsml")
# Methods accepting a TIA image as input (they need no time point nor Teff)
TIA_METHODS = ("svalue", "ldm")


def compute_dose_generic(
    image_arr: np.ndarray,
    labels_arr: np.ndarray,
    isotope: str,
    voxel_volume_ml: float,
    time_from_injection_h: float,
    labels_map: Dict[int, str],
    ct_arr: Optional[np.ndarray] = None,
    rois: Optional[List[str]] = None,
    method_params: Optional[Dict[str, Any]] = None,
    verbose: bool = False,
    method: str = "hanscheid2017",
) -> dict:
    """
    Generic loop to calculate dosimetry using strategy functions.
    """
    method_key = method.lower()
    strategy = DOSIMETRY_STRATEGIES.get(method_key)
    if not strategy:
        raise ValueError(f"Unknown dosimetry method: {method}")

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

        # Common calculations: volume and mass
        roi_volume_ml = np.sum(mask) * voxel_volume_ml

        if ct_arr is not None:
            density = np.maximum(0.0, 1.0 + (ct_arr[mask] / 1000.0))
            roi_mass_g = np.sum(density) * voxel_volume_ml
        else:
            roi_mass_g = roi_volume_ml

        # Delegate calculation and logging to the strategy
        try:
            dose, teff = strategy(
                image_arr=image_arr,
                mask=mask,
                voxel_volume_ml=voxel_volume_ml,
                roi_name=roi_name,
                time_from_injection_h=time_from_injection_h,
                isotope=isotope,
                roi_mass_g=roi_mass_g,
                roi_volume_ml=roi_volume_ml,
                method_params=method_params,
                verbose=verbose,
                ct_arr=ct_arr,
            )
            results[roi_name] = {
                "mean_dose_gy": float(dose),
                "roi_mass_g": float(roi_mass_g),
                "roi_volume_ml": float(roi_volume_ml),
            }
            if teff is not None:
                results[roi_name]["roi_teff_h"] = float(teff)
            else:
                # TIA input: report the ROI time-integrated activity
                tia_mbq_s = float(np.sum(image_arr[mask]) / 1e6)
                results[roi_name]["tia_mbq_s"] = tia_mbq_s
                results[roi_name]["tia_mbq_h"] = tia_mbq_s / 3600.0
        except ValueError as e:
            if "not found in" in str(e) and (
                "S-value matrix" in str(e) or "phantom" in str(e)
            ):
                if rois is None:
                    logger.warning(f"Skipping dosimetry for ROI '{roi_name}': {e}")
                else:
                    raise ValueError(
                        f"Dosimetry Halt: {e}\n"
                        "Note: OpenDose lacks a single 'Skeleton' compartment and instead requires specific cortical/spongiosa targets.\n"
                        "To proceed, you must explicitly provide an external S-value override file via the CLI parameters (e.g., using --custom-s-values override.json) or adjust your input sidecar JSON mapping definitions to explicitly match target OpenDose bone/tissue definitions."
                    )
            else:
                raise e

    return results
