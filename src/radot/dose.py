import json
import logging
from typing import List, Optional, Dict, Any, Union, Tuple
import numpy as np
import SimpleITK as sitk
from pathlib import Path
from datetime import datetime

logger = logging.getLogger("radot")
from radot.db import get_physical_half_life_h
from radot.meta import get_time_from_injection_h


def dose_ldm(
    image_arr: np.ndarray,
    voxel_volume_ml: float,
    time_from_injection_h: float,
    isotope: str,
    ct_arr: Optional[np.ndarray] = None,
    hu_threshold: Optional[float] = -800.0,
    dose_rate: bool = False,
    cascade: bool = False,
) -> np.ndarray:
    """
    Voxel-level physical Local Deposition Model (LDM).
    - image_arr is physical activity in Bq per voxel
    - returns 3D float numpy array of absorbed dose (in Gy) or dose rate (in Gy/h) per voxel
    """
    from radot.db import get_cascade_energy_constant_mev

    # Get Delta_local in MeV (accounting for decay cascade if requested)
    delta_local_mev, _ = get_cascade_energy_constant_mev(isotope, cascade=cascade)

    # Convert MeV to Joules: 1 MeV = 1.60218e-13 Joules
    # Since mass will be in grams, we scale by 1000 to convert Joules/g to Gy (Joules/kg):
    # 1.60218e-13 * 1000 = 1.60218e-10
    mev_to_joules_per_g = 1.60218e-10

    # Compute density map (g/cm3 = g/mL) from Hounsfield Units (HU)
    if ct_arr is not None:
        density_arr = 1.0 + (ct_arr / 1000.0)
        # Ensure density is positive
        density_arr = np.clip(density_arr, 1e-3, None)
    else:
        density_arr = np.ones_like(image_arr)

    # Voxel mass in grams: Volume (mL) * Density (g/mL)
    voxel_mass_g = voxel_volume_ml * density_arr

    if dose_rate:
        # Instantaneous Dose Rate (Gy/h) at the time of acquisition:
        # activity (Bq) * Delta_local (MeV) * conversions (J/g/MeV) * 3600 (s/h) / mass (g)
        dose_arr = (
            image_arr * delta_local_mev * mev_to_joules_per_g * 3600.0
        ) / voxel_mass_g
    else:
        # Physical decay constant lambda = ln(2) / T_p (in hours)
        t_p_h = get_physical_half_life_h(isotope)
        decay_constant_h = np.log(2) / t_p_h

        # Compute Time Integrated Activity (TIA) per voxel in Bq.s:
        decay_constant_s = decay_constant_h / 3600.0
        tia_arr = (
            image_arr
            * np.exp(decay_constant_h * time_from_injection_h)
            / decay_constant_s
        )

        # Absorbed Dose (Gy) per voxel: TIA (Bq.s) * Delta_local (MeV) * conversions / mass (g)
        dose_arr = (tia_arr * delta_local_mev * mev_to_joules_per_g) / voxel_mass_g

    # Apply HU threshold if CT is present
    if ct_arr is not None and hu_threshold is not None:
        dose_arr[ct_arr < hu_threshold] = 0.0

    return dose_arr


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


def dose_hanscheid2018(
    spect_bq_a: np.ndarray,
    roi_a: np.ndarray,
    time_from_injection_h: float,
    s_value: float,
    mass_scaling: float,
    effective_time_h: Optional[float] = None,
    source_roi_a: Optional[np.ndarray] = None,
) -> float:
    """
    Input image and ROI must be numpy arrays
    - spect must be in Bq (not concentration)
    - acquisition time in hours
    - output is in Gray
    - effective_time_h: optional Teff in hours (defaults to 2 * t)
    - source_roi_a: optional mask for source organ activity (for cross-dose).
      If None, defaults to roi_a (self-dose).
    """
    # compute activity in source ROI (defaults to roi_a for self-dose)
    act_mask = (source_roi_a == 1) if source_roi_a is not None else (roi_a == 1)
    v = spect_bq_a[act_mask]
    At = np.sum(v) / 1e6

    teff = effective_time_h if effective_time_h is not None else (2.0 * time_from_injection_h)

    # S is in (mGy/MBq/s), so we get dose in mGy
    dose = (
        At
        * (teff * 3600.0)
        / np.log(2)
        * s_value
        / mass_scaling
        / 1000.0
    )

    return float(dose)


def dose_madsen2018(
    spect_bq_a: np.ndarray,
    roi_a: np.ndarray,
    time_from_injection_h: float,
    s_value: float,
    mass_scaling: float,
    effective_time_h: float,
    source_roi_a: Optional[np.ndarray] = None,
) -> float:
    """
    Input image and ROI must be numpy arrays
    - spect must be in Bq (not concentration)
    - output is in Gray
    - if source_roi_a is provided, activity is integrated in source_roi_a (source organ),
      computing cross-dose to target organ roi_a using s_value S(target <- source)
      and target mass_scaling.
    """
    # compute activity in the source ROI, in MBq
    source_mask = source_roi_a if source_roi_a is not None else roi_a
    v = spect_bq_a[source_mask == 1]
    At = np.sum(v) / 1e6

    # effective clearance rate in s^-1
    k = np.log(2) / (effective_time_h * 3600.0)

    # Integrated activity in MBq.s
    integrated_activity = At * np.exp(k * time_from_injection_h * 3600.0) / k

    # S is in (mGy/MBq/s), so we get dose in mGy. Scale and convert to Gy.
    dose = integrated_activity * s_value / mass_scaling / 1000.0

    return float(dose)


def fit_exp_linear(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Fit an exponential y = A * exp(k * x) via log-linear regression."""
    denegative = 1.0
    y_copy = np.array(y, dtype=float)
    if y_copy[0] < 0:
        denegative = -1.0
        y_copy[0] = -y_copy[0]
    y_copy = np.maximum(y_copy, 1e-12)
    log_y = np.log(y_copy)
    k, A_log = np.polyfit(x, log_y, 1)
    A = np.exp(A_log) * denegative
    if k > 0:
        k = 0.0
    return float(A), float(k)


def decay_corrected_tac(
    times: np.ndarray, activities: np.ndarray, decay_constant_h: float
) -> np.ndarray:
    """Decay-correct activities back to injection time (t=0)."""
    return activities * np.exp(decay_constant_h * times)


def triexpo_fit_jackson2013(
    times: Union[List[float], np.ndarray],
    activities: Union[List[float], np.ndarray],
    as_dict: bool = True,
) -> Union[Dict[str, float], np.ndarray]:
    """
    Fit a 3-component exponential biological clearance model from 3 time points.
    Based on Jackson et al., Med Phys. 2013 Nov;40(11):112503 (VRAK software).
    Times and decay-corrected activities must be ordered: t0 < t1 < t2.
    """
    times = np.asarray(times, dtype=float)
    activities = np.asarray(activities, dtype=float)
    if len(times) != 3 or len(activities) != 3:
        raise ValueError(
            f"triexpo_fit_jackson2013 requires exactly 3 timepoints and activities, got {len(times)} and {len(activities)}."
        )

    t0 = float(times[0])
    t1 = float(times[1])
    t3 = float(times[2])
    d0 = float(activities[0])
    d1 = float(activities[1])
    d3 = float(activities[2])

    if d3 == d1 or d1 == d0 or d0 == d3:
        d0 += 0.0001
        d1 += 0.0002
        d3 += 0.0003

    d1lin = ((d3 - d0) / (t3 - t0)) * t1 + (d0 - ((d3 - d0) / (t3 - t0) * t0))
    A3 = d3 / np.exp(-0.001 * t3)
    d1_lowslope = A3 * np.exp(-0.001 * t1)

    if d1 < d1_lowslope:
        if d0 < d3:
            if d1 < d1lin:
                d1 = d1lin
        else:
            d1 = d1_lowslope

    if d3 < 0.1 * d1:
        d3 = 0.1 * d1
    if d0 < 0.2 * d1:
        d0 = 0.2 * d1

    if d3 > d1:
        k3 = -0.001
        A3 = d3 / np.exp(k3 * t3)

        x = np.array([t1, t3])
        y = np.array([d1 - (A3 * np.exp(k3 * t1)), 0.01 * d3])
        A2, k2 = fit_exp_linear(x, y)
        if (A2 + A3) < 0:
            A2 = -A3
            A1 = 0.0
            k1 = 0.0
        else:
            A1 = -(A2 + A3)
            k1 = -1.3
    else:
        x = np.array([t1, t3])
        y = np.array([d1, d3])
        A3, k3 = fit_exp_linear(x, y)
        if k3 > -0.001:
            k3 = -0.001
            A3 = d3 / np.exp(k3 * t3)
        x = np.array([t0, t1])
        y = np.array([d0 - (A3 * np.exp(k3 * t0)), 0.01 * d1])
        A2, k2 = fit_exp_linear(x, y)
        if (A2 + A3) < 0:
            A2 = -A3
            A1 = 0.0
            k1 = 0.0
        else:
            A1 = -(A2 + A3)
            k1 = -1.3

    if as_dict:
        return {
            "A1": float(A1),
            "k1": float(k1),
            "A2": float(A2),
            "k2": float(k2),
            "A3": float(A3),
            "k3": float(k3),
        }
    return np.array([A1, k1, A2, k2, A3, k3], dtype=float)


def triexpo_integrate_jackson2013(
    params: Dict[str, float], decay_constant_h: float
) -> float:
    """
    Analytically integrate the physical time-activity curve:
    A(t) = sum_{i=1}^3 A_i * exp((k_i - lambda_p) * t) from t = 0 to infinity.
    Returns the time-integrated activity in activity_unit * seconds (scaled by 3600).
    """
    total_tia_h = 0.0
    for i in (1, 2, 3):
        a_i = params[f"A{i}"]
        k_i = params[f"k{i}"]
        eff_rate = decay_constant_h - k_i
        if eff_rate > 0:
            total_tia_h += a_i / eff_rate
    tia_s = total_tia_h * 3600.0
    return max(0.0, float(tia_s))


def dose_jackson2013(
    times_h: Union[List[float], np.ndarray],
    activities_mbq: Union[List[float], np.ndarray],
    s_value: float,
    mass_scaling: float,
    isotope: str,
) -> tuple[float, float, float, Dict[str, float]]:
    """
    Compute absorbed dose in Gy using Jackson 2013 tri-exponential MTP method.
    Returns:
        (dose_gy, tia_mbq_s, teff_terminal_h, fit_params)
    """
    times = np.asarray(times_h, dtype=float)
    acts = np.asarray(activities_mbq, dtype=float)

    # Sort timepoints in increasing order
    sort_idx = np.argsort(times)
    times = times[sort_idx]
    acts = acts[sort_idx]

    t_p_h = get_physical_half_life_h(isotope)
    lambda_p = np.log(2) / t_p_h

    # Decay-correct activities to t=0
    d_acts = decay_corrected_tac(times, acts, lambda_p)

    # Fit tri-exponential biological curve
    fit_params = triexpo_fit_jackson2013(times, d_acts, as_dict=True)

    # Analytically integrate to get TIA in MBq.s
    tia_mbq_s = triexpo_integrate_jackson2013(fit_params, lambda_p)

    # Dose in Gy
    dose_gy = (tia_mbq_s * s_value) / mass_scaling / 1000.0

    # Terminal effective half-life
    eff_terminal_rate = lambda_p - fit_params["k3"]
    teff_terminal_h = np.log(2) / eff_terminal_rate if eff_terminal_rate > 0 else t_p_h

    return float(dose_gy), float(tia_mbq_s), float(teff_terminal_h), fit_params


# Approximate maximum beta ranges in soft tissue (in mm)
ISOTOPE_BETA_RANGE_MM = {
    "lu177": 1.8,
    "y90": 11.0,
    "i131": 2.4,
    "tb161": 2.2,
    "re188": 10.0,
    "sm153": 3.1,
    "p32": 8.0,
    "cu67": 1.8,
}


def roi_dose_uncertainty_bounds(
    dose: Union[np.ndarray, sitk.Image],
    uncertainty: Union[np.ndarray, sitk.Image],
    mask: Union[np.ndarray, sitk.Image],
) -> tuple[float, float]:
    """
    Compute theoretical lower and upper bounds on the relative uncertainty of an ROI mean dose.

    Parameters:
        dose: 3D dose or dose rate map (SimpleITK Image or NumPy array).
        uncertainty: 3D relative uncertainty map (eps_i = sigma_i / D_i).
        mask: ROI binary mask (voxels > 0 are included).

    Returns:
        (unc_min, unc_max):
            - unc_min: Lower bound assuming independent voxels (Cov = 0).
                       unc_min = sqrt(sum((eps_i * D_i)^2)) / sum(D_i)
            - unc_max: Upper bound assuming 100% correlation across voxels.
                       unc_max = sum(eps_i * D_i) / sum(D_i)
    """
    d_arr = (
        sitk.GetArrayViewFromImage(dose)
        if isinstance(dose, sitk.Image)
        else np.asarray(dose)
    )
    u_arr = (
        sitk.GetArrayViewFromImage(uncertainty)
        if isinstance(uncertainty, sitk.Image)
        else np.asarray(uncertainty)
    )
    m_arr = (
        sitk.GetArrayViewFromImage(mask)
        if isinstance(mask, sitk.Image)
        else np.asarray(mask)
    )

    roi_filter = m_arr > 0
    if not np.any(roi_filter):
        raise ValueError("ROI mask contains no voxels.")

    d_roi = d_arr[roi_filter].astype(np.float64)
    u_roi = u_arr[roi_filter].astype(np.float64)

    total_dose = np.sum(d_roi)
    if total_dose <= 0:
        return 0.0, 0.0

    # Absolute standard deviation per voxel: sigma_i = eps_i * D_i
    sigma_vox = u_roi * d_roi

    # Lower bound: independent voxels, variance of sum = sum of variances
    unc_min = float(np.sqrt(np.sum(sigma_vox**2)) / total_dose)

    # Upper bound: perfectly correlated voxels, std of sum = sum of stds
    unc_max = float(np.sum(sigma_vox) / total_dose)

    return unc_min, unc_max


def roi_dose_uncertainty_estimate(
    dose: Union[np.ndarray, sitk.Image],
    uncertainty: Union[np.ndarray, sitk.Image],
    mask: Union[np.ndarray, sitk.Image],
    spacing_mm: Optional[tuple[float, float, float]] = None,
    isotope: Optional[str] = None,
    correlation_length_mm: Optional[float] = None,
) -> float:
    """
    Estimate realistic ROI relative uncertainty accounting for finite particle range.

    Estimates the number of voxels correlated by a single primary history (nu_cluster)
    using the particle range and voxel spacing:
        nu_cluster = prod(1 + range / spacing_d)
        unc_est = min(unc_max, unc_min * sqrt(nu_cluster))

    Parameters:
        dose: 3D dose or dose rate map (SimpleITK Image or NumPy array).
        uncertainty: 3D relative uncertainty map.
        mask: ROI binary mask.
        spacing_mm: Voxel spacing in mm (x, y, z). Auto-extracted if dose is sitk.Image.
        isotope: Radionuclide name (e.g. 'Lu177', 'Y90') to retrieve typical beta range.
        correlation_length_mm: Direct override for particle range / correlation length in mm.

    Returns:
        float: Estimated relative uncertainty as a fraction.
    """
    if spacing_mm is None and isinstance(dose, sitk.Image):
        spacing_mm = dose.GetSpacing()
    elif spacing_mm is None:
        spacing_mm = (1.0, 1.0, 1.0)

    if correlation_length_mm is None:
        iso_key = (isotope or "").lower().replace("-", "")
        correlation_length_mm = ISOTOPE_BETA_RANGE_MM.get(iso_key, 2.0)

    unc_min, unc_max = roi_dose_uncertainty_bounds(dose, uncertainty, mask)

    # Effective correlation cluster size across spatial dimensions
    nu_cluster = 1.0
    for s in spacing_mm:
        if s > 0:
            nu_cluster *= 1.0 + float(correlation_length_mm) / float(s)

    unc_est = float(min(unc_max, unc_min * np.sqrt(nu_cluster)))
    return unc_est


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


DOSIMETRY_STRATEGIES = {
    "hanscheid2017": strategy_hanscheid2017,
    "hanscheid2018": strategy_hanscheid2018,
    "madsen2018": strategy_madsen2018,
}


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
                "roi_teff_h": float(teff),
            }
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


def compute_roi_dose_mtp(
    images: List[Any],
    labels: Union[Any, List[Any]],
    method: str = "jackson2013",
    isotope: str = "Lu177",
    times_h: Optional[List[float]] = None,
    ct_image: Optional[Union[Any, List[Any]]] = None,
    rois: Optional[List[str]] = None,
    labels_map: Optional[Dict[int, str]] = None,
    image_paths: Optional[List[Optional[Path]]] = None,
    labels_path: Optional[Union[Path, List[Optional[Path]]]] = None,
    ct_path: Optional[Union[Path, List[Optional[Path]]]] = None,
    method_params: Optional[Dict[str, Any]] = None,
    verbose: bool = False,
) -> Dict[str, Any]:
    """
    Computes multiple timepoint (MTP) absorbed doses for specified regions of interest.
    Currently supports the Jackson 2013 tri-exponential method.
    Supports either a single segmentation (registered across all timepoints)
    or independent segmentations per timepoint.
    """
    from radot.stats import has_same_grid
    from radot.resampling import resample_image
    from radot.meta import (
        get_metadata_case_insensitive,
        get_time_from_injection_h,
        get_sidecar_path,
    )
    from radot.db import resolve_organ_name, get_s_value, get_reference_mass

    # 1. Resolve and read SimpleITK images
    itk_images = []
    resolved_paths = []
    for i, img in enumerate(images):
        if isinstance(img, (str, Path)):
            p = Path(img)
            resolved_paths.append(p)
            itk_images.append(sitk.ReadImage(str(p)))
        else:
            itk_images.append(img)
            if image_paths and i < len(image_paths) and image_paths[i]:
                resolved_paths.append(Path(image_paths[i]))
            else:
                resolved_paths.append(None)

    num_tp = len(itk_images)
    method_key = method.lower()
    if method_key == "jackson2013" and num_tp != 3:
        raise ValueError(
            f"Method 'jackson2013' requires exactly 3 timepoints, but received {num_tp}."
        )

    # 2. Resolve and read labels (can be list or single)
    if isinstance(labels, (list, tuple)):
        raw_labels = list(labels)
    else:
        raw_labels = [labels] * num_tp

    if isinstance(labels_path, (list, tuple)):
        raw_labels_paths = list(labels_path)
    else:
        raw_labels_paths = [labels_path] * num_tp

    itk_labels = []
    resolved_labels_paths = []
    for i, lbl in enumerate(raw_labels):
        p = raw_labels_paths[i] if i < len(raw_labels_paths) else None
        if isinstance(lbl, (str, Path)):
            p = Path(lbl)
            itk_labels.append(sitk.ReadImage(str(p)))
            resolved_labels_paths.append(p)
        elif lbl is not None:
            itk_labels.append(lbl)
            resolved_labels_paths.append(Path(p) if p else None)
        else:
            itk_labels.append(None)
            resolved_labels_paths.append(None)

    if any(l is None for l in itk_labels):
        raise ValueError("Labels segmentation image is required for all MTP timepoints.")

    # 3. Resolve and read companion CT (can be list or single)
    if isinstance(ct_image, (list, tuple)):
        raw_cts = list(ct_image)
    else:
        raw_cts = [ct_image] * num_tp

    if isinstance(ct_path, (list, tuple)):
        raw_ct_paths = list(ct_path)
    else:
        raw_ct_paths = [ct_path] * num_tp

    itk_cts = []
    resolved_ct_paths = []
    for i, ct in enumerate(raw_cts):
        p = raw_ct_paths[i] if i < len(raw_ct_paths) else None
        if ct is None and p is not None:
            ct = p
        if isinstance(ct, (str, Path)):
            p = Path(ct)
            itk_cts.append(sitk.ReadImage(str(p)))
            resolved_ct_paths.append(p)
        elif ct is not None:
            itk_cts.append(ct)
            resolved_ct_paths.append(Path(p) if p else None)
        else:
            itk_cts.append(None)
            resolved_ct_paths.append(None)

    # 4. Resolve acquisition times
    if times_h is None:
        times_h = []
        for i, (img, p) in enumerate(zip(itk_images, resolved_paths)):
            t = get_time_from_injection_h(img, p)
            if t is None:
                name_str = p.name if p else f"image #{i+1}"
                raise ValueError(
                    f"Time from injection could not be found for timepoint {i+1} ({name_str}). "
                    "Please ensure injection_date and acquisition_date are recorded in the image metadata sidecar."
                )
            times_h.append(float(t))
    else:
        times_h = [float(t) for t in times_h]

    # Sort timepoints chronologically
    times_arr = np.array(times_h, dtype=float)
    sort_idx = np.argsort(times_arr)
    itk_images = [itk_images[i] for i in sort_idx]
    resolved_paths = [resolved_paths[i] for i in sort_idx]
    itk_labels = [itk_labels[i] for i in sort_idx]
    resolved_labels_paths = [resolved_labels_paths[i] for i in sort_idx]
    itk_cts = [itk_cts[i] for i in sort_idx]
    resolved_ct_paths = [resolved_ct_paths[i] for i in sort_idx]
    times_h = [times_h[i] for i in sort_idx]

    # 5. Resolve phantom
    phantom = "male"
    if method_params and "phantom" in method_params:
        phantom = method_params["phantom"]

    # 6. Process each timepoint: alignment, mass, and activity extraction
    tp_data = []
    all_known_rois = set()

    for i in range(num_tp):
        img = itk_images[i]
        p_img = resolved_paths[i]
        lbl = itk_labels[i]
        p_lbl = resolved_labels_paths[i]
        ct = itk_cts[i]

        # Load labels mapping for this timepoint
        lbl_map_i = {}
        if p_lbl is not None:
            json_path = get_sidecar_path(p_lbl)
            if json_path.exists():
                try:
                    with open(json_path, "r") as f:
                        raw_map = json.load(f)
                        lbl_map_i = {int(k): v for k, v in raw_map.items()}
                except Exception:
                    pass
        if not lbl_map_i and labels_map:
            lbl_map_i = dict(labels_map)

        name_to_lbl = {v.lower(): k for k, v in lbl_map_i.items()}
        for r_name in lbl_map_i.values():
            all_known_rois.add(r_name)

        # Resample SPECT to labels grid if necessary
        unit = get_metadata_case_insensitive(img, "unit", p_img) or "bqml"
        unit = unit.lower()
        preserve_total = unit == "bq"

        if not has_same_grid(img, lbl):
            logger.info(
                f"Timepoint {i+1} image grid does not match its labels grid. Resampling to labels grid..."
            )
            img_aligned = resample_image(
                img,
                image_type="spect",
                preserve_total=preserve_total,
                gaussian=True,
                reference_image=lbl,
            )
        else:
            img_aligned = img

        # Resample CT to labels grid if needed
        ct_aligned = None
        if ct is not None:
            if not has_same_grid(ct, lbl):
                logger.info(
                    f"Timepoint {i+1} CT grid does not match its labels grid. Resampling CT to labels grid..."
                )
                ct_aligned = resample_image(
                    ct,
                    image_type="ct",
                    preserve_total=False,
                    gaussian=True,
                    reference_image=lbl,
                )
            else:
                ct_aligned = ct

        # Voxel volumes and masses
        sp = lbl.GetSpacing()
        vox_vol_ml = (sp[0] * sp[1] * sp[2]) / 1000.0
        lbl_arr = sitk.GetArrayFromImage(lbl)

        if ct_aligned is not None:
            ct_arr = sitk.GetArrayFromImage(ct_aligned)
            density_arr = np.maximum(0.0, 1.0 + (ct_arr / 1000.0))
        else:
            density_arr = np.ones_like(lbl_arr, dtype=float)

        vox_mass_g = vox_vol_ml * density_arr

        # Activity array in MBq
        img_arr = sitk.GetArrayFromImage(img_aligned).astype(np.float64)
        if unit == "bqml":
            mbq_arr = (img_arr * vox_vol_ml) / 1e6
        elif unit == "bq":
            mbq_arr = img_arr / 1e6
        else:
            raise ValueError(
                f"Unsupported SPECT unit '{unit}' for timepoint {i+1}. Expected 'bqml' or 'bq'."
            )

        tp_data.append({
            "lbl_arr": lbl_arr,
            "vox_vol_ml": vox_vol_ml,
            "vox_mass_g": vox_mass_g,
            "mbq_arr": mbq_arr,
            "name_to_lbl": name_to_lbl,
            "lbl_map": lbl_map_i,
        })

    # 7. Compute dose per ROI
    target_rois = list(rois) if rois is not None else sorted(list(all_known_rois))
    results = {}

    for roi_name in target_rois:
        roi_acts_mbq = []
        roi_masses_g = []
        roi_volumes_ml = []

        for tp in tp_data:
            lbl_id = tp["name_to_lbl"].get(roi_name.lower())
            if lbl_id is not None:
                mask = tp["lbl_arr"] == lbl_id
                act = float(np.sum(tp["mbq_arr"][mask]))
                mass = float(np.sum(tp["vox_mass_g"][mask]))
                vol = float(np.sum(mask) * tp["vox_vol_ml"])
            else:
                act, mass, vol = 0.0, 0.0, 0.0
                try:
                    num_id = int(roi_name)
                    if num_id in tp["lbl_map"]:
                        mask = tp["lbl_arr"] == num_id
                        act = float(np.sum(tp["mbq_arr"][mask]))
                        mass = float(np.sum(tp["vox_mass_g"][mask]))
                        vol = float(np.sum(mask) * tp["vox_vol_ml"])
                except ValueError:
                    pass

            roi_acts_mbq.append(act)
            roi_masses_g.append(mass)
            roi_volumes_ml.append(vol)

        # Skip if organ was never segmented / has zero volume
        if not any(v > 0 for v in roi_volumes_ml):
            continue

        # Baseline organ mass and volume (from TP1 if non-zero, else first non-zero)
        non_zero_masses = [m for m in roi_masses_g if m > 0]
        roi_mass_g = roi_masses_g[0] if roi_masses_g[0] > 0 else non_zero_masses[0]
        non_zero_vols = [v for v in roi_volumes_ml if v > 0]
        roi_volume_ml = roi_volumes_ml[0] if roi_volumes_ml[0] > 0 else non_zero_vols[0]

        # Resolve OpenDose S-value and mass scaling
        opendose_organ = resolve_organ_name(roi_name, phantom)
        s_val = get_s_value(phantom, isotope, opendose_organ, opendose_organ)
        ref_mass = get_reference_mass(phantom, opendose_organ)
        mass_scaling = roi_mass_g / ref_mass

        # Compute dose using Jackson 2013
        dose_gy, tia_mbq_s, teff_h, fit_params = dose_jackson2013(
            times_h=times_h,
            activities_mbq=roi_acts_mbq,
            s_value=s_val,
            mass_scaling=mass_scaling,
            isotope=isotope,
        )

        results[roi_name] = {
            "mean_dose_gy": float(dose_gy),
            "roi_mass_g": float(roi_mass_g),
            "roi_volume_ml": float(roi_volume_ml),
            "roi_teff_h": float(teff_h),
            "tia_mbq_s": float(tia_mbq_s),
            "tia_mbq_h": float(tia_mbq_s / 3600.0),
            "activities_mbq": roi_acts_mbq,
            "times_h": times_h,
            "fit_params": fit_params,
        }

        if verbose or logger.isEnabledFor(logging.DEBUG):
            log_msg = (
                f"\nDosimetry calculation log for ROI '{roi_name}' using Jackson 2013 (MTP):\n"
                f"  - Radionuclide: {isotope}\n"
                f"  - Phantom model: {phantom}\n"
                f"  - OpenDose Organ Name: {opendose_organ}\n"
                f"  - Reference Organ Mass: {ref_mass:.3f} g\n"
                f"  - ROI Volume (baseline): {roi_volume_ml:.3f} mL\n"
                f"  - ROI Mass (baseline): {roi_mass_g:.3f} g (mass scaling: {mass_scaling:.4f})\n"
                f"  - S-value (self): {s_val:.3e} mGy/MBq.s\n"
                f"  - Timepoints (h): {[round(t, 2) for t in times_h]}\n"
                f"  - Activities (MBq): {[round(a, 2) for a in roi_acts_mbq]}\n"
                f"  - Fitted tri-exponential: A1={fit_params['A1']:.1f}, k1={fit_params['k1']:.3f}, "
                f"A2={fit_params['A2']:.1f}, k2={fit_params['k2']:.3f}, "
                f"A3={fit_params['A3']:.1f}, k3={fit_params['k3']:.3f}\n"
                f"  - Time-Integrated Activity (TIA): {tia_mbq_s / 3600.0:.2f} MBq.h ({tia_mbq_s:.3e} MBq.s)\n"
                f"  - Terminal Teff: {teff_h:.1f} h\n"
                f"  - Absorbed Dose: {dose_gy:.4f} Gy"
            )
            if verbose:
                logger.info(log_msg)
            else:
                logger.debug(log_msg)

    return results


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
) -> Any:
    """
    Computes average/voxel absorbed doses for specified regions of interest.
    """
    # Check if multiple images (MTP) were passed
    if isinstance(image, (list, tuple)):
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

    # 1. Fast-Fail Boundary Checks
    method_key = method.lower()
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

    # 3. Timing Extraction
    time_from_injection_h = get_time_from_injection_h(image, image_path)

    # 4. Calibrated Activity Scaling
    unit = get_metadata_case_insensitive(image, "unit", image_path)
    if unit is None:
        raise ValueError("Image unit is not specified and was not found in metadata.")

    if unit.lower() != "bq":
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


def load_dose_config(
    config: Union[str, Path, Dict[str, Any]],
) -> Dict[str, Any]:
    """
    Load and normalize a dosimetry configuration file or dictionary.
    Resolves relative file paths relative to the configuration file location.
    """
    if isinstance(config, (str, Path)):
        config_path = Path(config).resolve()
        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_path}")
        with open(config_path, "r") as f:
            raw_data = json.load(f)
        base_dir = config_path.parent
    elif isinstance(config, dict):
        raw_data = dict(config)
        base_dir = Path.cwd()
    else:
        raise TypeError(f"Expected path or dict for config, got {type(config).__name__}")

    def resolve_path(p: Optional[Union[str, Path]]) -> Optional[Path]:
        if p is None:
            return None
        p_obj = Path(p)
        if not p_obj.is_absolute():
            return (base_dir / p_obj).resolve()
        return p_obj.resolve()

    res: Dict[str, Any] = {}

    # Activity Image path
    img = raw_data.get("image_path") or raw_data.get("image")
    if img is not None:
        if isinstance(img, (list, tuple)):
            res["timepoints"] = [
                resolve_path(p) if isinstance(p, (str, Path)) else p for p in img
            ]
        else:
            res["image_path"] = resolve_path(img)

    # Timepoints (MTP)
    timepoints = raw_data.get("timepoints") or raw_data.get("images")
    if timepoints is not None and isinstance(timepoints, (list, tuple)):
        resolved_tps = []
        for tp in timepoints:
            if isinstance(tp, (str, Path)):
                resolved_tps.append(resolve_path(tp))
            elif isinstance(tp, dict):
                tp_dict = dict(tp)
                for k in ["image", "image_path", "path"]:
                    if k in tp_dict and tp_dict[k] is not None:
                        tp_dict[k] = resolve_path(tp_dict[k])
                for k in ["labels", "labels_path"]:
                    if k in tp_dict and tp_dict[k] is not None:
                        tp_dict[k] = resolve_path(tp_dict[k])
                for k in ["ct", "ct_path"]:
                    if k in tp_dict and tp_dict[k] is not None:
                        tp_dict[k] = resolve_path(tp_dict[k])
                resolved_tps.append(tp_dict)
            else:
                resolved_tps.append(tp)
        res["timepoints"] = resolved_tps

    # Labels path
    lbl = raw_data.get("labels_path") or raw_data.get("labels")
    if lbl is not None:
        res["labels_path"] = resolve_path(lbl)

    # CT path
    ct = raw_data.get("ct_path") or raw_data.get("ct")
    if ct is not None:
        res["ct_path"] = resolve_path(ct)

    # Output path
    out = raw_data.get("output_path") or raw_data.get("output")
    if out is not None:
        res["output_path"] = resolve_path(out)

    # Method / Model
    method = raw_data.get("method") or raw_data.get("model")
    if method is not None:
        res["method"] = str(method)

    # Radionuclide / Isotope
    isotope = raw_data.get("isotope") or raw_data.get("radionuclide")
    if isotope is not None:
        res["isotope"] = str(isotope)

    # Total activity
    act = (
        raw_data.get("total_activity_mbq")
        if "total_activity_mbq" in raw_data
        else raw_data.get("total_activity")
    )
    if act is not None:
        res["total_activity_mbq"] = float(act)

    # Method params & Teff
    method_params: Dict[str, Any] = {}
    if "method_params" in raw_data and isinstance(raw_data["method_params"], dict):
        method_params = dict(raw_data["method_params"])

    # Phantom
    phantom = raw_data.get("phantom")
    if phantom is not None:
        res["phantom"] = str(phantom)
        if "phantom" not in method_params:
            method_params["phantom"] = str(phantom)

    # Resample target
    resample = raw_data.get("resample_target") or raw_data.get("resample")
    if resample is not None:
        res["resample_target"] = str(resample)

    # Cascade
    if "cascade" in raw_data:
        res["cascade"] = bool(raw_data["cascade"])

    # Dose rate
    if "dose_rate" in raw_data:
        res["dose_rate"] = bool(raw_data["dose_rate"])

    # HU threshold
    if "hu_threshold" in raw_data:
        res["hu_threshold"] = float(raw_data["hu_threshold"])

    # Verbose
    if "verbose" in raw_data:
        res["verbose"] = bool(raw_data["verbose"])

    # Quiet
    if "quiet" in raw_data:
        res["quiet"] = bool(raw_data["quiet"])

    # ROIs
    rois = raw_data.get("rois")
    if rois is not None:
        if isinstance(rois, str):
            rois_list = [r.strip() for r in rois.split(",") if r.strip()]
        elif isinstance(rois, dict):
            rois_list = []
            for roi_name, roi_val in rois.items():
                rois_list.append(str(roi_name))
                if isinstance(roi_val, dict):
                    if roi_name not in method_params:
                        method_params[roi_name] = {}
                    method_params[roi_name].update(roi_val)
                elif roi_val is not None:
                    try:
                        if roi_name not in method_params:
                            method_params[roi_name] = {}
                        method_params[roi_name]["teff"] = float(roi_val)
                    except (ValueError, TypeError):
                        pass
        elif isinstance(rois, (list, tuple)):
            rois_list = []
            for r in rois:
                if isinstance(r, str):
                    for item in r.split(","):
                        item_stripped = item.strip()
                        if item_stripped:
                            rois_list.append(item_stripped)
                elif isinstance(r, dict):
                    name = r.get("name") or r.get("roi")
                    if name:
                        rois_list.append(str(name))
                        params = {k: v for k, v in r.items() if k not in ("name", "roi")}
                        if params:
                            if name not in method_params:
                                method_params[name] = {}
                            method_params[name].update(params)
                else:
                    rois_list.append(str(r))
        else:
            rois_list = list(rois)
        res["rois"] = rois_list if rois_list else None

    teff = raw_data.get("teff")
    if teff is not None:
        if isinstance(teff, dict):
            for roi_name, val in teff.items():
                if roi_name not in method_params:
                    method_params[roi_name] = {}
                if isinstance(val, dict):
                    method_params[roi_name].update(val)
                else:
                    method_params[roi_name]["teff"] = float(val)
        elif isinstance(teff, (list, tuple)):
            for item in teff:
                if isinstance(item, str):
                    for t in item.split(","):
                        t_strip = t.strip()
                        if not t_strip:
                            continue
                        parts = t_strip.split(":")
                        if len(parts) == 2:
                            r_name = parts[0].strip()
                            if r_name not in method_params:
                                method_params[r_name] = {}
                            method_params[r_name]["teff"] = float(parts[1].strip())

    if method_params:
        res["method_params"] = method_params

    return res


def compute_roi_dose_file(
    image_path: Optional[Union[str, Path, List[Any]]] = None,
    labels_path: Optional[Union[str, Path]] = None,
    method: Optional[str] = None,
    isotope: Optional[str] = None,
    ct_path: Optional[Union[str, Path]] = None,
    total_activity_mbq: Optional[float] = None,
    rois: Optional[List[str]] = None,
    method_params: Optional[Dict[str, Any]] = None,
    output_path: Optional[Union[str, Path]] = None,
    resample_target: Optional[str] = None,
    verbose: bool = False,
    hu_threshold: Optional[float] = -800.0,
    dose_rate: bool = False,
    cascade: bool = False,
    config: Optional[Union[str, Path, Dict[str, Any]]] = None,
    timepoints: Optional[List[Any]] = None,
) -> Any:
    # 0. Check if image_path is actually a JSON config file
    if (
        image_path is not None
        and isinstance(image_path, (str, Path))
        and str(image_path).lower().endswith(".json")
        and config is None
        and method is None
        and isotope is None
    ):
        config = image_path
        image_path = None

    # Merge config file if provided
    if config is not None:
        cfg = load_dose_config(config)
        if timepoints is None:
            timepoints = cfg.get("timepoints")
        if image_path is None:
            image_path = cfg.get("image_path")
        if labels_path is None:
            labels_path = cfg.get("labels_path")
        if ct_path is None:
            ct_path = cfg.get("ct_path")
        if method is None:
            method = cfg.get("method")
        if isotope is None:
            isotope = cfg.get("isotope")
        if total_activity_mbq is None:
            total_activity_mbq = cfg.get("total_activity_mbq")
        if rois is None:
            rois = cfg.get("rois")
        if output_path is None:
            output_path = cfg.get("output_path")
        if resample_target is None:
            resample_target = cfg.get("resample_target")
        if not verbose and cfg.get("verbose"):
            verbose = True
        if hu_threshold == -800.0 and "hu_threshold" in cfg:
            hu_threshold = cfg["hu_threshold"]
        if not dose_rate and cfg.get("dose_rate"):
            dose_rate = True
        if not cascade and cfg.get("cascade"):
            cascade = True

        cfg_mp = cfg.get("method_params", {})
        if method_params is not None:
            cfg_mp.update(method_params)
        method_params = cfg_mp

    # Multi-timepoint (MTP) dispatch
    if timepoints is not None or isinstance(image_path, (list, tuple)):
        tp_list = timepoints if timepoints is not None else image_path
        extracted_paths = []
        extracted_labels = []
        extracted_cts = []
        extracted_times = []
        has_times = False
        for tp in tp_list:
            if isinstance(tp, (str, Path)):
                extracted_paths.append(Path(tp))
                extracted_labels.append(None)
                extracted_cts.append(None)
            elif isinstance(tp, dict):
                p = tp.get("image_path") or tp.get("image") or tp.get("path")
                if p:
                    extracted_paths.append(Path(p))
                lbl_tp = tp.get("labels_path") or tp.get("labels")
                extracted_labels.append(Path(lbl_tp) if lbl_tp else None)
                ct_tp = tp.get("ct_path") or tp.get("ct")
                extracted_cts.append(Path(ct_tp) if ct_tp else None)
                t = tp.get("time_h") if "time_h" in tp else tp.get("time")
                if t is not None:
                    has_times = True
                    extracted_times.append(float(t))
                else:
                    extracted_times.append(None)
            else:
                extracted_paths.append(tp)
                extracted_labels.append(None)
                extracted_cts.append(None)

        times_h = (
            extracted_times
            if (has_times and all(t is not None for t in extracted_times))
            else None
        )

        # Merge per-timepoint and global labels / CT:
        # If both are provided (global labels and one per tp, the one per tp is used)
        resolved_labels = [
            (lbl if lbl is not None else (Path(labels_path) if labels_path else None))
            for lbl in extracted_labels
        ]
        resolved_cts = [
            (c if c is not None else (Path(ct_path) if ct_path else None))
            for c in extracted_cts
        ]

        if not method:
            raise ValueError("Dosimetry method ('method') must be provided.")
        if not isotope:
            raise ValueError("Radionuclide ('isotope' / 'radionuclide') must be provided.")

        for idx, lbl in enumerate(resolved_labels):
            if lbl is None:
                raise ValueError(
                    f"Labels segmentation missing for timepoint {idx + 1}. "
                    "Provide a global 'labels' or specify 'labels' per timepoint."
                )

        ref_labels_path = resolved_labels[0]

        if output_path is not None:
            output_path = Path(output_path)

        if method_params is None:
            method_params = {}

        labels_map = check_dosimetry_inputs(
            method=method,
            isotope=isotope,
            labels_path=ref_labels_path,
            rois=rois,
            method_params=method_params,
        )

        for p in extracted_paths:
            if isinstance(p, Path) and not p.exists():
                raise FileNotFoundError(f"Activity image not found: {p}")
        for lbl in resolved_labels:
            if isinstance(lbl, Path) and not lbl.exists():
                raise FileNotFoundError(f"Labels image not found: {lbl}")
        for c in resolved_cts:
            if isinstance(c, Path) and not c.exists():
                raise FileNotFoundError(f"CT image not found: {c}")

        results = compute_roi_dose_mtp(
            images=extracted_paths,
            labels=resolved_labels,
            method=method,
            isotope=isotope,
            times_h=times_h,
            ct_image=resolved_cts if any(c is not None for c in resolved_cts) else None,
            rois=rois,
            labels_map=labels_map,
            labels_path=resolved_labels,
            ct_path=resolved_cts if any(c is not None for c in resolved_cts) else None,
            method_params=method_params,
            verbose=verbose,
        )

        if output_path:
            logger.info(f"Saving dose results to: {output_path}")
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, "w") as f:
                json.dump(results, f, indent=4)
        else:
            pretty_print_dose(results)

        return results

    if not image_path:
        raise ValueError("Activity image path ('image') must be provided.")
    if not method:
        raise ValueError("Dosimetry method ('method') must be provided.")
    if not isotope:
        raise ValueError("Radionuclide ('isotope' / 'radionuclide') must be provided.")

    image_path = Path(image_path)
    if labels_path is not None:
        labels_path = Path(labels_path)
    if ct_path is not None:
        ct_path = Path(ct_path)
    if output_path is not None:
        output_path = Path(output_path)

    # 1. Fast boundary/input checks before reading heavy images
    if method_params is None:
        method_params = {}

    method_key = method.lower()
    labels_map = {}
    if method_key != "ldm":
        if not labels_path:
            raise ValueError("Labels path is required for S-value dosimetry methods.")
        labels_map = check_dosimetry_inputs(
            method=method,
            isotope=isotope,
            labels_path=labels_path,
            rois=rois,
            method_params=method_params,
        )
    else:
        if labels_path is not None:
            labels_map = check_dosimetry_inputs(
                method=method,
                isotope=isotope,
                labels_path=labels_path,
                rois=rois,
                method_params=method_params,
            )
        else:
            check_dosimetry_inputs(
                method=method,
                isotope=isotope,
                labels_path=None,
                rois=None,
                method_params=None,
            )

    if resample_target == "like_ct" and not ct_path:
        raise ValueError(
            "Cannot resample like CT image because companion CT path was not provided."
        )

    if not image_path.exists():
        raise FileNotFoundError(f"Activity image not found: {image_path}")
    if labels_path is not None and not labels_path.exists():
        raise FileNotFoundError(f"Labels image not found: {labels_path}")

    logger.info(f"Reading activity image: {image_path}")
    image = sitk.ReadImage(str(image_path))

    labels = None
    if labels_path is not None:
        logger.info(f"Reading labels image: {labels_path}")
        labels = sitk.ReadImage(str(labels_path))

    ct_image = None
    if ct_path:
        if not ct_path.exists():
            raise FileNotFoundError(f"CT image not found: {ct_path}")
        logger.info(f"Reading companion CT image: {ct_path}")
        ct_image = sitk.ReadImage(str(ct_path))

    # Perform resampling if target is specified
    if resample_target:
        logger.info(f"Resampling images to align grids: target is {resample_target}")
        from radot.resampling import resample_image
        from radot.meta import get_metadata_case_insensitive

        modality = (
            get_metadata_case_insensitive(image, "modality", image_path) or "spect"
        )
        unit = get_metadata_case_insensitive(image, "unit", image_path)
        preserve_total = unit is not None and unit.lower() == "bq"

        if resample_target == "like_activity":
            if labels is not None:
                labels = resample_image(
                    labels, image_type="labels", reference_image=image
                )
            if ct_image is not None:
                ct_image = resample_image(
                    ct_image,
                    image_type="ct",
                    preserve_total=False,
                    gaussian=True,
                    reference_image=image,
                )
        elif resample_target == "like_labels":
            if labels is None:
                raise ValueError(
                    "Cannot resample like labels because labels image is not provided."
                )
            image = resample_image(
                image,
                image_type=modality.lower(),
                preserve_total=preserve_total,
                gaussian=True,
                reference_image=labels,
            )
            if ct_image is not None:
                ct_image = resample_image(
                    ct_image,
                    image_type="ct",
                    preserve_total=False,
                    gaussian=True,
                    reference_image=labels,
                )
        elif resample_target == "like_ct":
            if ct_image is None:
                raise ValueError(
                    "Cannot resample like CT because CT image is not provided."
                )
            image = resample_image(
                image,
                image_type=modality.lower(),
                preserve_total=preserve_total,
                gaussian=True,
                reference_image=ct_image,
            )
            if labels is not None:
                labels = resample_image(
                    labels, image_type="labels", reference_image=ct_image
                )
        else:
            raise ValueError(f"Unknown resample target: {resample_target}")

    results = compute_roi_dose(
        image=image,
        labels=labels,
        method=method,
        isotope=isotope,
        rois=rois,
        ct_image=ct_image,
        total_activity_mbq=total_activity_mbq,
        image_path=image_path,
        labels_path=labels_path,
        verbose=verbose,
        method_params=method_params,
        hu_threshold=hu_threshold,
        dose_rate=dose_rate,
        cascade=cascade,
    )

    if method_key == "ldm":
        from radot.meta import read_metadata, write_metadata

        if labels_path is not None:
            results_dict, dose_image = results
            if output_path:
                if output_path.name.endswith((".nii", ".nii.gz", ".nrrd", ".mhd")):
                    logger.info(f"Saving dose image to: {output_path}")
                    output_path.parent.mkdir(parents=True, exist_ok=True)
                    sitk.WriteImage(dose_image, str(output_path))
                    try:
                        orig_meta = read_metadata(image_path).get("metadata", {})
                    except Exception:
                        orig_meta = {}
                    write_metadata(
                        image_path=output_path,
                        modality="dose",
                        unit="gy/h" if dose_rate else "gy",
                        weight=orig_meta.get("weight"),
                        activity=orig_meta.get("activity"),
                        injection_date=orig_meta.get("injection_date"),
                        acquisition_date=orig_meta.get("acquisition_date"),
                    )
                else:
                    logger.info(f"Saving dose results to: {output_path}")
                    output_path.parent.mkdir(parents=True, exist_ok=True)
                    with open(output_path, "w") as f:
                        json.dump(results_dict, f, indent=4)
            if not output_path or output_path.name.endswith(
                (".nii", ".nii.gz", ".nrrd", ".mhd")
            ):
                pretty_print_dose(results_dict)
            return results
        else:
            if not output_path:
                output_path = Path("dose.nii.gz")
            logger.info(f"Saving dose image to: {output_path}")
            output_path.parent.mkdir(parents=True, exist_ok=True)
            sitk.WriteImage(results, str(output_path))
            try:
                orig_meta = read_metadata(image_path).get("metadata", {})
            except Exception:
                orig_meta = {}
            write_metadata(
                image_path=output_path,
                modality="dose",
                unit="gy/h" if dose_rate else "gy",
                weight=orig_meta.get("weight"),
                activity=orig_meta.get("activity"),
                injection_date=orig_meta.get("injection_date"),
                acquisition_date=orig_meta.get("acquisition_date"),
            )
            return results

    # Print results to stdout
    if output_path:
        logger.info(f"Saving dose results to: {output_path}")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(results, f, indent=4)
    else:
        # Pretty print results as a table
        pretty_print_dose(results)

    return results


def pretty_print_dose(results: dict) -> None:
    """Pretty print dosimetry results as a table."""
    import sys

    is_rate = any("mean_dose_rate_gyh" in data for data in results.values())
    has_tia = any("tia_mbq_h" in data for data in results.values())
    if is_rate:
        headers = ["ROI", "Mean Dose Rate (Gy/h)", "Mass (g)", "Volume (mL)"]
    elif has_tia:
        headers = [
            "ROI",
            "Mean Dose (Gy)",
            "Mass (g)",
            "Volume (mL)",
            "Teff (h)",
            "TIA (MBq.h)",
        ]
    else:
        headers = ["ROI", "Mean Dose (Gy)", "Mass (g)", "Volume (mL)", "Teff (h)"]

    rows = []
    for roi, data in results.items():
        if is_rate:
            mean_val = data.get("mean_dose_rate_gyh", 0.0)
            rows.append(
                [
                    roi,
                    f"{mean_val:.3f}",
                    f"{data['roi_mass_g']:.3f}",
                    f"{data['roi_volume_ml']:.3f}",
                ]
            )
        elif has_tia:
            rows.append(
                [
                    roi,
                    f"{data['mean_dose_gy']:.3f}",
                    f"{data['roi_mass_g']:.3f}",
                    f"{data['roi_volume_ml']:.3f}",
                    f"{data['roi_teff_h']:.1f}",
                    f"{data.get('tia_mbq_h', 0.0):.2f}",
                ]
            )
        else:
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
        sys.stdout.write("No dose computed.\n")
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

    sys.stdout.write("\n" + header_line + "\n")
    sys.stdout.write(separator + "\n")
    for row in rows:
        sys.stdout.write(format_row(row) + "\n")
    sys.stdout.write("\n")
