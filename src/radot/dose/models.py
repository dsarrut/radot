"""Dose formulas of the single-time-point models (LDM, Hanscheid, Madsen, S-value)."""

from typing import Optional
import numpy as np

from radot.db import get_physical_half_life_h


def dose_ldm(
    image_arr: np.ndarray,
    voxel_volume_ml: float,
    time_from_injection_h: float,
    isotope: str,
    ct_arr: Optional[np.ndarray] = None,
    hu_threshold: Optional[float] = -800.0,
    dose_rate: bool = False,
    cascade: bool = False,
    tia: bool = False,
) -> np.ndarray:
    """
    Voxel-level physical Local Deposition Model (LDM).
    - image_arr is physical activity in Bq per voxel,
      or time-integrated activity (TIA) in Bq.s per voxel if tia is True
      (time_from_injection_h is then ignored)
    - returns 3D float numpy array of absorbed dose (in Gy) or dose rate (in Gy/h) per voxel
    """
    from radot.db import get_cascade_energy_constant_mev

    if tia and dose_rate:
        raise ValueError("Dose rate cannot be computed from a time-integrated activity image.")

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
    elif tia:
        # Absorbed Dose (Gy) per voxel: TIA (Bq.s) * Delta_local (MeV) * conversions / mass (g)
        dose_arr = (image_arr * delta_local_mev * mev_to_joules_per_g) / voxel_mass_g
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


def dose_svalue(
    tia_bqs_a: np.ndarray,
    roi_a: np.ndarray,
    s_value: float,
    mass_scaling: float,
    source_roi_a: Optional[np.ndarray] = None,
) -> float:
    """
    Absorbed dose from a time-integrated activity (TIA) image and an S-value.
    - tia_bqs_a must be in Bq.s per voxel
    - output is in Gray
    - if source_roi_a is provided, TIA is summed in source_roi_a (source organ),
      computing cross-dose to target organ roi_a using s_value S(target <- source)
      and target mass_scaling.
    """
    source_mask = source_roi_a if source_roi_a is not None else roi_a
    tia_mbq_s = np.sum(tia_bqs_a[source_mask == 1]) / 1e6

    # S is in (mGy/MBq/s), so we get dose in mGy. Scale and convert to Gy.
    dose = tia_mbq_s * s_value / mass_scaling / 1000.0

    return float(dose)
