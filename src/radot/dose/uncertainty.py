"""ROI dose uncertainty from a Monte Carlo dose map and its voxel uncertainty."""

from typing import Optional, Union
import numpy as np
import SimpleITK as sitk


# Approximate CSDA range in water at the mean beta energy (in mm).
# Most electrons stop well before the maximum beta range: the mean range sets the
# size of the voxel cluster correlated by one electron history.
ISOTOPE_MEAN_BETA_RANGE_MM = {
    "lu177": 0.23,  # mean beta energy 134 keV
    "y90": 4.0,  # 934 keV
    "i131": 0.38,  # 182 keV
    "tb161": 0.29,  # 154 keV
    "re188": 3.2,  # 763 keV
    "sm153": 0.54,  # 225 keV
    "p32": 2.8,  # 695 keV
    "cu67": 0.25,  # 141 keV
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


def roi_mean_chord_mm(
    mask: Union[np.ndarray, sitk.Image],
    spacing_mm: Optional[tuple[float, float, float]] = None,
) -> float:
    """
    Mean chord length of an ROI for isotropic random lines (Cauchy formula): 4 V / S.

    The surface S is the area of the voxel faces on the ROI border, times 2/3 (voxel faces
    overestimate the area of a smooth surface by about 3/2).

    Parameters:
        mask: ROI binary mask (voxels > 0 are included).
        spacing_mm: Voxel spacing in mm (x, y, z). Auto-extracted if mask is sitk.Image.

    Returns:
        float: Mean chord length in mm.
    """
    if spacing_mm is None and isinstance(mask, sitk.Image):
        spacing_mm = mask.GetSpacing()
    elif spacing_mm is None:
        spacing_mm = (1.0, 1.0, 1.0)
    m_arr = (
        sitk.GetArrayViewFromImage(mask)
        if isinstance(mask, sitk.Image)
        else np.asarray(mask)
    )
    m = (m_arr > 0).astype(np.int8)
    if not np.any(m):
        raise ValueError("ROI mask contains no voxels.")

    sx, sy, sz = (float(s) for s in spacing_mm)
    # Array axes are (z, y, x): faces normal to each axis
    face_areas = (sx * sy, sx * sz, sy * sz)
    surface = 0.0
    for axis, area in enumerate(face_areas):
        pad = [(1, 1) if a == axis else (0, 0) for a in range(3)]
        surface += np.abs(np.diff(np.pad(m, pad), axis=axis)).sum() * area
    surface *= 2.0 / 3.0
    volume = m.sum() * sx * sy * sz
    return float(4.0 * volume / surface)


def roi_dose_uncertainty_estimate(
    dose: Union[np.ndarray, sitk.Image],
    uncertainty: Union[np.ndarray, sitk.Image],
    mask: Union[np.ndarray, sitk.Image],
    spacing_mm: Optional[tuple[float, float, float]] = None,
    isotope: Optional[str] = None,
    correlation_length_mm: Optional[float] = None,
    photons: bool = False,
) -> float:
    """
    Estimate realistic ROI relative uncertainty accounting for the voxels correlated by one history.

    The voxel uncertainties of a Monte Carlo dose map are history-by-history, but one history
    deposits in several voxels, so the voxels are correlated. If one history deposits in about
    nu_cluster voxels of the ROI:
        unc_est = min(unc_max, unc_min * sqrt(nu_cluster))

    - electrons (default): a cluster of the mean beta range,
        nu_cluster = prod(1 + range / spacing_d)
    - photons=True (e.g. track-length estimator dose map): a photon deposits along its path
      across the ROI, nu_cluster = 1 + mean_chord / mean_spacing (see roi_mean_chord_mm).

    For a map mixing electrons and photons, estimate each part separately when possible
    (e.g. the e/ and gamma_tle/ outputs of an OpenGATE VRT simulation) and add the variances:
    with photons, the electron estimate of the merged map underestimates the uncertainty.

    Parameters:
        dose: 3D dose or dose rate map (SimpleITK Image or NumPy array).
        uncertainty: 3D relative uncertainty map.
        mask: ROI binary mask.
        spacing_mm: Voxel spacing in mm (x, y, z). Auto-extracted if dose is sitk.Image.
        isotope: Radionuclide name (e.g. 'Lu177', 'Y90') to retrieve the mean beta range.
        correlation_length_mm: Direct override for the correlation length in mm.
        photons: Use the ROI mean chord as correlation length (photon dose map).

    Returns:
        float: Estimated relative uncertainty as a fraction.
    """
    if spacing_mm is None and isinstance(dose, sitk.Image):
        spacing_mm = dose.GetSpacing()
    elif spacing_mm is None:
        spacing_mm = (1.0, 1.0, 1.0)

    unc_min, unc_max = roi_dose_uncertainty_bounds(dose, uncertainty, mask)

    if photons:
        length = correlation_length_mm
        if length is None:
            length = roi_mean_chord_mm(mask, spacing_mm)
        nu_cluster = 1.0 + float(length) / float(np.mean(spacing_mm))
    else:
        if correlation_length_mm is None:
            iso_key = (isotope or "").lower().replace("-", "")
            correlation_length_mm = ISOTOPE_MEAN_BETA_RANGE_MM.get(iso_key, 2.0)
        # Effective correlation cluster size across spatial dimensions
        nu_cluster = 1.0
        for s in spacing_mm:
            if s > 0:
                nu_cluster *= 1.0 + float(correlation_length_mm) / float(s)

    unc_est = float(min(unc_max, unc_min * np.sqrt(nu_cluster)))
    return unc_est
