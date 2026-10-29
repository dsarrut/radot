"""Multi-time-point dosimetry: Jackson 2013 tri-exponential fit (3 time points)."""

import json
import logging
from typing import List, Optional, Dict, Any, Union
import numpy as np
import SimpleITK as sitk
from pathlib import Path

from radot.db import get_physical_half_life_h

logger = logging.getLogger("radot")


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
