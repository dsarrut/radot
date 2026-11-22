"""Dose from a ROI TIA table (`radot integrate` output, fit first at ROI level), with OpenDose
self S-values scaled to the ROI mass (same formula as `svalue` on a TIA image):

    D = TIA (MBq.s) * S (mGy/MBq.s) / (m_ROI / m_phantom) / 1000

The ROI mass comes from a label image (and a CT: stats.ct_density_g_ml, voxel-wise; without a CT, water
density 1 g/mL, logged). docs/6e_time_integration.md §6."""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np
import SimpleITK as sitk

logger = logging.getLogger("radot")

TIA_TABLE_METHODS = ("svalue",)


def check_tia_table_options(
    method: str,
    cascade: Optional[bool] = None,
    hu_threshold: Optional[float] = None,
    dose_rate: Optional[bool] = None,
    total_activity_mbq: Optional[float] = None,
    dose_rate_map: bool = False,
    resample: Optional[str] = None,
    method_params: Optional[Dict[str, Any]] = None,
) -> None:
    """The method must be one of TIA_TABLE_METHODS; an option not used on a TIA table is an error."""
    if (method or "").lower() not in TIA_TABLE_METHODS:
        raise ValueError(
            f"A ROI TIA table (radot integrate) takes method {' or '.join(TIA_TABLE_METHODS)} "
            f"(got {method!r}): the time integration is already done."
        )
    unused = []
    if cascade:
        unused.append("cascade")
    if hu_threshold is not None:
        unused.append("hu_threshold")
    if dose_rate:
        unused.append("dose_rate")
    if total_activity_mbq is not None:
        unused.append("total_activity_mbq")
    if dose_rate_map:
        unused.append("dose_rate_map")
    if resample not in (None, "like_labels"):
        unused.append(f"resample '{resample}' (like_labels only: the CT to the labels grid)")
    teff = [r for r, v in (method_params or {}).items() if isinstance(v, dict) and "teff" in v]
    if teff:
        unused.append(f"teff for {', '.join(teff)}")
    if unused:
        raise ValueError("Not used with a ROI TIA table: " + "; ".join(unused) + ". Remove it.")


def compute_roi_dose_tia_table(
    table: Union[str, Path, dict],
    labels: Union[str, Path, sitk.Image],
    method: str = "svalue",
    ct_image: Optional[Union[str, Path, sitk.Image]] = None,
    isotope: Optional[str] = None,
    rois: Optional[List[str]] = None,
    method_params: Optional[Dict[str, Any]] = None,
    labels_path: Optional[Union[str, Path]] = None,
    labels_map: Optional[Dict[int, str]] = None,
    resample: Optional[str] = None,
    verbose: bool = False,
) -> dict:
    """Mean absorbed dose per ROI of a TIA table (radot.integrate). Returns the ROI results table.

    - radionuclide: the table; an `isotope` that contradicts it is an error;
    - phantom: method_params['phantom'], else the patient sex recorded in the table, else an error;
    - ROIs: all the ROIs of the table, or `rois` (a requested ROI missing from the table or the labels,
      or without OpenDose data, is an error; without `rois`, skipped with a WARNING);
    - mean_dose_fit_se_gy: the part of the dose uncertainty from the TIA fit only (S-value and mass not
      included).
    """
    from radot.db import (
        OrganNotFoundError,
        get_reference_mass,
        get_s_value,
        normalize_isotope,
        normalize_phantom,
        resolve_organ_name,
    )
    from radot.dose.compute import load_labels_map
    from radot.integrate import read_integrate
    from radot.resampling import resample_image
    from radot.stats import ct_density_g_ml, has_same_grid

    check_tia_table_options(method, resample=resample, method_params=method_params)
    if isinstance(table, (str, Path)):
        table_name = Path(table).name
        doc = read_integrate(table)
    else:
        table_name = "TIA table"
        doc = table
    if doc.get("quantity", "activity") != "activity":
        raise ValueError(f"{table_name}: quantity '{doc.get('quantity')}', a TIA table is 'activity'.")

    # radionuclide: the table (an option that contradicts it is an error)
    table_iso = doc.get("radionuclide")
    if not table_iso:
        raise ValueError(f"{table_name} has no radionuclide.")
    if isotope and normalize_isotope(isotope) != normalize_isotope(table_iso):
        raise ValueError(
            f"Radionuclide mismatch: the option says '{isotope}', {table_name} says '{table_iso}'. Drop the option."
        )
    isotope = table_iso

    # phantom: option/config, else the patient sex of the table
    phantom = (method_params or {}).get("phantom")
    if phantom is not None:
        normalize_phantom(phantom)
    else:
        sex = doc.get("sex")
        if sex not in ("male", "female"):
            raise ValueError(
                f"OpenDose phantom not provided: give --phantom male|female ({table_name} has no patient sex: "
                f"set it in the image sidecars, radot meta IMAGE --sex, before radot integrate)."
            )
        phantom = sex
        logger.info(f"OpenDose phantom: {phantom} (patient sex recorded in {table_name})")

    # labels (and CT) for the ROI masses
    if isinstance(labels, (str, Path)):
        labels_path = Path(labels)
        labels = sitk.ReadImage(str(labels_path))
    if labels_path is not None:
        labels_map = load_labels_map(labels_path)
    if not labels_map:
        raise ValueError("No ROI names for the labels (no label image path with a sidecar, no labels_map).")
    lbl_arr = sitk.GetArrayFromImage(labels)
    sp = labels.GetSpacing()
    vox_ml = sp[0] * sp[1] * sp[2] / 1000.0
    if ct_image is not None:
        if isinstance(ct_image, (str, Path)):
            ct_image = sitk.ReadImage(str(ct_image))
        if not has_same_grid(ct_image, labels):
            if resample != "like_labels":
                raise ValueError(
                    "The CT grid does not match the labels grid. Resample first, or give --resample like_labels."
                )
            logger.info("CT resampled to the labels grid.")
            ct_image = resample_image(ct_image, image_type="ct", preserve_total=False, gaussian=True,
                                      reference_image=labels)
        density = ct_density_g_ml(sitk.GetArrayFromImage(ct_image))
        mass_source = "CT"
    else:
        density = None
        mass_source = "volume, water density 1 g/mL"
        logger.info("No CT: ROI masses from the label volume with water density 1 g/mL.")
    name_to_id = {v.lower(): k for k, v in labels_map.items()}

    table_rois = doc["rois"]
    rois = list(rois) if rois else None
    if rois is not None:
        absent = [r for r in rois if r not in table_rois]
        if absent:
            raise ValueError(f"Requested ROI(s) not in {table_name}: {', '.join(absent)}.")
    targets = rois if rois is not None else list(table_rois)

    results = {}
    for roi in targets:
        entry = table_rois[roi]
        lid = name_to_id.get(roi.lower())
        mask = (lbl_arr == lid) if lid is not None else None
        nvox = int(mask.sum()) if mask is not None else 0
        try:
            if nvox == 0:
                raise LookupError(f"ROI '{roi}' has no voxel in the labels (needed for its mass)")
            organ = resolve_organ_name(roi, phantom)
            s_value = get_s_value(phantom, isotope, organ, organ)
            ref_mass = get_reference_mass(phantom, organ)
        except LookupError as e:
            if rois is not None:
                raise ValueError(f"{e}.") from e
            logger.warning(f"Skipping dosimetry for ROI '{roi}': {e}")
            continue
        except OrganNotFoundError as e:
            if rois is not None:
                raise
            logger.warning(f"Skipping dosimetry for ROI '{roi}': {e}")
            continue
        volume = nvox * vox_ml
        mass = float(np.sum(density[mask]) * vox_ml) if density is not None else volume
        mass_scaling = mass / ref_mass
        gy_per_mbq_s = s_value / mass_scaling / 1000.0
        tia_mbq_s = float(entry["tia_mbq_s"])
        se_h = entry.get("tia_se_mbq_h")
        results[roi] = {
            "mean_dose_gy": tia_mbq_s * gy_per_mbq_s,
            "mean_dose_fit_se_gy": se_h * 3600.0 * gy_per_mbq_s if se_h is not None else None,
            "roi_mass_g": mass,
            "roi_volume_ml": volume,
            "tia_mbq_s": tia_mbq_s,
            "tia_mbq_h": tia_mbq_s / 3600.0,
            "tia_se_mbq_h": se_h,
            "tia_model": entry.get("model"),
        }
        msg = (
            f"ROI '{roi}': TIA {tia_mbq_s / 3600.0:.4g} MBq.h ({entry.get('model')}), OpenDose organ '{organ}' "
            f"({phantom}), S {s_value:.4e} mGy/MBq.s, mass {mass:.4g} g ({mass_source}), reference "
            f"{ref_mass:.4g} g, dose {results[roi]['mean_dose_gy']:.4g} Gy"
        )
        logger.info(msg) if verbose else logger.debug(msg)
    return results
