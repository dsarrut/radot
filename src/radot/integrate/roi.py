"""ROI-level time integration: ROI curves from images and labels, then one fit per ROI
(radot.integrate.models). docs/6e_time_integration.md.

Two quantities, from the image units:
- activity (bq, bqml): fit first, ROI activity (MBq) -> TIA (MBq.h), then `radot dose` on the table;
- dose rate (gy/h, gy/s, e.g. `radot dose -m ldm --dose-rate` at each time point): fit after, ROI mean
  dose rate (Gy/h, mass-weighted with the CT of each time point) -> mean absorbed dose (Gy)."""

import json
import logging
import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import numpy as np
import SimpleITK as sitk

from radot.integrate.models import parse_error_model
from radot.integrate.simple import TEFF_MODELS, check_integration_options, integrate_curve

logger = logging.getLogger("radot")

INTEGRATE_FORMAT = 1
ACTIVITY_UNITS = ("bq", "bqml")
DOSE_RATE_UNITS = {"gy/h": 1.0, "gy/s": 3600.0}  # unit -> factor to Gy/h
QUANTITIES = ("activity", "dose_rate")


def quantity_of_unit(unit: Optional[str]) -> Optional[str]:
    """'activity' (bq, bqml), 'dose_rate' (gy/h, gy/s), else None."""
    u = (unit or "").strip().lower()
    if u in ACTIVITY_UNITS:
        return "activity"
    if u in DOSE_RATE_UNITS:
        return "dose_rate"
    return None

ImageLike = Union[str, Path, sitk.Image]


def _read(item, path=None):
    if isinstance(item, (str, Path)):
        p = Path(item)
        return sitk.ReadImage(str(p)), p
    return item, (Path(path) if path else None)


def roi_time_activity_curves(
    images: Sequence[ImageLike],
    labels: Union[ImageLike, Sequence[ImageLike]],
    isotope: Optional[str] = None,
    rois: Optional[Sequence[str]] = None,
    image_paths: Optional[Sequence[Optional[Union[str, Path]]]] = None,
    labels_paths: Optional[Sequence[Optional[Union[str, Path]]]] = None,
    labels_map: Optional[Dict[int, str]] = None,
    resample: Optional[str] = None,
    ct: Optional[Union[ImageLike, Sequence[ImageLike]]] = None,
) -> dict:
    """ROI curve at the acquisition time of each image, sorted by time from injection: the activity
    (MBq) of activity images (bq, bqml), or the mean dose rate (Gy/h) of dose-rate images (gy/h, gy/s).

    labels: one label image (registered time points) or one per image. ROI names from the label
    sidecar (or labels_map). Activities are undone from an 'admin' decay correction
    (meta.activity_at_acquisition_factor): never decay-corrected for the fit.
    A ROI absent at some time point is an error when requested (rois), else skipped with a WARNING.
    ct (dose rate only, one image or one per image): the ROI mean dose rate is mass-weighted with the CT
    density of that time point (stats.density_map_g_ml; without a CT, water density with a WARNING).
    Returns {"quantity", "times_h", "isotope", "decay_correction", "image_names", "sex",
    "rois": {name: [MBq or Gy/h per tp]}, "volumes_ml": {name: [mL per tp]}, "masses_g" (dose rate)}.
    """
    from radot.dose.compute import load_labels_map, resolve_isotope
    from radot.meta import (
        activity_at_acquisition_factor,
        get_decay_correction,
        get_metadata_case_insensitive,
        get_time_from_injection_h,
    )
    from radot.resampling import resample_image
    from radot.stats import density_map_g_ml, has_same_grid

    n = len(images)
    if n < 1:
        raise ValueError("Time integration needs at least one image.")
    image_paths = list(image_paths) if image_paths else [None] * n
    imgs, paths = zip(*[_read(im, image_paths[i] if i < len(image_paths) else None) for i, im in enumerate(images)])
    if isinstance(labels, (list, tuple)):
        if len(labels) not in (1, n):
            raise ValueError(f"{len(labels)} label images for {n} images: give one, or one per image.")
        labels = list(labels) * (n if len(labels) == 1 else 1)
    else:
        labels = [labels] * n
    labels_paths = list(labels_paths) if labels_paths else [None] * n
    if len(labels_paths) == 1:
        labels_paths = labels_paths * n
    lbls, lpaths = zip(*[_read(lb, labels_paths[i]) for i, lb in enumerate(labels)])
    if ct is None:
        cts = [None] * n
    else:
        cts = list(ct) if isinstance(ct, (list, tuple)) else [ct]
        if len(cts) not in (1, n):
            raise ValueError(f"{len(cts)} CT images for {n} images: give one, or one per image.")
        cts = [_read(c)[0] for c in cts] * (n if len(cts) == 1 else 1)
    if resample not in (None, "like_labels"):
        raise ValueError(f"resample '{resample}': only 'like_labels' (each image to the grid of its labels).")

    isotope = resolve_isotope(isotope, imgs, paths)
    rois = list(rois) if rois else None
    names = [p.name if p else f"image #{i + 1}" for i, p in enumerate(paths)]

    # time, unit, decay correction of each image
    times, units, dcs = [], [], []
    for i, (img, p) in enumerate(zip(imgs, paths)):
        try:
            times.append(float(get_time_from_injection_h(img, p)))
            dcs.append(get_decay_correction(img, p))
        except ValueError as e:
            raise ValueError(f"Time point {i + 1} ({names[i]}): {e}") from e
        unit = get_metadata_case_insensitive(img, "unit", p)
        if quantity_of_unit(unit) is None:
            raise ValueError(
                f"Time point {i + 1} ({names[i]}): unit {unit!r}; time integration needs activity images "
                f"(bq, bqml; convert first: radot convert-units) or dose-rate images (gy/h, gy/s; radot dose "
                f"-m ldm --dose-rate)."
            )
        units.append(unit.lower())
    quantities = {quantity_of_unit(u) for u in units}
    if len(quantities) > 1:
        raise ValueError(
            "Activity and dose-rate images cannot be mixed: "
            + ", ".join(f"{a}: {u}" for a, u in zip(names, units)) + "."
        )
    quantity = quantities.pop()
    if quantity == "activity" and any(c is not None for c in cts):
        raise ValueError(
            "A CT is used only with dose-rate images (mass-weighted mean dose rate); for activity images the "
            "masses come with `radot dose --ct` on the table."
        )
    if quantity == "dose_rate" and all(c is None for c in cts):
        logger.warning(
            "No CT: the ROI mean dose rates are weighted with water density 1 g/mL (give the CT of each "
            "time point, --ct, for the CT density)."
        )
    if len(set(times)) != n:
        raise ValueError(f"Two images have the same time from injection: {times} h.")
    # patient sex (for the OpenDose phantom of a later `radot dose` on the table): the image sidecars
    sexes = {str(v).strip().lower() for v in (get_metadata_case_insensitive(im, "sex", p) for im, p in zip(imgs, paths))
             if v is not None and str(v).strip()}
    if len(sexes) > 1:
        raise ValueError(f"The image sidecars give different patient sexes: {', '.join(sorted(sexes))}.")
    sex = sexes.pop() if sexes else None
    if len(set(dcs)) > 1:
        logger.warning(
            "The images have different decay corrections ("
            + ", ".join(f"{a}: {d}" for a, d in zip(names, dcs))
            + "): each is brought back to the activity at its acquisition time."
        )

    order = np.argsort(times)
    curves: Dict[str, List[float]] = {}
    volumes: Dict[str, List[float]] = {}
    masses: Dict[str, List[float]] = {}
    all_names: List[str] = []
    per_tp = []
    for i in order:
        img, p, lbl, lp = imgs[i], paths[i], lbls[i], lpaths[i]
        if lp is not None:
            lmap = load_labels_map(lp)
        elif labels_map:
            lmap = dict(labels_map)
        else:
            raise ValueError(f"Time point {i + 1}: no ROI names (no label image path with a sidecar, no labels_map).")
        if not has_same_grid(img, lbl):
            if resample != "like_labels":
                raise ValueError(
                    f"Time point {i + 1} ({names[i]}): the image grid does not match its labels grid. Resample "
                    f"first, or give --resample like_labels."
                )
            logger.info(f"Time point {i + 1}: image resampled to the labels grid.")
            img = resample_image(img, image_type="spect" if quantity == "activity" else "dose",
                                 preserve_total=units[i] == "bq", gaussian=True, reference_image=lbl)
        sp = lbl.GetSpacing()
        vox_ml = sp[0] * sp[1] * sp[2] / 1000.0
        arr = sitk.GetArrayFromImage(img).astype(np.float64)
        factor = activity_at_acquisition_factor(dcs[i], isotope, times[i])
        density = None
        if quantity == "activity":
            values = (arr * vox_ml / 1e6 if units[i] == "bqml" else arr / 1e6) * factor
        else:
            values = arr * DOSE_RATE_UNITS[units[i]] * factor
            ct_img = cts[i]
            if ct_img is not None and not has_same_grid(ct_img, lbl):
                if resample != "like_labels":
                    raise ValueError(
                        f"Time point {i + 1}: the CT grid does not match its labels grid. Resample first, or "
                        f"give --resample like_labels."
                    )
                ct_img = resample_image(ct_img, image_type="ct", gaussian=True, reference_image=lbl)
            ct_arr = sitk.GetArrayFromImage(ct_img) if ct_img is not None else None
            density = density_map_g_ml(ct_arr, arr.shape, log=False)  # water: logged once above
        lbl_arr = sitk.GetArrayFromImage(lbl)
        per_tp.append((lbl_arr, values, vox_ml, density, {v.lower(): k for k, v in lmap.items()}))
        for v in lmap.values():
            if v not in all_names:
                all_names.append(v)

    targets = rois if rois is not None else all_names
    for roi in targets:
        acts, vols, ms = [], [], []
        for lbl_arr, values, vox_ml, density, name_to_id in per_tp:
            lid = name_to_id.get(roi.lower())
            mask = (lbl_arr == lid) if lid is not None else None
            nvox = int(mask.sum()) if mask is not None else 0
            vols.append(nvox * vox_ml)
            if not nvox:
                acts.append(0.0)
                ms.append(0.0)
            elif density is None:
                acts.append(float(values[mask].sum()))
            else:  # mass-weighted mean dose rate, as the ROI mean dose of `radot dose`
                w = density[mask]
                acts.append(float(np.sum(values[mask] * w) / np.sum(w)))
                ms.append(float(np.sum(w) * vox_ml))
        missing = [str(k + 1) for k, v in enumerate(vols) if v <= 0]
        if missing:
            msg = f"ROI '{roi}' has no voxel at time point(s) {', '.join(missing)} (sorted by time)"
            if rois is not None:
                raise ValueError(msg + ": the fit needs it at every time point.")
            logger.warning(msg + ": skipped.")
            continue
        curves[roi] = acts
        volumes[roi] = vols
        if quantity == "dose_rate":
            masses[roi] = ms
    return {
        "quantity": quantity,
        "times_h": [times[i] for i in order],
        "isotope": isotope,
        "decay_correction": [dcs[i] for i in order],
        "image_names": [names[i] for i in order],
        "sex": sex,
        "rois": curves,
        "volumes_ml": volumes,
        "masses_g": masses if quantity == "dose_rate" else None,
    }


def integrate_rois(
    images: Sequence[ImageLike],
    labels: Union[ImageLike, Sequence[ImageLike]],
    model: str,
    error_model: Optional[str] = None,
    candidates: Optional[Sequence[str]] = None,
    select: Optional[str] = None,
    isotope: Optional[str] = None,
    rois: Optional[Sequence[str]] = None,
    image_paths=None,
    labels_paths=None,
    labels_map: Optional[Dict[int, str]] = None,
    resample: Optional[str] = None,
    head: Optional[str] = None,
    tail: Optional[str] = None,
    tail_points: Optional[int] = None,
    teff: Optional[Dict[str, float]] = None,
    ct: Optional[Union[ImageLike, Sequence[ImageLike]]] = None,
) -> dict:
    """Integral per ROI, ROI level. Returns the integrate document (see integrate_document).
    Activity images (fit first): {"radot_integrate", "quantity": "activity", ..., "rois": {name:
    {"tia_mbq_h", "tia_mbq_s", "tia_se_mbq_h", "model", "fits", "times_h", "activities_mbq", ...}}}.
    Dose-rate images (fit after, ct: the CT of each time point for the mass-weighted mean): "quantity":
    "dose_rate", "rois": {name: {"mean_dose_gy", "mean_dose_se_gy", "model", "fits", "times_h",
    "dose_rates_gy_h", "masses_g", ...}}; the fits give "dose_gy" / "dose_se_gy" instead of the TIA keys.

    model: a fitted model (f1 ... f6a, auto: error_model, candidates, select), 'trapezoid' (head, tail,
    tail_points) or a single time point model (physical, hanscheid2018, madsen2018: one image; teff
    {roi: hours}, required for madsen2018, else hanscheid2018 uses 2t). The options are checked once
    (check_integration_options). A ROI whose integration fails is an error when requested (rois), else
    skipped with a WARNING."""
    from radot.db import get_physical_half_life_h
    from radot.dose.models import log_teff_sources, resolve_teff_h

    key = (model or "").strip().lower()
    # the options depend only on the model and the number of time points: one error, not one per ROI
    check_integration_options(model, len(images), error_model, candidates, select, head, tail, tail_points,
                              teff_given=bool(teff))
    tac = roi_time_activity_curves(images, labels, isotope=isotope, rois=rois, image_paths=image_paths,
                                   labels_paths=labels_paths, labels_map=labels_map, resample=resample, ct=ct)
    dose_rate = tac["quantity"] == "dose_rate"
    lam_p = math.log(2.0) / get_physical_half_life_h(tac["isotope"])
    teff_params = {r: {"teff": float(v)} for r, v in (teff or {}).items()}
    out = {}
    at_bound = []
    for roi, acts in tac["rois"].items():
        teff_h = teff_source = None
        if key in TEFF_MODELS:  # a missing Teff (madsen2018) is an error, not a skipped ROI
            teff_h, teff_source = resolve_teff_h(key, roi, teff_params, tac["times_h"][0], tac["isotope"])
        try:
            r = integrate_curve(tac["times_h"], acts, lam_p, model, error_model=error_model,
                                candidates=candidates, select=select, head=head, tail=tail,
                                tail_points=tail_points, teff_h=teff_h, teff_source=teff_source)
        except ValueError as e:
            if rois:
                raise ValueError(f"ROI '{roi}': {e}") from e
            logger.warning(f"ROI '{roi}' skipped: {e}")
            continue
        for name, f in r["fits"].items():
            for issue in f["issues"]:
                if name == r["model"] or r["model"] == "average":
                    logger.info(f"ROI '{roi}', model {name}: {issue}")
            if f.get("at_bound") and name == r["model"]:
                at_bound.append(roi)
        if dose_rate:
            out[roi] = {
                "mean_dose_gy": r["tia_mbq_h"],  # integral of the dose rate (Gy/h) over h
                "mean_dose_se_gy": r["tia_se_mbq_h"],
                "model": r["model"],
                "times_h": tac["times_h"],
                "dose_rates_gy_h": acts,
                "volumes_ml": tac["volumes_ml"][roi],
                "masses_g": tac["masses_g"][roi],
                "fits": {k: dose_fit_keys(f) for k, f in r["fits"].items()},
                "excluded": r["excluded"],
            }
            continue
        out[roi] = {
            "tia_mbq_h": r["tia_mbq_h"],
            "tia_mbq_s": r["tia_mbq_h"] * 3600.0,
            "tia_se_mbq_h": r["tia_se_mbq_h"],
            "model": r["model"],
            "times_h": tac["times_h"],
            "activities_mbq": acts,
            "volumes_ml": tac["volumes_ml"][roi],
            "fits": r["fits"],
            "excluded": r["excluded"],
        }
    if at_bound:
        logger.warning(
            f"Biological clearance at the bound 0 (pure physical decay) for: {', '.join(at_bound)}."
        )
    if key in TEFF_MODELS:
        log_teff_sources(key, out.keys(), teff_params, tac["times_h"][0], tac["isotope"])
    settings = {
        "quantity": tac["quantity"],
        "radionuclide": tac["isotope"],
        "lambda_phys_h": lam_p,
        "model": model,
        "candidates": list(candidates) if candidates else None,
        "select": select,
        "error_model": parse_error_model(error_model).name if error_model else None,
        "head": head,
        "tail": tail,
        "tail_points": tail_points,
        "teff_h": dict(teff) if teff else None,
        "times_h": tac["times_h"],
        "images": tac["image_names"],
        "decay_correction": tac["decay_correction"],
        "sex": tac["sex"],
    }
    return integrate_document(out, settings)


# integral keys of a fit (activity, MBq.h) -> the same for a dose-rate curve (Gy)
DOSE_FIT_KEYS = {"tia_mbq_h": "dose_gy", "tia_se_mbq_h": "dose_se_gy", "parts_mbq_h": "parts_gy"}


def dose_fit_keys(fit: dict) -> dict:
    """A fit dict of a dose-rate curve: the integral keys named as a dose (DOSE_FIT_KEYS)."""
    return {DOSE_FIT_KEYS.get(k, k): v for k, v in fit.items()}


def integrate_document(rois: dict, settings: dict) -> dict:
    from radot import __version__

    return {"radot_integrate": INTEGRATE_FORMAT, "radot_version": __version__, **settings, "rois": rois}


def write_integrate(doc: dict, path: Union[str, Path]) -> None:
    path = Path(path)
    logger.info(f"Saving the time integration to: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(_finite(doc), f, indent=4, allow_nan=False)


def _finite(o):
    """JSON values: numpy scalars as Python numbers, a non-finite number (e.g. the CV of a parameter
    at 0) as null."""
    if isinstance(o, dict):
        return {k: _finite(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_finite(v) for v in o]
    if isinstance(o, (np.floating, np.integer)):
        o = o.item()
    if isinstance(o, float) and not math.isfinite(o):
        return None
    return o


def is_integrate_file(path) -> bool:
    """True for a radot integrate file (a JSON object with the 'radot_integrate' key)."""
    p = Path(path)
    if p.suffix.lower() != ".json" or not p.is_file():
        return False
    try:
        with open(p) as f:
            doc = json.load(f)
    except (OSError, ValueError):
        return False
    return isinstance(doc, dict) and "radot_integrate" in doc


def read_integrate(path: Union[str, Path]) -> dict:
    """A radot integrate file (the whole document). Not one, or an unknown version, is an error."""
    path = Path(path)
    with open(path) as f:
        doc = json.load(f)
    if not isinstance(doc, dict) or "radot_integrate" not in doc or not isinstance(doc.get("rois"), dict):
        raise ValueError(f"{path.name} is not a radot integrate file (no 'radot_integrate' and 'rois' keys).")
    if doc["radot_integrate"] != INTEGRATE_FORMAT:
        raise ValueError(
            f"{path.name}: radot integrate format {doc['radot_integrate']!r}, this radot reads format "
            f"{INTEGRATE_FORMAT}."
        )
    return doc
