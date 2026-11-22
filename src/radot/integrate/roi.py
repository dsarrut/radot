"""ROI-level time integration (fit first): ROI time-activity curves from activity images and labels,
then one fit per ROI (radot.integrate.models). docs/6e_time_integration.md."""

import json
import logging
import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

import numpy as np
import SimpleITK as sitk

from radot.integrate.models import check_model_choice, integrate_tac, parse_error_model

logger = logging.getLogger("radot")

INTEGRATE_FORMAT = 1
ACTIVITY_UNITS = ("bq", "bqml")

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
) -> dict:
    """ROI activity (MBq) at the acquisition time of each image, sorted by time from injection.

    labels: one label image (registered time points) or one per image. ROI names from the label
    sidecar (or labels_map). Activities are undone from an 'admin' decay correction
    (meta.activity_at_acquisition_factor): never decay-corrected for the fit.
    A ROI absent at some time point is an error when requested (rois), else skipped with a WARNING.
    Returns {"times_h", "isotope", "decay_correction", "image_names", "rois": {name: [MBq per tp]},
    "volumes_ml": {name: [mL per tp]}}.
    """
    from radot.dose.compute import load_labels_map, resolve_isotope
    from radot.meta import (
        activity_at_acquisition_factor,
        get_decay_correction,
        get_metadata_case_insensitive,
        get_time_from_injection_h,
    )
    from radot.resampling import resample_image
    from radot.stats import has_same_grid

    n = len(images)
    if n < 2:
        raise ValueError(f"Time integration needs at least 2 images (time points), got {n}.")
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
        if unit is None or unit.lower() not in ACTIVITY_UNITS:
            raise ValueError(
                f"Time point {i + 1} ({names[i]}): unit {unit!r}; time integration of activity images needs "
                f"bq or bqml (convert first: radot convert-units). Dose-rate images are not supported yet."
            )
        units.append(unit.lower())
    if len(set(units)) > 1:
        raise ValueError(f"The images have different units: {', '.join(f'{a}: {u}' for a, u in zip(names, units))}.")
    if len(set(times)) != n:
        raise ValueError(f"Two images have the same time from injection: {times} h.")
    if len(set(dcs)) > 1:
        logger.warning(
            "The images have different decay corrections ("
            + ", ".join(f"{a}: {d}" for a, d in zip(names, dcs))
            + "): each is brought back to the activity at its acquisition time."
        )

    order = np.argsort(times)
    curves: Dict[str, List[float]] = {}
    volumes: Dict[str, List[float]] = {}
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
            img = resample_image(img, image_type="spect", preserve_total=units[i] == "bq", gaussian=True,
                                 reference_image=lbl)
        sp = lbl.GetSpacing()
        vox_ml = sp[0] * sp[1] * sp[2] / 1000.0
        arr = sitk.GetArrayFromImage(img).astype(np.float64)
        mbq = arr * vox_ml / 1e6 if units[i] == "bqml" else arr / 1e6
        mbq = mbq * activity_at_acquisition_factor(dcs[i], isotope, times[i])
        lbl_arr = sitk.GetArrayFromImage(lbl)
        per_tp.append((lbl_arr, mbq, vox_ml, {v.lower(): k for k, v in lmap.items()}))
        for v in lmap.values():
            if v not in all_names:
                all_names.append(v)

    targets = rois if rois is not None else all_names
    for roi in targets:
        acts, vols = [], []
        for lbl_arr, mbq, vox_ml, name_to_id in per_tp:
            lid = name_to_id.get(roi.lower())
            mask = (lbl_arr == lid) if lid is not None else None
            nvox = int(mask.sum()) if mask is not None else 0
            acts.append(float(mbq[mask].sum()) if nvox else 0.0)
            vols.append(nvox * vox_ml)
        missing = [str(k + 1) for k, v in enumerate(vols) if v <= 0]
        if missing:
            msg = f"ROI '{roi}' has no voxel at time point(s) {', '.join(missing)} (sorted by time)"
            if rois is not None:
                raise ValueError(msg + ": the fit needs it at every time point.")
            logger.warning(msg + ": skipped.")
            continue
        curves[roi] = acts
        volumes[roi] = vols
    return {
        "times_h": [times[i] for i in order],
        "isotope": isotope,
        "decay_correction": [dcs[i] for i in order],
        "image_names": [names[i] for i in order],
        "rois": curves,
        "volumes_ml": volumes,
    }


def integrate_rois(
    images: Sequence[ImageLike],
    labels: Union[ImageLike, Sequence[ImageLike]],
    model: str,
    error_model: str,
    candidates: Optional[Sequence[str]] = None,
    select: Optional[str] = None,
    isotope: Optional[str] = None,
    rois: Optional[Sequence[str]] = None,
    image_paths=None,
    labels_paths=None,
    labels_map: Optional[Dict[int, str]] = None,
    resample: Optional[str] = None,
) -> dict:
    """TIA per ROI (fit first, ROI level). Returns the integrate document (see integrate_document):
    {"radot_integrate", ..., "rois": {name: {"tia_mbq_h", "tia_mbq_s", "tia_se_mbq_h", "model", "fits",
    "times_h", "activities_mbq", ...}}}. A ROI whose fit fails is an error when requested (rois), else
    skipped with a WARNING."""
    from radot.db import get_physical_half_life_h

    em = parse_error_model(error_model)
    tac = roi_time_activity_curves(images, labels, isotope=isotope, rois=rois, image_paths=image_paths,
                                   labels_paths=labels_paths, labels_map=labels_map, resample=resample)
    lam_p = math.log(2.0) / get_physical_half_life_h(tac["isotope"])
    # the model choice depends only on the number of time points: one error, not one per ROI
    check_model_choice(model, em, len(tac["times_h"]), candidates, select)
    out = {}
    at_bound = []
    for roi, acts in tac["rois"].items():
        try:
            r = integrate_tac(tac["times_h"], acts, lam_p, model, em, candidates=candidates, select=select)
        except ValueError as e:
            if rois:
                raise ValueError(f"ROI '{roi}': {e}") from e
            logger.warning(f"ROI '{roi}' skipped: {e}")
            continue
        for name, f in r["fits"].items():
            for issue in f["issues"]:
                if name == r["model"] or r["model"] == "average":
                    logger.info(f"ROI '{roi}', model {name}: {issue}")
            if f["at_bound"] and name == r["model"]:
                at_bound.append(roi)
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
    settings = {
        "quantity": "activity",
        "radionuclide": tac["isotope"],
        "lambda_phys_h": lam_p,
        "model": model,
        "candidates": list(candidates) if candidates else None,
        "select": select,
        "error_model": em.name,
        "times_h": tac["times_h"],
        "images": tac["image_names"],
        "decay_correction": tac["decay_correction"],
    }
    return integrate_document(out, settings)


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
