"""
Quality control of the inputs of a dosimetry (radot check): metadata, units, geometry, ROIs.

The checks only read the files and report findings; they never modify the files and never raise on a
suspicious value. Each check is registered with the inputs it needs, and only the checks whose inputs
are given run. Thresholds are read from data/check_thresholds.json (overridable).
"""

import copy
import importlib.resources
import json
import logging
import re
from dataclasses import dataclass, field
from datetime import datetime
from enum import IntEnum
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence

import numpy as np
import SimpleITK as sitk

from radot.meta import get_sidecar_path, read_metadata

logger = logging.getLogger("radot")

DATE_FORMAT = "%Y-%m-%d %H:%M:%S"
ACTIVITY_UNITS = ("bq", "bqml", "suv", "sul", "bqs", "bqsml")
DOSE_UNITS = ("gy", "gy/h", "gy/s")
# Order of the categories in the reports
CATEGORIES = ("meta", "units", "activity", "geom", "roi")


class Level(IntEnum):
    OK = 0
    INFO = 1
    WARN = 2
    FAIL = 3


@dataclass
class Finding:
    code: str
    level: Level
    message: str
    subject: Optional[str] = None
    value: Any = None
    expected: Any = None

    @property
    def category(self) -> str:
        return self.code.split(".")[0]

    def to_dict(self) -> dict:
        return {
            "code": self.code,
            "level": self.level.name,
            "subject": self.subject,
            "message": self.message,
            "value": self.value,
            "expected": self.expected,
        }


@dataclass
class Check:
    code: str
    inputs: tuple
    description: str
    func: Callable


# Registry of the checks, by code (e.g. "meta.dates_order")
CHECKS: Dict[str, Check] = {}


def register_check(code: str, inputs: Sequence[str], description: str):
    """Decorator: register a check function f(ctx) -> list of Finding."""

    def decorator(func):
        CHECKS[code] = Check(code, tuple(inputs), description, func)
        return func

    return decorator


def _deep_merge(base: dict, other: dict) -> dict:
    for k, v in other.items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            _deep_merge(base[k], v)
        else:
            base[k] = v
    return base


def load_thresholds(path: Optional[Path] = None) -> dict:
    """Default thresholds, merged with a user JSON file (same sections) if given."""
    default = importlib.resources.files("radot").joinpath("data", "check_thresholds.json")
    with default.open("r", encoding="utf-8") as f:
        thresholds = json.load(f)
    if path is not None:
        with open(path, "r", encoding="utf-8") as f:
            _deep_merge(thresholds, json.load(f))
    return thresholds


@dataclass
class CheckContext:
    """Inputs of the checks; metadata read once and shared."""

    activity: List[Path] = field(default_factory=list)
    ct: Optional[Path] = None
    labels: Optional[Path] = None
    dose: Optional[Path] = None
    method: Optional[str] = None
    thresholds: dict = field(default_factory=dict)
    _meta: Dict[Path, dict] = field(default_factory=dict)
    _arrays: Dict[Path, tuple] = field(default_factory=dict)
    _images: Dict[Path, sitk.Image] = field(default_factory=dict)
    _roi_stats: Optional[Dict[int, dict]] = None

    def available(self) -> set:
        names = {"activity": bool(self.activity), "ct": self.ct, "labels": self.labels, "dose": self.dose}
        return {k for k, v in names.items() if v}

    def meta(self, path: Path) -> dict:
        """Sidecar metadata of an image (read_metadata), with the raw sidecar under '_sidecar'."""
        if path not in self._meta:
            m = read_metadata(path)["metadata"]
            sidecar = get_sidecar_path(path)
            m["_sidecar"] = sidecar if sidecar.exists() else None
            self._meta[path] = m
        return self._meta[path]

    def sidecar(self, path: Path) -> dict:
        """Raw sidecar content (keys not read by read_metadata, e.g. decay_correction_reference)."""
        s = self.meta(path)["_sidecar"]
        try:
            return json.loads(s.read_text()) if s else {}
        except (OSError, ValueError):
            return {}

    def image(self, path: Path) -> sitk.Image:
        """SimpleITK image (float32), read once."""
        if path not in self._images:
            self._images[path] = sitk.ReadImage(str(path), sitk.sitkFloat32)
        return self._images[path]

    def array(self, path: Path) -> tuple:
        """Voxel values (float64 numpy array) and voxel volume (mL), read once."""
        if path not in self._arrays:
            img = self.image(path)
            voxel_ml = float(np.prod(img.GetSpacing())) / 1000.0
            self._arrays[path] = (sitk.GetArrayViewFromImage(img).astype(np.float64), voxel_ml)
        return self._arrays[path]

    def labels_map(self) -> Dict[int, str]:
        """Label id -> ROI name, from the labels sidecar ({} if none)."""
        s = get_sidecar_path(self.labels)
        try:
            return {int(k): str(v) for k, v in json.loads(s.read_text()).items()}
        except (OSError, ValueError, AttributeError):
            return {}

    def roi_stats(self) -> Dict[int, dict]:
        """
        Per label id: voxels, volume_ml, centroid (physical, LPS mm), border (touches the image border),
        mean_hu (if a CT is given, CT resampled on the labels grid). Computed once.
        """
        if self._roi_stats is not None:
            return self._roi_stats
        img = self.image(self.labels)
        lab = sitk.GetArrayViewFromImage(img).astype(np.int64)
        voxel_ml = float(np.prod(img.GetSpacing())) / 1000.0
        idx = np.nonzero(lab)
        ids = lab[idx]
        n = int(ids.max()) + 1 if ids.size else 1
        counts = np.bincount(ids, minlength=n)
        shape = lab.shape
        border = np.zeros(ids.shape, dtype=bool)
        for axis, i in enumerate(idx):
            border |= (i == 0) | (i == shape[axis] - 1)
        n_border = np.bincount(ids, weights=border, minlength=n)
        centroid_idx = [np.bincount(ids, weights=i, minlength=n) for i in idx]  # z, y, x
        hu_sum = None
        if self.ct:
            ct = sitk.Resample(self.image(self.ct), img, sitk.Transform(), sitk.sitkLinear, -1000.0)
            hu_sum = np.bincount(ids, weights=sitk.GetArrayViewFromImage(ct)[idx], minlength=n)
        stats = {}
        for i in np.nonzero(counts)[0]:
            c = counts[i]
            zyx = [float(ci[i]) / c for ci in centroid_idx]
            point = img.TransformContinuousIndexToPhysicalPoint(zyx[::-1])
            stats[int(i)] = {
                "voxels": int(c),
                "volume_ml": c * voxel_ml,
                "centroid": point,
                "border": bool(n_border[i] > 0),
                "mean_hu": float(hu_sum[i] / c) if hu_sum is not None else None,
            }
        self._roi_stats = stats
        return stats

    def sex(self) -> Optional[str]:
        """Patient sex from the activity sidecars (None if unknown)."""
        for p in self.activity:
            s = (self.meta(p).get("sex") or "").lower()
            if s in ("male", "m"):
                return "male"
            if s in ("female", "f"):
                return "female"
        return None


def _parse_date(value) -> Optional[datetime]:
    try:
        return datetime.strptime(str(value), DATE_FORMAT)
    except (TypeError, ValueError):
        return None


def _float(value) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _time_from_injection_h(m: dict) -> Optional[float]:
    inj, acq = _parse_date(m.get("injection_date")), _parse_date(m.get("acquisition_date"))
    if inj and acq:
        return (acq - inj).total_seconds() / 3600.0
    return None


def _half_life_h(radionuclide) -> Optional[float]:
    from radot.db import get_physical_half_life_h

    try:
        return get_physical_half_life_h(str(radionuclide))
    except Exception:
        return None


# ----------------------------------------------------------------------------------------------
# A. Metadata (sidecars only, no image read)
# ----------------------------------------------------------------------------------------------


def required_fields(unit: Optional[str], method: Optional[str]) -> Dict[str, tuple]:
    """
    Sidecar fields needed for an activity image: {field: (level if missing, reason)}.
    Without a method, the fields of a generic dosimetry are expected (missing = WARN).
    """
    from radot.dose.strategies import TIA_UNITS

    need = Level.FAIL if method else Level.WARN
    unit = (unit or "").lower()
    fields = {
        "unit": (Level.FAIL, "unit of the voxel values"),
        "radionuclide": (Level.WARN, "else give --isotope to radot dose"),
    }
    if unit not in TIA_UNITS:
        reason = "time from injection" + (f", for {method}" if method else "")
        fields["injection_date"] = (need, reason)
        fields["acquisition_date"] = (need, reason)
    if unit in ("suv", "sul"):
        fields["weight"] = (Level.FAIL, f"conversion of {unit} to activity")
        fields["activity"] = (Level.FAIL, f"conversion of {unit} to activity")
    if unit == "sul":
        fields["height"] = (Level.WARN, "lean body mass (james, janmahasatian)")
        fields["sex"] = (Level.WARN, "lean body mass (james, janmahasatian)")
    fields.setdefault("activity", (Level.INFO, "injected activity: activity balance, SUV"))
    fields.setdefault("weight", (Level.INFO, "SUV conversion"))
    fields.setdefault("sex", (Level.INFO, "phantom of the S-values (male by default)"))
    return fields


@register_check("meta.required", ["activity"], "Sidecar fields needed by the dosimetry are present")
def check_meta_required(ctx: CheckContext) -> List[Finding]:
    from radot.dose.strategies import TIA_METHODS, TIA_UNITS

    findings = []
    for p in ctx.activity:
        m = ctx.meta(p)
        if m["_sidecar"] is None:
            findings.append(Finding("meta.required", Level.FAIL, "no sidecar JSON", p.name))
            continue
        unit = (m.get("unit") or "").lower()
        missing = []
        for name, (level, reason) in required_fields(unit, ctx.method).items():
            if m.get(name) in (None, ""):
                missing.append(name)
                findings.append(Finding("meta.required", level, f"{name} missing ({reason})", p.name))
        if ctx.method and unit in TIA_UNITS and ctx.method.lower() not in TIA_METHODS:
            findings.append(
                Finding(
                    "meta.required",
                    Level.FAIL,
                    f"{ctx.method} does not accept a time-integrated activity image ({unit}); use {', '.join(TIA_METHODS)}",
                    p.name,
                )
            )
        if not missing:
            findings.append(Finding("meta.required", Level.OK, "all fields present", p.name))
    return findings


@register_check("meta.modality", [], "Modality and unit of each image are consistent with its role")
def check_meta_modality(ctx: CheckContext) -> List[Finding]:
    findings = []
    roles = [(p, "activity", ("pet", "spect"), ACTIVITY_UNITS) for p in ctx.activity]
    if ctx.ct:
        roles.append((ctx.ct, "ct", ("ct",), ("hu",)))
    if ctx.dose:
        roles.append((ctx.dose, "dose", ("dose",), DOSE_UNITS))
    for p, role, modalities, units in roles:
        m = ctx.meta(p)
        modality, unit = (m.get("modality") or "").lower(), (m.get("unit") or "").lower()
        if modality and modality not in modalities:
            findings.append(
                Finding("meta.modality", Level.FAIL, f"modality '{modality}' for the {role} image", p.name, modality, list(modalities))
            )
        elif unit and unit not in units:
            findings.append(
                Finding("meta.modality", Level.FAIL, f"unit '{unit}' for the {role} image", p.name, unit, list(units))
            )
        else:
            findings.append(Finding("meta.modality", Level.OK, f"{role}: {modality or '?'}, {unit or '?'}", p.name))
    return findings


@register_check("meta.dates_order", ["activity"], "Dates are valid and the injection precedes the acquisition")
def check_meta_dates_order(ctx: CheckContext) -> List[Finding]:
    findings = []
    for p in ctx.activity:
        m = ctx.meta(p)
        bad = [k for k in ("injection_date", "acquisition_date") if m.get(k) and _parse_date(m[k]) is None]
        for k in bad:
            findings.append(
                Finding("meta.dates_order", Level.FAIL, f"{k} '{m[k]}' is not a date", p.name, m[k], "YYYY-MM-DD HH:MM:SS")
            )
        t = _time_from_injection_h(m)
        if t is None:
            continue
        if t < 0:
            findings.append(
                Finding("meta.dates_order", Level.FAIL, f"acquisition {-t:.1f} h before the injection", p.name, round(t, 2), "> 0")
            )
        else:
            findings.append(Finding("meta.dates_order", Level.OK, f"{t:.1f} h after injection", p.name, round(t, 2)))
    return findings


@register_check("meta.time_from_injection", ["activity"], "Time from injection is plausible for the half-life; time points distinct")
def check_meta_time_from_injection(ctx: CheckContext) -> List[Finding]:
    from radot.dose.strategies import TIA_UNITS

    findings = []
    max_hl = ctx.thresholds["time_from_injection"]["max_half_lives"]
    times = []
    for p in ctx.activity:
        m = ctx.meta(p)
        if (m.get("unit") or "").lower() in TIA_UNITS:
            continue
        t = _time_from_injection_h(m)
        if t is None or t < 0:
            continue  # reported by meta.required / meta.dates_order
        times.append((t, p))
        hl = _half_life_h(m.get("radionuclide")) if m.get("radionuclide") else None
        if t == 0:
            findings.append(Finding("meta.time_from_injection", Level.WARN, "acquisition date = injection date", p.name, 0.0))
        elif hl and t > max_hl * hl:
            findings.append(
                Finding(
                    "meta.time_from_injection",
                    Level.WARN,
                    f"{t:.1f} h after injection = {t / hl:.1f} half-lives of {m['radionuclide']} (T½ {hl:.1f} h)",
                    p.name,
                    round(t, 2),
                    f"< {max_hl} half-lives",
                )
            )
        else:
            msg = f"{t:.1f} h" + (f" = {t / hl:.2f} T½" if hl else "")
            findings.append(Finding("meta.time_from_injection", Level.OK, msg, p.name, round(t, 2)))
    seen = {}
    for t, p in times:
        key = round(t, 2)
        if key in seen:
            findings.append(
                Finding("meta.time_from_injection", Level.WARN, f"same time point as {seen[key].name} ({t:.1f} h)", p.name, key)
            )
        seen[key] = p
    return findings


@register_check("meta.patient", ["activity"], "Weight, height, BMI and sex in plausible ranges")
def check_meta_patient(ctx: CheckContext) -> List[Finding]:
    th = ctx.thresholds["patient"]
    findings = []
    for p in ctx.activity:
        m = ctx.meta(p)
        w, h = _float(m.get("weight")), _float(m.get("height"))
        values = [("weight", w, "kg", th["weight_kg"]), ("height", h, "cm", th["height_cm"])]
        if w and h:
            values.append(("BMI", w / (h / 100.0) ** 2, "kg/m²", th["bmi"]))
        bad = False
        for name, v, u, (lo, hi) in values:
            if v is None:
                if m.get(name.lower()) not in (None, ""):
                    findings.append(Finding("meta.patient", Level.WARN, f"{name} '{m[name.lower()]}' is not a number", p.name))
                    bad = True
                continue
            if not lo <= v <= hi:
                findings.append(
                    Finding("meta.patient", Level.WARN, f"{name} {v:.1f} {u} outside [{lo}, {hi}]", p.name, round(v, 2), [lo, hi])
                )
                bad = True
        sex = m.get("sex")
        if sex and sex not in ("male", "female", "m", "f"):
            findings.append(Finding("meta.patient", Level.WARN, f"sex '{sex}' (male or female expected)", p.name, sex))
            bad = True
        if not bad and (w or h or sex):
            given = ", ".join(s for s in (w and f"{w:g} kg", h and f"{h:g} cm", sex) if s)
            findings.append(Finding("meta.patient", Level.OK, given, p.name))
    return findings


@register_check("meta.activity", ["activity"], "Injected activity in the usual range of the radionuclide (MBq)")
def check_meta_activity(ctx: CheckContext) -> List[Finding]:
    from radot.db import normalize_isotope

    ranges = ctx.thresholds["injected_activity_mbq"]
    findings = []
    for p in ctx.activity:
        m = ctx.meta(p)
        a = _float(m.get("activity"))
        if a is None:
            if m.get("activity") not in (None, ""):
                findings.append(Finding("meta.activity", Level.WARN, f"activity '{m['activity']}' is not a number", p.name))
            continue
        rad = normalize_isotope(str(m["radionuclide"])) if m.get("radionuclide") else None
        lo, hi = ranges.get(rad, ranges["default"])
        name = rad if rad in ranges else "any radionuclide"
        if a <= 0 or not lo <= a <= hi:
            findings.append(
                Finding(
                    "meta.activity",
                    Level.WARN,
                    f"injected activity {a:g} MBq outside [{lo:g}, {hi:g}] for {name} (MBq expected, not GBq or Bq)",
                    p.name,
                    a,
                    [lo, hi],
                )
            )
        else:
            findings.append(Finding("meta.activity", Level.OK, f"{a:g} MBq", p.name, a))
    return findings


@register_check("meta.radionuclide", ["activity"], "Radionuclide known to radot (ICRP 107)")
def check_meta_radionuclide(ctx: CheckContext) -> List[Finding]:
    findings = []
    for p in ctx.activity:
        rad = ctx.meta(p).get("radionuclide")
        if not rad:
            continue  # reported by meta.required
        hl = _half_life_h(rad)
        if hl is None:
            findings.append(Finding("meta.radionuclide", Level.FAIL, f"radionuclide '{rad}' unknown (ICRP 107)", p.name, rad))
        else:
            findings.append(Finding("meta.radionuclide", Level.OK, f"{rad}, T½ {hl:.2f} h", p.name, rad))
    return findings


@register_check("meta.consistency", ["activity"], "Same patient data in all the activity images (multi time points)")
def check_meta_consistency(ctx: CheckContext) -> List[Finding]:
    from radot.db import normalize_isotope

    if len(ctx.activity) < 2:
        return []
    findings = []
    keys = ("radionuclide", "injection_date", "activity", "weight", "height", "sex")
    for k in keys:
        values = {}
        for p in ctx.activity:
            v = ctx.meta(p).get(k)
            if v in (None, ""):
                continue
            if k == "radionuclide":
                v = normalize_isotope(str(v))
            values.setdefault(str(v), []).append(p.name)
        if len(values) > 1:
            level = Level.FAIL if k in ("radionuclide", "injection_date") else Level.WARN
            detail = "; ".join(f"{v}: {', '.join(names)}" for v, names in values.items())
            findings.append(Finding("meta.consistency", level, f"{k} differs between the images ({detail})", None, list(values)))
    if not findings:
        findings.append(Finding("meta.consistency", Level.OK, f"{len(ctx.activity)} images, same patient data"))
    return findings


# ----------------------------------------------------------------------------------------------
# B. Units and intensities
# ----------------------------------------------------------------------------------------------


@register_check("units.ct_hu", ["ct"], "CT values in HU: air near -1000, soft tissue near 0 (catches mu maps, +1024 offset)")
def check_units_ct_hu(ctx: CheckContext) -> List[Finding]:
    th = ctx.thresholds["ct"]
    full = sitk.GetArrayViewFromImage(ctx.image(ctx.ct))
    name = ctx.ct.name
    vmax = float(full.max())
    arr = full[::2, ::2, ::2]  # every second voxel: enough for the histogram and the median
    valid = arr[arr >= th["padding_below_hu"]]  # values below = padding outside the field of view
    if valid.size == 0 or vmax < th["padding_below_hu"]:
        return [Finding("units.ct_hu", Level.FAIL, "no value above the padding level", name)]
    vmin, median = float(valid.min()), float(np.median(valid))
    if vmax < th["mu_map_max"]:
        return [
            Finding("units.ct_hu", Level.FAIL, f"values in [{vmin:.3g}, {vmax:.3g}]: not HU (attenuation map?)", name, [vmin, vmax])
        ]
    if vmin > -200 and median > 500:
        return [
            Finding("units.ct_hu", Level.FAIL, f"min {vmin:.0f}, median {median:.0f}: HU with a +1024 offset?", name, [vmin, median])
        ]
    findings = []
    lo, hi = th["air_hu"]
    if not (valid < hi).any():
        findings.append(Finding("units.ct_hu", Level.WARN, f"no air in the image (min {vmin:.0f} HU)", name, vmin, [lo, hi]))
    soft = valid[(valid > -500) & (valid < 500)]
    peak = None
    if soft.size:
        counts, edges = np.histogram(soft, bins=100, range=(-500, 500))
        peak = float(edges[np.argmax(counts)] + 5)
        plo, phi = th["soft_tissue_peak_hu"]
        if not plo <= peak <= phi:
            findings.append(
                Finding("units.ct_hu", Level.WARN, f"soft tissue peak at {peak:.0f} HU", name, peak, [plo, phi])
            )
    if vmax > th["max_hu"]:
        findings.append(Finding("units.ct_hu", Level.WARN, f"max {vmax:.0f} HU (metal, or wrong scaling)", name, vmax, th["max_hu"]))
    if not findings:
        msg = f"min {vmin:.0f}, soft tissue peak {peak:.0f}, max {vmax:.0f} HU" if peak is not None else f"min {vmin:.0f}, max {vmax:.0f} HU"
        findings.append(Finding("units.ct_hu", Level.OK, msg, name))
    return findings


@register_check("units.activity_values", ["activity"], "Activity values: not empty, few negative values, not integer counts")
def check_units_activity_values(ctx: CheckContext) -> List[Finding]:
    th = ctx.thresholds["activity_image"]
    findings = []
    for p in ctx.activity:
        arr, _ = ctx.array(p)
        unit = (ctx.meta(p).get("unit") or "").lower()
        if not arr.any():
            findings.append(Finding("units.activity_values", Level.FAIL, "image is empty (all zeros)", p.name))
            continue
        neg = float((arr < 0).mean())
        bad = False
        if neg > th["max_negative_fraction"]:
            findings.append(
                Finding("units.activity_values", Level.WARN, f"{100 * neg:.1f} % negative voxels", p.name, round(neg, 4), th["max_negative_fraction"])
            )
            bad = True
        # integers with a small max: counts (activity concentrations stored as integers are much larger)
        if unit in ("bq", "bqml") and arr.max() <= th["counts_max"] and np.array_equal(arr, np.round(arr)):
            findings.append(
                Finding("units.activity_values", Level.WARN, f"integer values only (max {arr.max():g}): counts, not {unit}?", p.name, float(arr.max()))
            )
            bad = True
        if not bad:
            findings.append(Finding("units.activity_values", Level.OK, f"max {arr.max():.4g} {unit}", p.name, float(arr.max())))
    return findings


@register_check("units.suv_range", ["activity"], "SUV/SUL images: max in a plausible range")
def check_units_suv_range(ctx: CheckContext) -> List[Finding]:
    lo, hi = ctx.thresholds["suv"]["max_suv"]
    findings = []
    for p in ctx.activity:
        unit = (ctx.meta(p).get("unit") or "").lower()
        if unit not in ("suv", "sul"):
            continue
        vmax = float(ctx.array(p)[0].max())
        if not lo <= vmax <= hi:
            findings.append(Finding("units.suv_range", Level.WARN, f"max {unit} {vmax:.3g} outside [{lo}, {hi}]", p.name, vmax, [lo, hi]))
        else:
            findings.append(Finding("units.suv_range", Level.OK, f"max {unit} {vmax:.3g}", p.name, vmax))
    return findings


@register_check("units.dose_range", ["dose"], "Dose image (Gy): not empty, max plausible, no negative values")
def check_units_dose_range(ctx: CheckContext) -> List[Finding]:
    unit = (ctx.meta(ctx.dose).get("unit") or "").lower()
    name = ctx.dose.name
    if unit != "gy":
        return [Finding("units.dose_range", Level.INFO, f"skipped: unit '{unit or '?'}' (only gy is checked)", name)]
    arr, _ = ctx.array(ctx.dose)
    max_gy = ctx.thresholds["dose"]["max_gy"]
    vmax = float(arr.max())
    findings = []
    if vmax <= 0:
        return [Finding("units.dose_range", Level.FAIL, "no positive dose", name, vmax)]
    if vmax > max_gy:
        findings.append(Finding("units.dose_range", Level.WARN, f"max {vmax:.4g} Gy", name, vmax, max_gy))
    if (arr < 0).any():
        findings.append(Finding("units.dose_range", Level.WARN, f"negative doses (min {arr.min():.3g} Gy)", name, float(arr.min())))
    if not findings:
        findings.append(Finding("units.dose_range", Level.OK, f"max {vmax:.4g} Gy", name, vmax))
    return findings


# ----------------------------------------------------------------------------------------------
# C. Activity balance
# ----------------------------------------------------------------------------------------------


def total_activity_mbq(ctx: CheckContext, p: Path) -> Optional[float]:
    """Total activity in the image (MBq; MBq.s for a TIA image), None for suv/sul or unknown units."""
    arr, voxel_ml = ctx.array(p)
    unit = (ctx.meta(p).get("unit") or "").lower()
    if unit in ("bq", "bqs"):
        return float(arr.sum()) / 1e6
    if unit in ("bqml", "bqsml"):
        return float(arr.sum()) * voxel_ml / 1e6
    return None


def _decay_to_acquisition(ctx: CheckContext, p: Path) -> bool:
    """False if the image is decay-corrected to the injection time (DICOM ADMIN, radot convert)."""
    s = ctx.sidecar(p)
    return not (
        str(s.get("decay_correction", "")).lower() == "admin"
        or str(s.get("decay_correction_reference", "")).lower() == "injection_time"
    )


@register_check("activity.fov_fraction", ["activity"], "Total activity in the image vs injected activity (decay-corrected)")
def check_activity_fov_fraction(ctx: CheckContext) -> List[Finding]:
    from radot.dose.strategies import TIA_UNITS

    th = ctx.thresholds["activity_balance"]
    findings = []
    for p in ctx.activity:
        m = ctx.meta(p)
        unit = (m.get("unit") or "").lower()
        a_inj = _float(m.get("activity"))
        if unit == "suv":
            w = _float(m.get("weight"))
            if not w:
                findings.append(Finding("activity.fov_fraction", Level.INFO, "skipped: no weight", p.name))
                continue
            arr, voxel_ml = ctx.array(p)
            fraction = float(arr.sum()) * voxel_ml / (w * 1000.0)
            ref = "of the injected activity (SUV)"
        else:
            total = total_activity_mbq(ctx, p)
            if total is None or not a_inj:
                why = f"unit '{unit or '?'}'" if total is None else "no injected activity"
                findings.append(Finding("activity.fov_fraction", Level.INFO, f"skipped: {why}", p.name))
                continue
            hl = _half_life_h(m.get("radionuclide")) if m.get("radionuclide") else None
            if hl is None:
                findings.append(Finding("activity.fov_fraction", Level.INFO, "skipped: no radionuclide", p.name))
                continue
            if unit in TIA_UNITS:
                # upper bound: all the injected activity, decaying physically, forever
                tia_max = a_inj * hl * 3600.0 / np.log(2)
                fraction = total / tia_max
                level = Level.FAIL if fraction > th["max_fraction"] else Level.OK
                findings.append(
                    Finding(
                        "activity.fov_fraction",
                        level,
                        f"TIA = {100 * fraction:.1f} % of injected activity x physical mean life (max {100 * th['max_fraction']:.0f} %)",
                        p.name,
                        round(fraction, 4),
                    )
                )
                continue
            ref_activity = a_inj
            t = _time_from_injection_h(m)
            if _decay_to_acquisition(ctx, p):
                if t is None:
                    findings.append(Finding("activity.fov_fraction", Level.INFO, "skipped: no time from injection", p.name))
                    continue
                ref_activity = a_inj * 0.5 ** (t / hl)
                ref = f"of the injected activity decayed to {t:.1f} h"
            else:
                ref = "of the injected activity (image decay-corrected to injection)"
            fraction = total / ref_activity
        msg = f"{100 * fraction:.1f} % {ref}"
        if fraction > th["max_fraction"]:
            findings.append(Finding("activity.fov_fraction", Level.FAIL, msg + ": units or calibration?", p.name, round(fraction, 4), [0, th["max_fraction"]]))
        elif fraction < th["min_fraction"]:
            findings.append(Finding("activity.fov_fraction", Level.FAIL, msg + ": uncalibrated counts?", p.name, round(fraction, 4), [th["min_fraction"], th["max_fraction"]]))
        elif fraction < th["warn_below"]:
            findings.append(Finding("activity.fov_fraction", Level.WARN, msg, p.name, round(fraction, 4), [th["warn_below"], th["max_fraction"]]))
        else:
            findings.append(Finding("activity.fov_fraction", Level.OK, msg, p.name, round(fraction, 4)))
    return findings


@register_check("activity.mtp_trend", ["activity"], "Multi time points: total activity decreases at least as fast as the physical decay")
def check_activity_mtp_trend(ctx: CheckContext) -> List[Finding]:
    if len(ctx.activity) < 2:
        return []
    tol = ctx.thresholds["activity_balance"]["mtp_tolerance"]
    points = []
    for p in ctx.activity:
        m = ctx.meta(p)
        t, total = _time_from_injection_h(m), total_activity_mbq(ctx, p)
        hl = _half_life_h(m.get("radionuclide")) if m.get("radionuclide") else None
        if t is None or total is None or hl is None or not _decay_to_acquisition(ctx, p):
            return [Finding("activity.mtp_trend", Level.INFO, "skipped: time, unit (bq, bqml), radionuclide or decay correction", p.name)]
        points.append((t, total, hl, p))
    points.sort(key=lambda x: x[0])
    findings = []
    for (t1, a1, hl, p1), (t2, a2, _, p2) in zip(points, points[1:]):
        physical = 0.5 ** ((t2 - t1) / hl)
        ratio = a2 / a1 if a1 > 0 else float("inf")
        teff = (t2 - t1) * np.log(2) / np.log(a1 / a2) if 0 < a2 < a1 else None
        pair = f"{p1.name} -> {p2.name}"
        if ratio > physical * (1 + tol):
            findings.append(
                Finding(
                    "activity.mtp_trend",
                    Level.WARN,
                    f"total activity x{ratio:.2f} in {t2 - t1:.1f} h, physical decay alone gives x{physical:.2f}: "
                    "different calibration, or images decay-corrected?",
                    pair,
                    round(ratio, 4),
                    round(physical, 4),
                )
            )
        else:
            teff_txt = f", Teff {teff:.1f} h (T½ {hl:.1f} h)" if teff else ""
            findings.append(Finding("activity.mtp_trend", Level.OK, f"total activity x{ratio:.2f} in {t2 - t1:.1f} h{teff_txt}", pair, round(ratio, 4)))
    return findings


# ----------------------------------------------------------------------------------------------
# D. Geometry
# ----------------------------------------------------------------------------------------------


def _on_grid(image: sitk.Image, reference: sitk.Image, nearest: bool = False, default: float = 0.0) -> np.ndarray:
    """Image resampled on the grid of the reference (physical space), as a numpy array."""
    interp = sitk.sitkNearestNeighbor if nearest else sitk.sitkLinear
    return sitk.GetArrayFromImage(sitk.Resample(image, reference, sitk.Transform(), interp, default))


def _inside(image: sitk.Image, reference: sitk.Image, factor: int = 8) -> np.ndarray:
    """
    Boolean array on the reference grid: voxels inside the field of view of image. The box is built
    with voxels `factor` times larger (same physical extent), to keep it cheap for large CTs.
    """
    size = np.array(image.GetSize())
    spacing = np.array(image.GetSpacing())
    coarse_size = np.ceil(size / factor).astype(int)
    coarse_spacing = spacing * size / coarse_size
    direction = np.array(image.GetDirection()).reshape(3, 3)
    origin = np.array(image.GetOrigin()) + direction @ ((coarse_spacing - spacing) / 2)
    ones = sitk.Image([int(s) for s in coarse_size], sitk.sitkUInt8) + 1
    ones.SetSpacing([float(s) for s in coarse_spacing])
    ones.SetOrigin([float(o) for o in origin])
    ones.SetDirection(image.GetDirection())
    return _on_grid(ones, reference, nearest=True) > 0


def _describe_shift(v) -> str:
    """LPS vector (mm) -> '12 mm posterior, 5 mm left'."""
    names = (("left", "right"), ("posterior", "anterior"), ("superior", "inferior"))
    parts = sorted(((abs(c), n[0] if c > 0 else n[1]) for c, n in zip(v, names) if abs(c) >= 1), reverse=True)
    return ", ".join(f"{a:.0f} mm {n}" for a, n in parts) or "0 mm"


@register_check("geom.direction", [], "Image axes aligned with the patient axes (no oblique images)")
def check_geom_direction(ctx: CheckContext) -> List[Finding]:
    findings = []
    for p in ctx.activity + [x for x in (ctx.ct, ctx.labels, ctx.dose) if x]:
        r = sitk.ImageFileReader()
        r.SetFileName(str(p))
        r.ReadImageInformation()
        d = np.abs(np.array(r.GetDirection()).reshape(3, 3))
        if not np.allclose(np.sort(d, axis=1)[:, -1], 1.0, atol=1e-3):
            findings.append(Finding("geom.direction", Level.WARN, "oblique image (direction not axis-aligned)", p.name, np.round(d, 3).tolist()))
    if not findings:
        findings.append(Finding("geom.direction", Level.OK, "all images axis-aligned"))
    return findings


@register_check("geom.overlap", [], "Activity inside the CT field of view; organs inside the activity field of view")
def check_geom_overlap(ctx: CheckContext) -> List[Finding]:
    th = ctx.thresholds["geom"]
    findings = []
    for p in ctx.activity:
        act = ctx.image(p)
        a = np.clip(ctx.array(p)[0], 0, None)
        total = a.sum()
        if ctx.ct and total > 0:
            frac = float(a[_inside(ctx.image(ctx.ct), act)].sum() / total)
            level = Level.FAIL if frac < th["overlap_fail_below"] else Level.WARN if frac < th["overlap_warn_below"] else Level.OK
            findings.append(Finding("geom.overlap", level, f"{100 * frac:.1f} % of the activity inside the CT field of view", p.name, round(frac, 4)))
        if ctx.labels:
            # organs partly outside the activity image: activity (dose) underestimated
            lab = _on_grid(ctx.image(ctx.labels), act, nearest=True, default=-1).astype(np.int64)
            labels_map, stats = ctx.labels_map(), ctx.roi_stats()
            cut = []
            for i, roi in sorted(labels_map.items()):
                if i not in stats or _roi_class(roi, {"organ": th["organs"]}) is None:
                    continue
                n_in = int((lab == i).sum())
                full = stats[i]["voxels"] * ctx.array(ctx.labels)[1] / ctx.array(p)[1]  # voxels on the activity grid
                if n_in < th["organ_inside_warn_below"] * full:
                    cut.append(f"{roi} {100 * n_in / full:.0f} %")
            if cut:
                findings.append(
                    Finding("geom.overlap", Level.WARN, f"organs partly outside the activity image (fraction inside): {', '.join(cut)}", p.name, cut)
                )
    return findings


@register_check("geom.activity_in_body", ["activity", "ct"], "Activity inside the CT body outline (catches large shifts between SPECT/PET and CT)")
def check_geom_activity_in_body(ctx: CheckContext) -> List[Finding]:
    th = ctx.thresholds["geom"]
    findings = []
    for p in ctx.activity:
        act = ctx.image(p)
        a = np.clip(ctx.array(p)[0], 0, None)
        ct = _on_grid(ctx.image(ctx.ct), act, default=-1000.0)
        body = sitk.GetImageFromArray((ct > -500).astype(np.uint8))
        body.CopyInformation(act)
        body = sitk.RelabelComponent(sitk.ConnectedComponent(body), sortByObjectSize=True) == 1
        radius = [int(np.ceil(th["body_dilation_mm"] / s)) for s in act.GetSpacing()]
        body = sitk.GetArrayFromImage(sitk.BinaryDilate(body, radius)) > 0
        inside_ct = _inside(ctx.image(ctx.ct), act)
        total = a[inside_ct].sum()
        if total <= 0:
            continue
        frac = float(a[body & inside_ct].sum() / total)
        level = Level.FAIL if frac < th["body_fail_below"] else Level.WARN if frac < th["body_warn_below"] else Level.OK
        msg = f"{100 * frac:.1f} % of the activity inside the CT body (+{th['body_dilation_mm']:g} mm)"
        if level > Level.OK:
            msg += ": SPECT/PET and CT not registered?"
        findings.append(Finding("geom.activity_in_body", level, msg, p.name, round(frac, 4)))
    return findings


def best_shift(a: np.ndarray, mask: np.ndarray, n: Sequence[int]) -> tuple:
    """
    Shift (index, z y x) of the mask, within +-n voxels, that maximises the activity inside it
    (cross-correlation by FFT on the bounding box of the mask padded by n). Returns
    (shift, counts at the best shift, counts at no shift).
    """
    idx = np.nonzero(mask)
    lo = [max(int(i.min()) - k, 0) for i, k in zip(idx, n)]
    hi = [min(int(i.max()) + k + 1, s) for i, k, s in zip(idx, n, a.shape)]
    box = tuple(slice(l, h) for l, h in zip(lo, hi))
    ab, mb = a[box], mask[box].astype(np.float64)
    shape = [s + k for s, k in zip(ab.shape, n)]  # extra zero padding: no wrap-around
    axes = (0, 1, 2)
    corr = np.fft.irfftn(np.fft.rfftn(ab, shape, axes) * np.conj(np.fft.rfftn(mb, shape, axes)), shape, axes)
    best, best_shift_idx = -np.inf, (0, 0, 0)
    for dz in range(-n[0], n[0] + 1):
        for dy in range(-n[1], n[1] + 1):
            for dx in range(-n[2], n[2] + 1):
                v = corr[dz % shape[0], dy % shape[1], dx % shape[2]]
                if v > best:
                    best, best_shift_idx = v, (dz, dy, dx)
    return best_shift_idx, float(best), float(corr[0, 0, 0])


@register_check("geom.organ_shift", ["activity", "labels"], "Uptake organs: shift of the label that maximises its activity (misregistration)")
def check_geom_organ_shift(ctx: CheckContext) -> List[Finding]:
    th = ctx.thresholds["geom"]
    findings = []
    labels_map = ctx.labels_map()
    for p in ctx.activity:
        act = ctx.image(p)
        a = np.clip(ctx.array(p)[0], 0, None)
        if a.max() <= 0:
            continue
        lab = _on_grid(ctx.image(ctx.labels), act, nearest=True).astype(np.int64)
        background = a[a > 0.01 * a.max()].mean()
        spacing = act.GetSpacing()  # x y z
        n = [int(th["shift_search_mm"] // s) for s in spacing[::-1]]  # z y x
        direction = np.array(act.GetDirection()).reshape(3, 3)
        shifted, ok = [], []
        for i, roi in sorted(labels_map.items()):
            if _roi_class(roi, {"organ": th["organs"]}) is None:
                continue
            mask = lab == i
            if mask.sum() < 8:
                continue
            if a[mask].mean() < th["min_uptake"] * background:
                continue  # low uptake: the best shift would be noise
            s_idx, best, zero = best_shift(a, mask, n)
            mm = direction @ (np.array(s_idx[::-1]) * spacing)  # LPS
            dist = float(np.linalg.norm(mm))
            gain = best / zero - 1 if zero > 0 else float("inf")
            if dist > th["shift_warn_mm"] and gain > th["shift_min_gain"]:
                edge = " (search limit)" if any(abs(s) == k for s, k in zip(s_idx, n)) else ""
                shifted.append(f"{roi} {_describe_shift(mm)}{edge}, +{100 * gain:.0f} % counts")
            else:
                ok.append(f"{roi} {dist:.0f} mm")
        if shifted:
            findings.append(
                Finding(
                    "geom.organ_shift",
                    Level.WARN,
                    f"labels not on the activity, best match: {'; '.join(shifted)}. Labels and activity image not registered?",
                    p.name,
                    shifted,
                )
            )
        if ok:
            findings.append(Finding("geom.organ_shift", Level.OK, f"best match shift: {', '.join(ok)}", p.name, ok))
    return findings


@register_check("activity.lr_symmetry", ["activity", "labels"], "Left/right uptake organs (kidneys): similar activity concentration")
def check_activity_lr_symmetry(ctx: CheckContext) -> List[Finding]:
    th = ctx.thresholds["geom"]
    min_ratio = ctx.thresholds["activity_balance"]["lr_min_ratio"]
    labels_map = ctx.labels_map()
    pairs = {}
    for i, roi in labels_map.items():
        side = roi_side(roi)
        if side and _roi_class(roi, {"organ": th["organs"]}):
            pairs.setdefault(_roi_base(roi), {})[side] = (i, roi)
    pairs = {b: s for b, s in pairs.items() if len(s) == 2}
    if not pairs:
        return []
    findings = []
    for p in ctx.activity:
        a = ctx.array(p)[0]
        lab = _on_grid(ctx.image(ctx.labels), ctx.image(p), nearest=True).astype(np.int64)
        for base, sides in sorted(pairs.items()):
            (il, rl), (ir, rr) = sides["L"], sides["R"]
            ml, mr = (lab == il), (lab == ir)
            if ml.sum() == 0 or mr.sum() == 0:
                continue
            cl, cr = float(a[ml].mean()), float(a[mr].mean())
            if max(cl, cr) <= 0:
                continue
            ratio = min(cl, cr) / max(cl, cr)
            low = rl if cl < cr else rr
            msg = f"{rl}/{rr} mean activity ratio {cl / cr:.2f}" if cr > 0 else f"{rr}: no activity"
            if ratio < min_ratio:
                findings.append(
                    Finding("activity.lr_symmetry", Level.WARN, f"{msg}: low uptake in {low} (non-functioning organ, or labels not registered?)", p.name, round(cl / cr, 3) if cr > 0 else None, min_ratio)
                )
            else:
                findings.append(Finding("activity.lr_symmetry", Level.OK, msg, p.name, round(cl / cr, 3)))
    return findings


# ----------------------------------------------------------------------------------------------
# E. ROIs
# ----------------------------------------------------------------------------------------------


def roi_side(name: str) -> Optional[str]:
    """'L' or 'R' from a ROI name (TG-263 'Kidney_L', TotalSegmentator 'kidney_left', 'rib_left_3')."""
    n = name.lower()
    if re.search(r"(^|[_\s-])(l|left)($|[_\s-])", n):
        return "L"
    if re.search(r"(^|[_\s-])(r|right)($|[_\s-])", n):
        return "R"
    return None


def _roi_base(name: str) -> str:
    """Name without the side, to pair 'Kidney_L' with 'Kidney_R'."""
    return re.sub(r"(^|[_\s-])(l|left|r|right)(?=$|[_\s-])", r"\1", name.lower()).strip("_ -")


def _roi_class(name: str, classes: dict) -> Optional[str]:
    n = name.lower()
    for cls, words in classes.items():
        if any(w in n for w in words):
            return cls
    return None


@register_check("roi.names", ["labels"], "Labels sidecar present; ROI names resolve to OpenDose organs")
def check_roi_names(ctx: CheckContext) -> List[Finding]:
    from radot.nomenclature import verify_sidecar_mappings

    name = ctx.labels.name
    labels_map = ctx.labels_map()
    if not labels_map:
        return [Finding("roi.names", Level.WARN, "no labels sidecar (label id -> ROI name)", name)]
    report = verify_sidecar_mappings(labels_map)
    findings = [Finding("roi.names", Level.OK, f"{len(report['mapped'])} of {len(labels_map)} ROIs mapped to OpenDose", name)]
    if report["unmapped"]:
        unmapped = sorted(report["unmapped"].values())
        findings.append(
            Finding(
                "roi.names",
                Level.INFO,
                f"{len(unmapped)} ROIs without OpenDose organ (no S-value, skipped by organ-level methods): {', '.join(unmapped)}",
                name,
                unmapped,
            )
        )
    return findings


@register_check("roi.empty", ["labels"], "ROIs of the sidecar present in the image, not tiny, not truncated by the image border")
def check_roi_empty(ctx: CheckContext) -> List[Finding]:
    min_ml = ctx.thresholds["roi"]["min_volume_ml"]
    labels_map, stats = ctx.labels_map(), ctx.roi_stats()
    name = ctx.labels.name
    if not stats:
        return [Finding("roi.empty", Level.FAIL, "labels image is empty", name)]
    # On a cropped image, empty or truncated ROIs of a full segmentation are normal: reported, grouped
    empty = [roi for i, roi in sorted(labels_map.items()) if i not in stats]
    tiny = [roi for i, roi in sorted(labels_map.items()) if i in stats and stats[i]["volume_ml"] < min_ml]
    border = [roi for i, roi in sorted(labels_map.items()) if i in stats and stats[i]["border"] and roi not in tiny]
    extra = sorted(set(stats) - set(labels_map)) if labels_map else []
    findings = [Finding("roi.empty", Level.OK, f"{len(stats)} labels in the image", name)]
    for rois, what in (
        (empty, "empty"),
        (tiny, f"smaller than {min_ml:g} mL"),
        (border, "touching the image border (truncated?)"),
    ):
        if rois:
            findings.append(Finding("roi.empty", Level.INFO, f"{len(rois)} ROIs {what}: {', '.join(rois)}", name, rois))
    if extra:
        findings.append(Finding("roi.empty", Level.INFO, f"label ids in the image but not in the sidecar: {extra}", name, extra))
    return findings


@register_check("roi.laterality", ["labels"], "Left/right ROIs on the correct side (LPS: patient left = +x); liver on the right")
def check_roi_laterality(ctx: CheckContext) -> List[Finding]:
    labels_map, stats = ctx.labels_map(), ctx.roi_stats()
    rois = {roi: stats[i] for i, roi in labels_map.items() if i in stats}
    if not rois:
        return []
    findings = []
    # pairs: x(left) > x(right)
    pairs = {}
    for roi in rois:
        side = roi_side(roi)
        if side:
            pairs.setdefault(_roi_base(roi), {})[side] = roi
    full = [s for s in pairs.values() if len(s) == 2]
    swapped = [s["L"] for s in full if rois[s["L"]]["centroid"][0] <= rois[s["R"]]["centroid"][0]]
    if swapped:
        findings.append(
            Finding(
                "roi.laterality",
                Level.FAIL,
                f"{len(swapped)} of {len(full)} left/right pairs swapped (left ROI on the patient right): {', '.join(sorted(swapped))}",
                ctx.labels.name,
                sorted(swapped),
            )
        )
    # organs on a known side vs the centre of all the labels (catches a mirrored image)
    centre = np.average([s["centroid"][0] for s in rois.values()], weights=[s["voxels"] for s in rois.values()])
    sides_known = ctx.thresholds["roi"]["organ_side"]
    wrong = []
    for roi, s in sorted(rois.items()):
        expected = sides_known.get(roi.lower())
        x = s["centroid"][0]
        if (expected == "R" and x >= centre) or (expected == "L" and x <= centre):
            wrong.append(f"{roi} on the patient {'left' if expected == 'R' else 'right'}")
    if wrong:
        findings.append(Finding("roi.laterality", Level.FAIL, f"{', '.join(wrong)}: image or labels mirrored?", ctx.labels.name, wrong))
    if not findings:
        findings.append(Finding("roi.laterality", Level.OK, f"{len(full)} left/right pairs consistent", ctx.labels.name))
    return findings


@register_check("roi.mass", ["labels", "ct"], "Mass (CT density) of whole organs vs the OpenDose reference mass")
def check_roi_mass(ctx: CheckContext) -> List[Finding]:
    from radot.db import get_reference_mass, resolve_organ_name

    # only whole organs listed in the thresholds (lobes, vessels, muscles are parts of an OpenDose organ)
    th = ctx.thresholds["roi"]["mass_ratio"]
    sex = ctx.sex()
    phantom = sex or "male"
    labels_map, stats = ctx.labels_map(), ctx.roi_stats()
    organs = {}
    for i, roi in labels_map.items():
        if i in stats:
            organs.setdefault(resolve_organ_name(roi, phantom).lower(), []).append((i, roi))
    findings, ratios = [], []
    for organ, rois in sorted(organs.items()):
        if organ not in th or len(rois) != 1:
            continue
        i, roi = rois[0]
        s = stats[i]
        try:
            ref = get_reference_mass(phantom, roi)
        except Exception:
            continue
        mass = s["volume_ml"] * max(0.0, 1.0 + s["mean_hu"] / 1000.0)
        ratio = mass / ref
        lo, hi = th[organ]
        if s["border"]:
            findings.append(Finding("roi.mass", Level.INFO, f"{mass:.0f} g, ratio {ratio:.2f}: not checked, touches the image border", roi, round(ratio, 3)))
        elif not lo <= ratio <= hi:
            findings.append(
                Finding("roi.mass", Level.WARN, f"{mass:.0f} g, reference {ref:.0f} g ({phantom}), ratio {ratio:.2f}", roi, round(ratio, 3), [lo, hi])
            )
        else:
            ratios.append(ratio)
    if ratios:
        msg = f"{len(ratios)} organs, mass {min(ratios):.2f}-{max(ratios):.2f} x OpenDose ({phantom}{'' if sex else ', sex unknown'})"
        findings.append(Finding("roi.mass", Level.OK, msg, ctx.labels.name))
    return findings


@register_check("roi.density", ["labels", "ct"], "Mean CT value of each ROI plausible for its tissue (lung, bone, soft tissue)")
def check_roi_density(ctx: CheckContext) -> List[Finding]:
    th = ctx.thresholds["roi"]
    labels_map, stats = ctx.labels_map(), ctx.roi_stats()
    n, bad = 0, {}
    for i, roi in sorted(labels_map.items()):
        s = stats.get(i)
        cls = _roi_class(roi, th["tissue_classes"])
        if s is None or cls is None:
            continue
        n += 1
        lo, hi = th["mean_hu"][cls]
        if not lo <= s["mean_hu"] <= hi:
            bad[roi] = round(s["mean_hu"], 1)
    if not n:
        return []
    if not bad:
        return [Finding("roi.density", Level.OK, f"{n} ROIs with a plausible mean HU", ctx.labels.name)]
    ranges = ", ".join(f"{c} {lo}..{hi}" for c, (lo, hi) in th["mean_hu"].items())
    detail = ", ".join(f"{roi} {hu:.0f}" for roi, hu in bad.items())
    return [
        Finding(
            "roi.density",
            Level.WARN,
            f"{len(bad)} of {n} ROIs with an implausible mean HU (expected {ranges}): labels misaligned with the CT? {detail}",
            ctx.labels.name,
            bad,
        )
    ]


# ----------------------------------------------------------------------------------------------
# Run
# ----------------------------------------------------------------------------------------------


def _selected(code: str, patterns: Optional[Sequence[str]]) -> bool:
    """A pattern is a code ('meta.patient') or a category ('meta')."""
    return any(code == p or code.startswith(p + ".") for p in patterns or [])


def run_checks(
    activity=None,
    ct=None,
    labels=None,
    dose=None,
    method: Optional[str] = None,
    thresholds=None,
    only: Optional[Sequence[str]] = None,
    skip: Optional[Sequence[str]] = None,
) -> List[Finding]:
    """
    Run the checks whose inputs are given and return the findings.
    activity: one image path or a list (multi time points); ct, labels, dose: image paths.
    method: dosimetry method planned (fields required); thresholds: dict or JSON path (merged with
    the defaults); only / skip: codes or categories.
    """
    if activity is None:
        activity = []
    elif isinstance(activity, (str, Path)):
        activity = [activity]
    if method:
        from radot.dose.strategies import DOSIMETRY_STRATEGIES

        known = list(DOSIMETRY_STRATEGIES) + ["ldm", "jackson2013"]
        if method.lower() not in known:
            raise ValueError(f"Unknown method '{method}'. Available: {known}")
    if isinstance(thresholds, dict):
        th = _deep_merge(load_thresholds(), copy.deepcopy(thresholds))
    else:
        th = load_thresholds(thresholds)
    ctx = CheckContext(
        activity=[Path(p) for p in activity],
        ct=Path(ct) if ct else None,
        labels=Path(labels) if labels else None,
        dose=Path(dose) if dose else None,
        method=method,
        thresholds=th,
    )
    for p in ctx.activity + [x for x in (ctx.ct, ctx.labels, ctx.dose) if x]:
        if not p.exists():
            raise FileNotFoundError(f"Image file not found: {p}")
    available = ctx.available()
    findings = []
    for code, check in CHECKS.items():
        if only and not _selected(code, only):
            continue
        if skip and _selected(code, skip):
            continue
        if not set(check.inputs) <= available:
            continue
        findings.extend(check.func(ctx))
    order = {c: i for i, c in enumerate(CATEGORIES)}
    findings.sort(key=lambda f: order.get(f.category, len(order)))
    return findings


def worst_level(findings: Sequence[Finding]) -> Level:
    return max((f.level for f in findings), default=Level.OK)
