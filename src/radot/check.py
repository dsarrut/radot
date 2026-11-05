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
from dataclasses import dataclass, field
from datetime import datetime
from enum import IntEnum
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence

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
