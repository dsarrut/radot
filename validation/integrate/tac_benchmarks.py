"""
Time integration (radot integrate) on synthetic ¹⁷⁷Lu time-activity curves with a known integral.

Single source of truth of the cases, the methods, the Monte Carlo and the acceptance ranges: used by
test_integrate_validation.py (CI, reduced number of realisations) and report.py (full run, stored in
results/integrate_results.json, page docs/9_validation_time_integration.md).

Each curve is sampled at the times of a schedule, with a proportional Gaussian noise (a fractional
standard deviation, the usual model of SPECT ROI activities). Every method sees the same noisy
realisations (common random numbers). Per configuration: relative error of the TIA (mean, SD, RMSE),
the failures, and for the methods with a standard error, its coverage (the fraction of realisations
with |TIA - true| <= SE, expected 68 % when the SE is right).
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np

from radot.db import get_physical_half_life_h
from radot.integrate import TAC_MODELS, integrate_curve

ISOTOPE = "Lu177"
LAMBDA_P = math.log(2.0) / get_physical_half_life_h(ISOTOPE)
LN2 = math.log(2.0)


def _rate(teff_h: float) -> float:
    """Biological rate of an effective half-life."""
    return LN2 / teff_h - LAMBDA_P


@dataclass(frozen=True)
class Shape:
    """A ¹⁷⁷Lu time-activity curve: a Kletting model with given parameters (MBq, 1/h)."""

    name: str
    model: str
    params: Tuple[float, ...]
    description: str

    def curve(self, t) -> np.ndarray:
        from radot.integrate import evaluate

        return evaluate(TAC_MODELS[self.model], np.asarray(self.params, float), np.asarray(t, float), LAMBDA_P)

    def tia_mbq_h(self) -> float:
        from radot.integrate.models import tia_and_gradient

        return tia_and_gradient(TAC_MODELS[self.model], np.asarray(self.params, float), LAMBDA_P)[0]


# Typical ¹⁷⁷Lu-PSMA / DOTATATE kinetics (effective half-lives as in the EANM 2022 recommendations):
# kidney uptake in ~1 h then Teff ≈ 40 h; liver: a fast blood-pool part then Teff ≈ 70 h; lesion: slow
# uptake (~7 h) then Teff ≈ 90 h.
SHAPES: Dict[str, Shape] = {
    s.name: s
    for s in (
        Shape("kidney", "f4a", (10.0, _rate(40.0), LN2 / 0.7), "uptake ~1 h, Teff 40 h (f4a)"),
        Shape("liver", "f4", (30.0, 70.0, 0.2, _rate(70.0)), "30 % fast (λ 0.2/h), 70 % Teff 70 h (f4)"),
        Shape("lesion", "f4a", (5.0, _rate(90.0), LN2 / 7.0), "uptake ~7 h, Teff 90 h (f4a)"),
    )
}

# Acquisition schedules (h after injection)
SCHEDULES: Dict[str, Tuple[float, ...]] = {
    "3tp": (4.0, 24.0, 96.0),
    "4tp": (4.0, 24.0, 96.0, 168.0),
    "6tp": (2.0, 20.0, 44.0, 68.0, 116.0, 164.0),  # Ivashchenko 2024
}
STP_TIME_H = 96.0  # single time point: the scan of the schedule closest to 96 h

NOISES = (0.05, 0.10)  # fractional standard deviation of each point


@dataclass(frozen=True)
class Method:
    name: str  # label in the report
    model: str
    options: Tuple[Tuple[str, object], ...] = ()
    known_noise: bool = False  # error model fsd:<the true noise> (absolute weighting)

    def kwargs(self, noise: float) -> dict:
        kw = dict(self.options)
        if self.known_noise:
            kw["error_model"] = f"fsd:{noise:g}"
        return kw


METHODS: Dict[str, Method] = {
    m.name: m
    for m in (
        Method("f2 proportional", "f2", (("error_model", "proportional"),)),
        Method("f2 fsd", "f2", known_noise=True),
        Method("f4a proportional", "f4a", (("error_model", "proportional"),)),
        Method("f4a fsd", "f4a", known_noise=True),
        Method("auto best", "auto", (("candidates", ("f2", "f4a")), ("select", "best")), known_noise=True),
        Method("auto average", "auto", (("candidates", ("f2", "f4a")), ("select", "average")), known_noise=True),
        Method("trapezoid ramp/phys", "trapezoid", (("head", "ramp"), ("tail", "phys"))),
        Method("trapezoid extrap/fit2", "trapezoid", (("head", "extrap"), ("tail", "fit"), ("tail_points", 2))),
        Method("hanscheid2018 (scan ~96 h)", "hanscheid2018"),
    )
}


def stp_index(times: Tuple[float, ...]) -> int:
    return int(np.argmin(np.abs(np.asarray(times) - STP_TIME_H)))


def applicable(method: Method, n_points: int) -> Optional[str]:
    """None when the method can run on n_points, else the reason (Kletting point rules)."""
    from radot.integrate import check_integration_options

    n = 1 if method.model == "hanscheid2018" else n_points
    try:
        check_integration_options(method.model, n, **method.kwargs(0.1))
    except ValueError as e:
        return str(e)
    return None


def realisations(shape: Shape, schedule: str, noise: float, n: int) -> np.ndarray:
    """(n, N) noisy activities; the seed depends only on the shape, schedule and noise."""
    keys = (list(SHAPES).index(shape.name), list(SCHEDULES).index(schedule), int(round(noise * 1000)))
    rng = np.random.default_rng([2026, *keys])
    y = shape.curve(SCHEDULES[schedule])
    return y[None, :] * (1.0 + noise * rng.standard_normal((n, len(y))))


def run_one(method: Method, times, values, noise: float) -> dict:
    kw = method.kwargs(noise)
    if method.model == "hanscheid2018":
        i = stp_index(times)
        return integrate_curve([times[i]], [values[i]], LAMBDA_P, "hanscheid2018",
                               teff_h=2.0 * times[i], teff_source="default")
    return integrate_curve(list(times), list(values), LAMBDA_P, method.model, **kw)


def run_config(shape: str, schedule: str, noise: float, method: str, n: int) -> dict:
    """Statistics of one configuration over n realisations (relative errors in %)."""
    s, m = SHAPES[shape], METHODS[method]
    times = SCHEDULES[schedule]
    reason = applicable(m, len(times))
    if reason:
        return {"applicable": False, "reason": reason}
    truth = s.tia_mbq_h()
    errors, z, failures = [], [], 0
    for values in realisations(s, schedule, noise, n):
        try:
            r = run_one(m, times, values, noise)
        except ValueError:
            failures += 1
            continue
        errors.append(r["tia_mbq_h"] / truth - 1.0)
        if r["tia_se_mbq_h"]:
            z.append((r["tia_mbq_h"] - truth) / r["tia_se_mbq_h"])
    e = 100.0 * np.asarray(errors)
    out = {
        "applicable": True,
        "n": n,
        "failures": failures,
        "bias_pct": float(np.mean(e)) if len(e) else None,
        "sd_pct": float(np.std(e, ddof=1)) if len(e) > 1 else None,
        "rmse_pct": float(np.sqrt(np.mean(e * e))) if len(e) else None,
        "median_pct": float(np.median(e)) if len(e) else None,
    }
    if z:
        za = np.abs(np.asarray(z))
        out["coverage_68"] = float(np.mean(za <= 1.0))
        out["coverage_95"] = float(np.mean(za <= 1.96))
        out["n_se"] = len(z)
    return out


def noise_free(shape: str, schedule: str, method: str):
    """Relative error (%) of the TIA on the exact points (no noise), None when not applicable, "failed"
    when the integration fails."""
    s, m = SHAPES[shape], METHODS[method]
    times = SCHEDULES[schedule]
    if applicable(m, len(times)):
        return None
    try:
        r = run_one(m, times, s.curve(times), 0.05)
    except ValueError:
        return "failed"
    return 100.0 * (r["tia_mbq_h"] / s.tia_mbq_h() - 1.0)


def run_all(n: int, shapes=None, schedules=None, noises=None, methods=None, progress=None) -> dict:
    """Every configuration: {"settings", "truth", "noise_free", "results": {shape: {schedule: {noise:
    {method: stats}}}}}."""
    shapes = shapes or list(SHAPES)
    schedules = schedules or list(SCHEDULES)
    noises = noises or list(NOISES)
    methods = methods or list(METHODS)
    res: dict = {}
    nf: dict = {}
    for sh in shapes:
        for sc in schedules:
            nf.setdefault(sh, {})[sc] = {m: noise_free(sh, sc, m) for m in methods}
            for no in noises:
                for m in methods:
                    res.setdefault(sh, {}).setdefault(sc, {}).setdefault(f"{no:g}", {})[m] = run_config(sh, sc, no, m, n)
                    if progress:
                        progress(sh, sc, no, m)
    return {
        "settings": {"n_realisations": n, "isotope": ISOTOPE, "stp_time_h": STP_TIME_H,
                     "schedules": {k: list(v) for k, v in SCHEDULES.items()},
                     "noises": list(noises),
                     "shapes": {k: {"model": v.model, "params": list(v.params), "description": v.description}
                                for k, v in SHAPES.items()}},
        "truth_mbq_h": {k: SHAPES[k].tia_mbq_h() for k in shapes},
        "noise_free_pct": nf,
        "results": res,
    }


# ----------------------------------------------------------------------------------------------
# acceptance (asserted by test_integrate_validation.py)
# ----------------------------------------------------------------------------------------------

# Noise-free: a model that contains the true shape recovers the integral (%)
NOISE_FREE_EXACT_PCT = 1e-4
# With noise, the true model with the known noise (fsd) at 6 time points: unbiased and an honest SE
TRUE_MODEL_BIAS_PCT = 2.0
COVERAGE_68_RANGE = (0.55, 0.80)
COVERAGE_95_RANGE = (0.88, 0.99)
