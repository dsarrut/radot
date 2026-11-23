"""Time integration without a fitted function (docs/6e_time_integration.md §3.3–3.4):

- `trapezoid`: through the measured points, with a head (0 -> t1) and a tail (tN -> infinity);
- single time point models (one image): `physical`, `hanscheid2018`, `madsen2018`, the time factors used
  inside `radot dose`, so that `integrate --model hanscheid2018` + `dose -m svalue` = `dose -m hanscheid2018`.

Also the option check of every model kind (check_integration_options) and the dispatch
(integrate_curve). Units: t in h, activities in MBq, TIA in MBq.h.
"""

import logging
import math
from typing import Any, Dict, Optional, Sequence

import numpy as np

from radot.integrate.models import (
    TAC_MODELS,
    check_model_choice,
    integrate_tac,
    parse_error_model,
)

logger = logging.getLogger("radot")

LN2 = math.log(2.0)
TRAPEZOID_HEADS = ("zero", "ramp", "const", "extrap")
TRAPEZOID_TAILS = ("phys", "fit")
STP_MODELS = ("physical", "hanscheid2018", "madsen2018")
TEFF_MODELS = ("hanscheid2018", "madsen2018")
ALL_MODELS = tuple(TAC_MODELS) + ("auto", "trapezoid") + STP_MODELS


def model_kind(model: Optional[str]) -> str:
    """'fit' (f1 ... f6a, auto), 'trapezoid' or 'stp'; an unknown or missing model is an error."""
    key = (model or "").strip().lower()
    if key in TAC_MODELS or key == "auto":
        return "fit"
    if key == "trapezoid":
        return "trapezoid"
    if key in STP_MODELS:
        return "stp"
    raise ValueError(f"Unknown model {model!r}. Available: {', '.join(ALL_MODELS)}.")


def check_integration_options(
    model: str,
    n_points: int,
    error_model: Optional[str] = None,
    candidates: Optional[Sequence[str]] = None,
    select: Optional[str] = None,
    head: Optional[str] = None,
    tail: Optional[str] = None,
    tail_points: Optional[int] = None,
    teff_given: bool = False,
) -> str:
    """The options against the model and the number of time points: the same for every ROI, checked
    once. A missing required option or an option the model does not use is an error. Returns the kind."""
    kind = model_kind(model)
    key = model.strip().lower()
    unused = []
    if kind != "fit":
        for name, v in (("error_model", error_model), ("candidates", candidates), ("select", select)):
            if v:
                unused.append(name)
    if kind != "trapezoid":
        for name, v in (("head", head), ("tail", tail), ("tail_points", tail_points)):
            if v is not None:
                unused.append(name)
    if teff_given and key not in TEFF_MODELS:
        unused.append(f"teff (only {', '.join(TEFF_MODELS)})")
    if unused:
        raise ValueError(f"Not used with model '{key}': {'; '.join(unused)}. Remove it.")

    if kind == "stp":
        if n_points != 1:
            raise ValueError(f"Model '{key}' is a single time point model: give one image (got {n_points}).")
        return kind
    if n_points < 2:
        raise ValueError(
            f"Model '{key}' needs at least 2 time points (got {n_points}); with one image use a single time "
            f"point model: {', '.join(STP_MODELS)}."
        )
    if kind == "fit":
        check_model_choice(key, parse_error_model(error_model), n_points, candidates, select)
        return kind
    # trapezoid: head and tail without default
    if head not in TRAPEZOID_HEADS:
        raise ValueError(f"Model 'trapezoid' needs head: {', '.join(TRAPEZOID_HEADS)} (got {head!r}).")
    if tail not in TRAPEZOID_TAILS:
        raise ValueError(f"Model 'trapezoid' needs tail: {', '.join(TRAPEZOID_TAILS)} (got {tail!r}).")
    if tail == "fit":
        if tail_points is None:
            raise ValueError("Tail 'fit' needs tail_points: the number of last points of the fit (>= 2).")
        if not 2 <= int(tail_points) <= n_points:
            raise ValueError(f"tail_points {tail_points}: between 2 and the number of time points ({n_points}).")
    elif tail_points is not None:
        raise ValueError("tail_points is only used with tail 'fit'. Remove it.")
    return kind


def integrate_trapezoid(times_h, values_mbq, lambda_phys_h: float, head: str, tail: str,
                        tail_points: Optional[int] = None) -> dict:
    """Trapezoids between the points, a head 0 -> t1 and a tail tN -> infinity (Kurkowska 2025)."""
    t = np.asarray(times_h, dtype=float)
    a = np.asarray(values_mbq, dtype=float)
    if len(t) < 2:
        raise ValueError(f"Trapezoid needs at least 2 time points, got {len(t)}.")
    t1, a1 = t[0], a[0]
    if head == "zero":
        head_tia = 0.0
    elif head == "ramp":
        head_tia = a1 * t1 / 2.0
    elif head == "const":
        head_tia = a1 * t1
    elif head == "extrap":
        if a[0] <= 0 or a[1] <= 0:
            raise ValueError(f"Head 'extrap' needs the first two values > 0 (got {a[0]:.4g}, {a[1]:.4g} MBq).")
        k = math.log(a[0] / a[1]) / (t[1] - t[0])  # mono-exponential through the first two points
        head_tia = a1 * t1 if abs(k * t1) < 1e-12 else a1 * math.expm1(k * t1) / k
    else:
        raise ValueError(f"Unknown head {head!r}: {', '.join(TRAPEZOID_HEADS)}.")
    body_tia = float(np.sum((a[1:] + a[:-1]) * np.diff(t) / 2.0))
    a_n = a[-1]
    if a_n <= 0:
        raise ValueError(f"The last value is {a_n:.4g} MBq: no tail can start from it.")
    if tail == "phys":
        lam_tail = lambda_phys_h
    elif tail == "fit":
        m = int(tail_points)
        tt, aa = t[-m:], a[-m:]
        if np.any(aa <= 0):
            raise ValueError(f"Tail 'fit' needs the last {m} values > 0.")
        slope = np.polyfit(tt, np.log(aa), 1)[0]  # log-linear least squares
        lam_tail = -float(slope)
        if lam_tail < lambda_phys_h:
            teff = f"{LN2 / lam_tail:.4g} h" if lam_tail > 0 else "no decrease"
            raise ValueError(
                f"Tail fit on the last {m} points slower than the physical decay (Teff {teff} > T_phys "
                f"{LN2 / lambda_phys_h:.4g} h): use tail 'phys'."
            )
    else:
        raise ValueError(f"Unknown tail {tail!r}: {', '.join(TRAPEZOID_TAILS)}.")
    tail_tia = a_n / lam_tail  # from the measured last value, with the tail rate
    tia = head_tia + body_tia + tail_tia
    return {
        "tia_mbq_h": float(tia),
        "tia_se_mbq_h": None,
        "model": "trapezoid",
        "fits": {"trapezoid": {
            "head": head, "tail": tail, "tail_points": tail_points,
            "parts_mbq_h": {"head": float(head_tia), "body": body_tia, "tail": float(tail_tia)},
            "tail_lambda_eff_h": float(lam_tail),
            "teff_terminal_h": float(LN2 / lam_tail),
            "issues": [],
        }},
        "excluded": {},
    }


def integrate_stp(time_h: float, value_mbq: float, lambda_phys_h: float, model: str,
                  teff_h: Optional[float] = None, teff_source: Optional[str] = None) -> dict:
    """Single time point: the curve shape is assumed. physical: A e^(λp t)/λp; hanscheid2018:
    A Teff / ln2; madsen2018: A e^(λeff t)/λeff, λeff = ln2/Teff."""
    key = model.strip().lower()
    t, a = float(time_h), float(value_mbq)
    if key == "physical":
        teff_h, teff_source = LN2 / lambda_phys_h, "physical"
        tia = a * math.exp(lambda_phys_h * t) / lambda_phys_h
    elif key == "hanscheid2018":
        tia = a * teff_h / LN2
    elif key == "madsen2018":
        lam = LN2 / teff_h
        tia = a * math.exp(lam * t) / lam
    else:
        raise ValueError(f"Unknown single time point model {model!r}: {', '.join(STP_MODELS)}.")
    return {
        "tia_mbq_h": float(tia),
        "tia_se_mbq_h": None,
        "model": key,
        "fits": {key: {"teff_h": float(teff_h), "teff_source": teff_source, "time_h": t,
                       "teff_terminal_h": float(teff_h), "issues": []}},
        "excluded": {},
    }


def integrate_curve(times_h, values_mbq, lambda_phys_h: float, model: str, *, error_model=None,
                    candidates=None, select=None, head=None, tail=None, tail_points=None,
                    teff_h: Optional[float] = None, teff_source: Optional[str] = None) -> dict:
    """TIA of one curve with any model (dispatch to integrate_tac, integrate_trapezoid, integrate_stp).
    The options are checked first (check_integration_options)."""
    kind = check_integration_options(model, len(times_h), error_model, candidates, select, head, tail,
                                     tail_points, teff_given=teff_source == "given")
    if kind == "fit":
        return integrate_tac(times_h, values_mbq, lambda_phys_h, model, error_model,
                             candidates=candidates, select=select)
    if kind == "trapezoid":
        return integrate_trapezoid(times_h, values_mbq, lambda_phys_h, head, tail, tail_points)
    key = model.strip().lower()
    if key in TEFF_MODELS and teff_h is None:
        raise ValueError(f"Model '{key}' needs the effective half-life (teff_h).")
    return integrate_stp(times_h[0], values_mbq[0], lambda_phys_h, key, teff_h, teff_source)
