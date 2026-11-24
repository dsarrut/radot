"""Time-activity curve models, fit, model selection (docs/6e_time_integration.md §3–4).

Notation of Kletting et al. 2013 (NUKFIT): sums of exponentials with the physical decay inside the
model, fitted to the activities at the acquisition times (never decay-corrected before the fit):

    f(t) = sum_k A_k exp(-(lambda_k + lambda_phys) t),   TIA = sum_k A_k / (lambda_k + lambda_phys)

Units: t in h, A in MBq, lambda in 1/h, TIA in MBq.h.

Fit: starting points from a log grid of the rates (lambda_k >= 0), the amplitudes being linear for given
rates (weighted least squares, closed form); then all the parameters refined together with
scipy.optimize.least_squares (bounded, analytic Jacobian) from the 3 best grid points.

Error model (variance of point i, data-based weighting, Kletting 2013 Eq. 2-3):
    sigma_i^2 = nu * theta_i^2
    'constant'      theta_i = 1,          nu estimated  (relative weighting, J = K + 1)
    'proportional'  theta_i = y_i,        nu estimated  (relative weighting, J = K + 1)
    'fsd:X'         theta_i = X * y_i,    nu = 1        (absolute weighting, J = K)
    'sd:X'          theta_i = X (MBq),    nu = 1        (absolute weighting, J = K)

Number of points (Kletting 2013, Table I): fit J <= N; parameter SE J + 1 <= N; AICc J + 2 <= N.
"""

import logging
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

logger = logging.getLogger("radot")

LN2 = math.log(2.0)
CV_ACCEPTABLE = 0.5  # parameter CV below 50 %: acceptable (Kletting 2013 Table I, Ivashchenko 2024)
CORRELATION_LIMIT = 0.8  # |correlation| below 0.8 (Kletting 2013 Table I)
RATE_BOUND_FRACTION = 1e-6  # a rate below 1e-6 lambda_phys is at the bound 0
WEIGHT_MIN = 0.01  # models with an Akaike weight above 1 % are used for inference (Kletting 2013)
SELECT_METHODS = ("best", "average")


@dataclass(frozen=True)
class TacModel:
    """A sum of exponentials. Each term is (coefficients over the amplitudes, rate index or None for
    the physical decay only)."""

    name: str
    n_amplitudes: int
    n_rates: int
    terms: Tuple[Tuple[Tuple[float, ...], Optional[int]], ...]
    description: str
    # rates that can be exchanged without changing the model (sorted in the grid search)
    exchangeable: Tuple[int, ...] = ()

    @property
    def n_params(self) -> int:
        return self.n_amplitudes + self.n_rates

    @property
    def param_names(self) -> List[str]:
        return [f"A{i + 1}" for i in range(self.n_amplitudes)] + [
            f"lambda{i + 1}" for i in range(self.n_rates)
        ]


TAC_MODELS: Dict[str, TacModel] = {
    m.name: m
    for m in (
        TacModel("f1", 1, 0, (((1.0,), None),), "A1 e^(-λp t): physical decay only"),
        TacModel("f2", 1, 1, (((1.0,), 0),), "A1 e^(-(λ1+λp) t): mono-exponential"),
        TacModel("f3", 2, 1, (((1.0, 0.0), 0), ((0.0, 1.0), None)),
                 "A1 e^(-(λ1+λp) t) + A2 e^(-λp t): clearance + retained fraction"),
        TacModel("f4", 2, 2, (((1.0, 0.0), 0), ((0.0, 1.0), 1)),
                 "A1 e^(-(λ1+λp) t) + A2 e^(-(λ2+λp) t): bi-exponential", (0, 1)),
        TacModel("f4a", 1, 2, (((1.0,), 0), ((-1.0,), 1)),
                 "A1 [e^(-(λ1+λp) t) - e^(-(λ2+λp) t)]: uptake then washout, f(0) = 0"),
        TacModel("f5", 3, 2, (((1.0, 0.0, 0.0), 0), ((0.0, 1.0, 0.0), 1), ((0.0, 0.0, 1.0), None)),
                 "f4 + A3 e^(-λp t): bi-exponential + retained fraction", (0, 1)),
        TacModel("f6", 3, 3, (((1.0, 0.0, 0.0), 0), ((0.0, 1.0, 0.0), 1), ((0.0, 0.0, 1.0), 2)),
                 "sum of 3 exponentials: tri-exponential", (0, 1, 2)),
        TacModel("f6a", 2, 3, (((1.0, 0.0), 0), ((0.0, 1.0), 1), ((-1.0, -1.0), 2)),
                 "f6 with A3 = -(A1 + A2): tri-exponential, f(0) = 0", (0, 1)),
    )
}


# ----------------------------------------------------------------------------------------------
# error model
# ----------------------------------------------------------------------------------------------


@dataclass(frozen=True)
class ErrorModel:
    name: str  # as given: constant, proportional, fsd:0.1, sd:2
    kind: str  # constant, proportional
    scale: Optional[float]  # None: nu estimated (relative weighting); else known (absolute weighting)

    @property
    def relative(self) -> bool:
        return self.scale is None

    def n_estimated(self, k: int) -> int:
        """J: the adjustable parameters, plus the variance scale nu for relative weighting."""
        return k + 1 if self.relative else k

    def theta(self, y: np.ndarray) -> np.ndarray:
        if self.kind == "constant":
            base = np.ones_like(y)
        else:
            if np.any(y <= 0):
                raise ValueError(
                    f"Error model '{self.name}': the variance is proportional to the measured value, "
                    f"which must be > 0 at every time point (got {', '.join(f'{v:.4g}' for v in y)})."
                )
            base = np.asarray(y, dtype=float).copy()
        return base * (self.scale if self.scale is not None else 1.0)


def parse_error_model(value: Optional[str]) -> ErrorModel:
    """'constant', 'proportional' (nu estimated), 'fsd:X' (known fractional SD), 'sd:X' (known SD in
    MBq). No default: it changes the weights of the points and the model selection."""
    if not value:
        raise ValueError(
            "Missing error model: give constant, proportional (variance scale estimated from the data), "
            "fsd:X (known fractional standard deviation, e.g. fsd:0.1) or sd:X (known standard deviation, "
            "MBq)."
        )
    v = value.strip().lower()
    if v in ("constant", "proportional"):
        return ErrorModel(v, v, None)
    if ":" in v:
        kind, _, num = v.partition(":")
        if kind in ("fsd", "sd"):
            try:
                x = float(num)
            except ValueError:
                x = float("nan")
            if not (x > 0 and math.isfinite(x)):
                raise ValueError(f"Invalid error model '{value}': {kind}:X needs a number X > 0.")
            return ErrorModel(v, "proportional" if kind == "fsd" else "constant", x)
    raise ValueError(
        f"Unknown error model '{value}'. Use constant, proportional, fsd:X or sd:X."
    )


# ----------------------------------------------------------------------------------------------
# model evaluation
# ----------------------------------------------------------------------------------------------


def get_tac_model(name: str) -> TacModel:
    key = name.strip().lower()
    if key not in TAC_MODELS:
        raise ValueError(f"Unknown model '{name}'. Available: {', '.join(TAC_MODELS)}.")
    return TAC_MODELS[key]


def _design(model: TacModel, t: np.ndarray, rates: np.ndarray, lam_p: float) -> np.ndarray:
    """X[i, j]: f(t_i) = X @ amplitudes for the given rates."""
    X = np.zeros((len(t), model.n_amplitudes))
    for coefs, r in model.terms:
        lam = lam_p + (rates[r] if r is not None else 0.0)
        X += np.outer(np.exp(-lam * t), np.asarray(coefs))
    return X


def evaluate(model: TacModel, params: np.ndarray, t, lam_p: float) -> np.ndarray:
    t = np.asarray(t, dtype=float)
    amps, rates = params[: model.n_amplitudes], params[model.n_amplitudes:]
    return _design(model, t, rates, lam_p) @ amps


def _jacobian(model: TacModel, params: np.ndarray, t: np.ndarray, lam_p: float) -> np.ndarray:
    amps, rates = params[: model.n_amplitudes], params[model.n_amplitudes:]
    J = np.zeros((len(t), model.n_params))
    J[:, : model.n_amplitudes] = _design(model, t, rates, lam_p)
    for coefs, r in model.terms:
        if r is None:
            continue
        a = float(np.dot(coefs, amps))
        J[:, model.n_amplitudes + r] += -a * t * np.exp(-(lam_p + rates[r]) * t)
    return J


def tia_and_gradient(model: TacModel, params: np.ndarray, lam_p: float) -> Tuple[float, np.ndarray]:
    """Integral from 0 to infinity (MBq.h) and its gradient with respect to the parameters."""
    amps, rates = params[: model.n_amplitudes], params[model.n_amplitudes:]
    tia = 0.0
    grad = np.zeros(model.n_params)
    for coefs, r in model.terms:
        lam = lam_p + (rates[r] if r is not None else 0.0)
        a = float(np.dot(coefs, amps))
        tia += a / lam
        grad[: model.n_amplitudes] += np.asarray(coefs) / lam
        if r is not None:
            grad[model.n_amplitudes + r] += -a / lam**2
    return tia, grad


# ----------------------------------------------------------------------------------------------
# fit
# ----------------------------------------------------------------------------------------------


def _amplitudes(model, t, y, w, rates, lam_p):
    X = _design(model, t, rates, lam_p)
    sw = np.sqrt(w)
    amps, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
    r = y - X @ amps
    return amps, float(np.sum(w * r * r))


def _grid_wrss(model: TacModel, t, y, w, mesh: np.ndarray, lam_p: float) -> np.ndarray:
    """Weighted RSS of every row of rates in mesh (M, n_rates), the amplitudes in closed form (batched
    weighted least squares through the pseudo-inverse)."""
    X = np.zeros((len(mesh), len(t), model.n_amplitudes))
    for coefs, r in model.terms:
        lam = lam_p + (mesh[:, r] if r is not None else np.zeros(len(mesh)))
        X += np.exp(-lam[:, None] * t[None, :])[:, :, None] * np.asarray(coefs)[None, None, :]
    sw = np.sqrt(w)
    Xw = X * sw[None, :, None]
    yw = y * sw
    amps = np.linalg.pinv(Xw) @ yw  # (M, n_amplitudes)
    r = yw[None, :] - np.einsum("mnk,mk->mn", Xw, amps)
    return np.sum(r * r, axis=1)


def _rate_grid(n_rates: int) -> np.ndarray:
    """Biological rates tried as starting points (1/h): 0, then half-lives from ~7000 h to ~4 min."""
    n = {1: 60, 2: 30, 3: 14}[n_rates]
    return np.concatenate([[0.0], np.logspace(-4, 1, n)])


def _fit_rates(model: TacModel, t, y, w, lam_p) -> Tuple[np.ndarray, np.ndarray, float]:
    """Rates (>= 0) minimising the weighted RSS, with the amplitudes in closed form."""
    if model.n_rates == 0:
        amps, wrss = _amplitudes(model, t, y, w, np.zeros(0), lam_p)
        return amps, np.zeros(0), wrss
    grid = _rate_grid(model.n_rates)
    mesh = np.array(np.meshgrid(*([grid] * model.n_rates), indexing="ij")).reshape(model.n_rates, -1).T
    if len(model.exchangeable) > 1:  # keep one ordering of the exchangeable rates
        ex = list(model.exchangeable)
        keep = np.all(np.diff(mesh[:, ex], axis=1) > 0, axis=1)
        mesh = mesh[keep]
    if model.name == "f4a":  # uptake (lambda2) faster than washout (lambda1), the other order is -A1
        mesh = mesh[mesh[:, 1] > mesh[:, 0]]
    scores = _grid_wrss(model, t, y, w, mesh, lam_p)
    starts = mesh[np.argsort(scores)[:3]]

    # refinement of all the parameters together: bounded least squares (trust region reflective,
    # rates >= 0), analytic Jacobian, from the best grid points
    from scipy.optimize import least_squares

    sw = np.sqrt(w)
    lower = np.r_[np.full(model.n_amplitudes, -np.inf), np.zeros(model.n_rates)]

    def residuals(p):
        return sw * (evaluate(model, p, t, lam_p) - y)

    def jac(p):
        return _jacobian(model, p, t, lam_p) * sw[:, None]

    best = None
    for rates0 in starts:
        amps0, _ = _amplitudes(model, t, y, w, rates0, lam_p)
        p0 = np.r_[amps0, rates0]
        res = least_squares(residuals, p0, jac=jac, bounds=(lower, np.inf), method="trf",
                            x_scale="jac", xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=2000)
        if best is None or res.cost < best.cost:
            best = res
    params = best.x
    amps, rates = params[: model.n_amplitudes], params[model.n_amplitudes:].copy()
    # a rate at the bound (active constraint, or a biological half-life > 1e6 physical half-lives) is
    # set to exactly 0, the amplitudes solved again
    active = best.active_mask[model.n_amplitudes:] != 0
    at_bound = active | (rates <= RATE_BOUND_FRACTION * lam_p)
    if np.any(at_bound):
        rates[at_bound] = 0.0
        amps, _ = _amplitudes(model, t, y, w, rates, lam_p)
    r = y - evaluate(model, np.r_[amps, rates], t, lam_p)
    return amps, rates, float(np.sum(w * r * r))


@dataclass
class FitResult:
    model: str
    params: Dict[str, float]
    tia_mbq_h: float
    n_points: int
    n_estimated: int  # J
    wrss: float
    nu: float
    minus2_log_likelihood: float
    r2: float
    params_se: Optional[Dict[str, float]] = None
    params_cv: Optional[Dict[str, float]] = None
    correlation_max: Optional[float] = None
    tia_se_mbq_h: Optional[float] = None
    aicc: Optional[float] = None
    weight: Optional[float] = None
    at_bound: List[str] = field(default_factory=list)
    teff_terminal_h: Optional[float] = None
    weighted_residuals: List[float] = field(default_factory=list)
    issues: List[str] = field(default_factory=list)  # quality criteria not met (Kletting 2013 Table I)

    def as_dict(self) -> dict:
        d = {
            "params": self.params,
            "params_se": self.params_se,
            "params_cv": self.params_cv,
            "correlation_max": self.correlation_max,
            "tia_mbq_h": self.tia_mbq_h,
            "tia_se_mbq_h": self.tia_se_mbq_h,
            "aicc": self.aicc,
            "weight": self.weight,
            "r2": self.r2,
            "nu": self.nu,
            "n_points": self.n_points,
            "n_estimated": self.n_estimated,
            "at_bound": self.at_bound,
            "teff_terminal_h": self.teff_terminal_h,
            "weighted_residuals": self.weighted_residuals,
            "issues": self.issues,
        }
        return d


def fit_tac(model_name: str, times_h, values_mbq, lambda_phys_h: float, error_model: ErrorModel) -> FitResult:
    """Fit one model. An unsuitable number of points or a non-physical curve is a ValueError; the
    quality criteria (CV, correlation, rates at the bound) are reported in `issues`, not raised."""
    model = get_tac_model(model_name)
    t = np.asarray(times_h, dtype=float)
    y = np.asarray(values_mbq, dtype=float)
    n = len(t)
    k = model.n_params
    j = error_model.n_estimated(k)
    if j > n:
        raise ValueError(
            f"Model {model.name}: {k} parameter(s){' + the variance scale' if error_model.relative else ''} "
            f"(J = {j}) for {n} time points: J <= N needed (Kletting 2013, Table I)."
        )
    theta = error_model.theta(y)
    w = 1.0 / theta**2
    amps, rates, wrss = _fit_rates(model, t, y, w, lambda_phys_h)
    params = np.concatenate([amps, rates])
    names = model.param_names

    tia, grad = tia_and_gradient(model, params, lambda_phys_h)
    # non-physical: the curve must stay >= 0 after injection, and the integral > 0
    t_check = np.concatenate([t, np.linspace(0.0, 10.0 * t.max(), 400)])
    f_check = evaluate(model, params, t_check, lambda_phys_h)
    scale = max(float(np.max(np.abs(y))), 1e-300)
    if not np.all(np.isfinite(params)) or not math.isfinite(tia):
        raise ValueError(f"Model {model.name}: the fit did not converge (non-finite parameters).")
    if tia <= 0 or np.min(f_check) < -1e-6 * scale:
        raise ValueError(
            f"Model {model.name}: non-physical fit (the curve is negative after injection, "
            f"TIA {tia:.4g} MBq.h): " + ", ".join(f"{nm}={v:.4g}" for nm, v in zip(names, params)) + "."
        )

    nu_ml = wrss / n if error_model.relative else 1.0
    nu_ml = max(nu_ml, 1e-300)
    m2ll = float(np.sum(np.log(2.0 * math.pi * theta**2)) + n * math.log(nu_ml) + wrss / nu_ml)
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    rss = float(np.sum((y - evaluate(model, params, t, lambda_phys_h)) ** 2))
    r2 = 1.0 - rss / ss_tot if ss_tot > 0 else float("nan")
    resid_w = (y - evaluate(model, params, t, lambda_phys_h)) / (theta * math.sqrt(nu_ml))

    res = FitResult(
        model=model.name,
        params={nm: float(v) for nm, v in zip(names, params)},
        tia_mbq_h=float(tia),
        n_points=n,
        n_estimated=j,
        wrss=wrss,
        nu=float(nu_ml),
        minus2_log_likelihood=m2ll,
        r2=float(r2),
        weighted_residuals=[float(v) for v in resid_w],
    )
    res.at_bound = [f"lambda{i + 1}" for i, r in enumerate(rates) if r == 0.0]
    if res.at_bound:
        res.issues.append(f"{', '.join(res.at_bound)} at the bound 0 (no biological clearance)")
    # terminal effective half-life: slowest term with a non-zero amplitude
    lams = [lambda_phys_h + (rates[r] if r is not None else 0.0) for coefs, r in model.terms
            if abs(float(np.dot(coefs, amps))) > 0]
    res.teff_terminal_h = float(LN2 / min(lams)) if lams else None

    # parameter SE: J + 1 <= N (Kletting 2013, Table I); covariance nu (J^T W J)^-1, nu = WRSS/(N-K)
    # for relative weighting (Ivashchenko 2024 Eq. 8), 1 for absolute weighting
    if j + 1 <= n:
        jac = _jacobian(model, params, t, lambda_phys_h)
        fim = jac.T @ (jac * w[:, None])
        nu_cov = wrss / (n - k) if error_model.relative else 1.0
        try:
            cov = nu_cov * np.linalg.inv(fim)
            se = np.sqrt(np.clip(np.diag(cov), 0, None))
            ok = np.all(np.isfinite(cov)) and np.linalg.cond(fim) < 1e14
        except np.linalg.LinAlgError:
            ok = False
        if ok:
            res.params_se = {nm: float(s) for nm, s in zip(names, se)}
            res.params_cv = {nm: (float(s / abs(p)) if p != 0 else float("inf"))
                             for nm, s, p in zip(names, se, params)}
            res.tia_se_mbq_h = float(math.sqrt(max(float(grad @ cov @ grad), 0.0)))
            d = np.sqrt(np.clip(np.diag(cov), 0, None))
            if k > 1 and np.all(d > 0):  # a parameter with a zero SE (degenerate fit): no correlation
                corr = cov / np.outer(d, d)
                off = np.abs(corr[~np.eye(k, dtype=bool)])
                res.correlation_max = float(np.max(off))
                if res.correlation_max >= CORRELATION_LIMIT:
                    res.issues.append(f"parameter correlation {res.correlation_max:.2f} >= {CORRELATION_LIMIT}")
            bad = [f"{nm} {cv:.0%}" for nm, cv in res.params_cv.items() if not cv < CV_ACCEPTABLE]
            if bad:
                res.issues.append(f"parameter CV >= {CV_ACCEPTABLE:.0%}: {', '.join(bad)}")
        else:
            res.issues.append("parameter covariance not computable (singular)")
    else:
        res.issues.append(f"parameter SE not computable: J + 1 = {j + 1} > N = {n}")

    if j + 2 <= n:
        res.aicc = m2ll + 2 * j + 2 * j * (j + 1) / (n - j - 1)
    return res


def passes_quality(res: FitResult) -> bool:
    """Used for model selection (Kletting 2013: only functions passing the criteria): parameter CV
    computable and below 50 %. Correlation and rates at the bound are reported, not excluded."""
    return res.params_cv is not None and all(cv < CV_ACCEPTABLE for cv in res.params_cv.values())


def akaike_weights(aiccs: Sequence[float]) -> np.ndarray:
    a = np.asarray(aiccs, dtype=float)
    d = a - np.min(a)
    e = np.exp(-d / 2.0)
    return e / np.sum(e)


def check_model_choice(model: str, em: ErrorModel, n: int, candidates=None, select=None):
    """The model choice against the number of time points N (the same for every ROI: checked once).
    Returns (model key, candidate names or None)."""
    key = get_tac_model(model).name if model.strip().lower() != "auto" else "auto"
    rel = " + 1 for the variance scale" if em.relative else ""
    if key != "auto":
        if candidates or select:
            raise ValueError(
                f"candidates / select are only used with model 'auto' (model '{model}' given)."
            )
        j = em.n_estimated(TAC_MODELS[key].n_params)
        if j > n:
            raise ValueError(
                f"Model {key}: J = {j} (parameters{rel}) for N = {n} time points: J <= N needed "
                f"(Kletting 2013, Table I)."
            )
        return key, None
    if not candidates:
        raise ValueError("Model 'auto' needs the candidate models (e.g. --candidates f2,f3,f4a).")
    if select not in SELECT_METHODS:
        raise ValueError(f"Model 'auto' needs select: {' or '.join(SELECT_METHODS)} (got {select!r}).")
    names = [get_tac_model(c).name for c in candidates]
    if len(set(names)) != len(names) or len(names) < 2:
        raise ValueError(f"Model 'auto' needs at least 2 different candidate models (got {list(candidates)}).")
    too_many = [c for c in names if em.n_estimated(TAC_MODELS[c].n_params) + 2 > n]
    if too_many:
        raise ValueError(
            f"Candidate model(s) {', '.join(too_many)}: the AICc needs J + 2 <= N (N = {n} time points, "
            f"J = parameters{rel}; Kletting 2013, Table I). Remove them from the candidates."
        )
    return key, names


def integrate_tac(
    times_h,
    values_mbq,
    lambda_phys_h: float,
    model: str,
    error_model,
    candidates: Optional[Sequence[str]] = None,
    select: Optional[str] = None,
) -> dict:
    """TIA of one time-activity curve.

    model: a model name (f1 ... f6a), or 'auto' with candidates and select ('best': lowest AICc;
    'average': TIA averaged over the models with an Akaike weight > 1 %, Kletting 2013 Eq. 10-11).
    Returns {"tia_mbq_h", "tia_se_mbq_h", "model", "fits": {name: FitResult.as_dict()}, "excluded"}.
    """
    em = error_model if isinstance(error_model, ErrorModel) else parse_error_model(error_model)
    t = np.asarray(times_h, dtype=float)
    if len(t) < 2:
        raise ValueError(f"Time integration needs at least 2 time points, got {len(t)}.")
    if np.any(np.diff(np.sort(t)) <= 0):
        raise ValueError(f"Two time points have the same time from injection: {list(times_h)} h.")
    if np.any(t < 0):
        raise ValueError(f"Negative time from injection: {list(times_h)} h.")
    key, names = check_model_choice(model, em, len(t), candidates, select)
    if key != "auto":
        res = fit_tac(key, t, values_mbq, lambda_phys_h, em)
        return {
            "tia_mbq_h": res.tia_mbq_h,
            "tia_se_mbq_h": res.tia_se_mbq_h,
            "model": res.model,
            "fits": {res.model: res.as_dict()},
            "excluded": {},
        }

    fits, excluded = {}, {}
    for c in names:
        try:
            res = fit_tac(c, t, values_mbq, lambda_phys_h, em)
        except ValueError as e:
            excluded[c] = str(e)
            continue
        if not passes_quality(res):
            excluded[c] = "; ".join(res.issues) or "quality criteria not met"
        fits[c] = res
    kept = [c for c in fits if c not in excluded]
    if not kept:
        raise ValueError(
            "No candidate model passes the quality criteria: "
            + "; ".join(f"{c}: {r}" for c, r in excluded.items())
        )
    w = akaike_weights([fits[c].aicc for c in kept])
    for c, wi in zip(kept, w):
        fits[c].weight = float(wi)
    if select == "best":
        best = kept[int(np.argmin([fits[c].aicc for c in kept]))]
        tia, se, label = fits[best].tia_mbq_h, fits[best].tia_se_mbq_h, best
    else:
        used = [c for c in kept if fits[c].weight > WEIGHT_MIN]
        wu = np.array([fits[c].weight for c in used])
        wu = wu / wu.sum()
        tias = np.array([fits[c].tia_mbq_h for c in used])
        ses = np.array([fits[c].tia_se_mbq_h for c in used])
        tia = float(np.sum(wu * tias))
        se = float(np.sum(wu * np.sqrt(ses**2 + (tias - tia) ** 2)))
        label = "average"
    return {
        "tia_mbq_h": float(tia),
        "tia_se_mbq_h": se,
        "model": label,
        "fits": {c: fits[c].as_dict() for c in fits},
        "excluded": excluded,
    }
