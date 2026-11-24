"""Time integration of a series of images (docs/6e_time_integration.md): `radot integrate`.

ROI level, fit first (E1): ROI time-activity curves, sums-of-exponentials models (Kletting 2013),
error model, AICc model selection. Always import from `radot.integrate`.
"""

from radot.integrate.models import (
    CV_ACCEPTABLE,
    SELECT_METHODS,
    TAC_MODELS,
    ErrorModel,
    FitResult,
    akaike_weights,
    evaluate,
    fit_tac,
    get_tac_model,
    integrate_tac,
    parse_error_model,
)
from radot.integrate.simple import (
    ALL_MODELS,
    STP_MODELS,
    TRAPEZOID_HEADS,
    TRAPEZOID_TAILS,
    check_integration_options,
    integrate_curve,
    integrate_stp,
    integrate_trapezoid,
)
from radot.integrate.roi import (
    ACTIVITY_UNITS,
    DOSE_RATE_UNITS,
    INTEGRATE_FORMAT,
    QUANTITIES,
    quantity_of_unit,
    integrate_rois,
    is_integrate_file,
    read_integrate,
    roi_time_activity_curves,
    write_integrate,
)

__all__ = [
    "ALL_MODELS",
    "STP_MODELS",
    "TRAPEZOID_HEADS",
    "TRAPEZOID_TAILS",
    "check_integration_options",
    "integrate_curve",
    "integrate_stp",
    "integrate_trapezoid",
    "CV_ACCEPTABLE",
    "SELECT_METHODS",
    "TAC_MODELS",
    "ErrorModel",
    "FitResult",
    "akaike_weights",
    "evaluate",
    "fit_tac",
    "get_tac_model",
    "integrate_tac",
    "parse_error_model",
    "ACTIVITY_UNITS",
    "DOSE_RATE_UNITS",
    "INTEGRATE_FORMAT",
    "QUANTITIES",
    "quantity_of_unit",
    "integrate_rois",
    "is_integrate_file",
    "read_integrate",
    "roi_time_activity_curves",
    "write_integrate",
]
