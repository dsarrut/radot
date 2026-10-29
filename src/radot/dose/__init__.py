"""
Dosimetry: models, ROI dose computation, dose-rate maps, multi-time-point and uncertainty.

The public API is re-exported here: import from radot.dose, not from the sub-modules.
"""

from radot.dose.models import (
    dose_ldm,
    dose_hanscheid2017,
    dose_hanscheid2018,
    dose_madsen2018,
    dose_svalue,
)
from radot.dose.uncertainty import (
    ISOTOPE_MEAN_BETA_RANGE_MM,
    roi_dose_uncertainty_bounds,
    roi_mean_chord_mm,
    roi_dose_uncertainty_estimate,
)
from radot.dose.mtp import (
    fit_exp_linear,
    decay_corrected_tac,
    triexpo_fit_jackson2013,
    triexpo_integrate_jackson2013,
    dose_jackson2013,
    compute_roi_dose_mtp,
)
from radot.dose.strategies import (
    strategy_hanscheid2017,
    strategy_hanscheid2018,
    strategy_madsen2018,
    strategy_svalue,
    DOSIMETRY_STRATEGIES,
    TIA_UNITS,
    TIA_METHODS,
    compute_dose_generic,
)
from radot.dose.rate_map import (
    DOSE_RATE_MAP_METHODS,
    DOSE_RATE_MAP_UNITS,
    dose_from_rate_map,
)
from radot.dose.compute import check_dosimetry_inputs, compute_roi_dose
from radot.dose.config import load_dose_config
from radot.dose.files import compute_roi_dose_file, pretty_print_dose
