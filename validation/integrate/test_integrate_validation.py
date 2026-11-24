"""
radot integrate on synthetic ¹⁷⁷Lu curves with a known integral (CI, no data).

The cases and the acceptance ranges are in tac_benchmarks.py; the full Monte Carlo and the page are made by
report.py (results/integrate_results.json, docs/9_validation_time_integration.md).
"""

import subprocess
import sys
import warnings
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import tac_benchmarks as B  # noqa: E402

N_CI = 300  # realisations in CI (the page uses report.N_FULL)

# the model contains the curve and the points constrain it: exact without noise
EXACT = [
    ("lesion", "3tp", "f4a fsd"),
    ("lesion", "4tp", "f4a fsd"),
    ("lesion", "4tp", "f4a proportional"),
    ("lesion", "6tp", "f4a fsd"),
    ("lesion", "6tp", "auto best"),
    ("kidney", "6tp", "f4a fsd"),
]


@pytest.mark.parametrize("shape, schedule, method", EXACT)
def test_noise_free_exact(shape, schedule, method):
    err = B.noise_free(shape, schedule, method)
    assert abs(err) < B.NOISE_FREE_EXACT_PCT, err


@pytest.mark.parametrize("noise", B.NOISES)
def test_true_model_unbiased_with_honest_standard_error(noise):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = B.run_config("lesion", "6tp", noise, "f4a fsd", N_CI)
    assert r["failures"] == 0 and r["n_se"] == N_CI
    assert abs(r["bias_pct"]) < B.TRUE_MODEL_BIAS_PCT
    assert B.COVERAGE_68_RANGE[0] <= r["coverage_68"] <= B.COVERAGE_68_RANGE[1]
    assert B.COVERAGE_95_RANGE[0] <= r["coverage_95"] <= B.COVERAGE_95_RANGE[1]


def test_truth_is_the_integral_of_the_curve():
    import numpy as np

    t = np.linspace(0.0, 20000.0, 2_000_001)
    for s in B.SHAPES.values():
        y = s.curve(t)
        assert float(np.sum((y[1:] + y[:-1]) * np.diff(t) / 2)) == pytest.approx(s.tia_mbq_h(), rel=1e-5)


def test_validation_page_is_up_to_date():
    # docs/9_validation_time_integration.md is generated from results/integrate_results.json
    report = Path(__file__).resolve().parent / "report.py"
    result = subprocess.run([sys.executable, str(report), "--check"], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
