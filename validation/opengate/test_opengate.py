import json
import subprocess
import sys
from pathlib import Path

import pytest

from benchmarks import (
    BENCH_DIR,
    BENCHMARK_FILES,
    CASES,
    Case,
    download_data,
    evaluate,
    gate_from_images,
    load_benchmark,
    organ_masks,
    run_opengate,
    SIM_DIR,
)

CHECKED_CASES = [c for c in CASES if c.accept is not None]


@pytest.fixture(scope="session")
def data_dir(request):
    """None for the default level; the downloaded images for --opengate=images."""
    if request.config.getoption("--opengate") == "json":
        return None
    return download_data()


@pytest.fixture(scope="session")
def gate_dir(request, data_dir):
    """None (published OpenGATE outputs), or new simulations for --opengate=simulate."""
    if request.config.getoption("--opengate") != "simulate":
        return None
    for isotope in BENCHMARK_FILES:
        run_opengate(isotope, data_dir, SIM_DIR)
    return SIM_DIR


@pytest.mark.parametrize("case", CHECKED_CASES, ids=str)
def test_radot_vs_opengate(case, data_dir, gate_dir):
    row = evaluate(case, data_dir=data_dir, gate_dir=gate_dir)
    lo, hi = case.accept
    assert row["radot"] > 0.0 and row["gate"] > 0.0
    assert row["passed"], (
        f"{case}: radot {row['radot']:.4g} vs GATE {row['gate']:.4g} mGy/MBq, "
        f"difference {row['diff']:+.2f} % (expected |diff| in [{lo}, {hi}) %)"
    )


def test_ac225_parent_only_and_cascade(data_dir, gate_dir):
    # Parent-only Madsen and LDM agree (both ignore the daughters) ...
    madsen = evaluate(Case("Ac225", "madsen2018", None), data_dir=data_dir, gate_dir=gate_dir)
    ldm = evaluate(Case("Ac225", "ldm", None), data_dir=data_dir, gate_dir=gate_dir)
    assert abs(madsen["radot"] - ldm["radot"]) / ldm["radot"] < 0.01

    # ... and GATE (full decay chain) is ~4.7 times higher than the parent alone
    ratio = ldm["gate"] / ldm["radot"]
    assert 4.5 < ratio < 5.0, f"GATE / LDM parent-only ratio: {ratio:.3f}"


@pytest.mark.parametrize("isotope", list(BENCHMARK_FILES))
def test_benchmark_json_matches_opengate_images(isotope, data_dir):
    # The stored OpenGATE values are those of the published images
    if data_dir is None:
        pytest.skip("needs --opengate=images")
    bench = load_benchmark(isotope)
    source, target = organ_masks(bench, data_dir)
    rate, _ = gate_from_images(isotope, source, data_dir)
    assert rate == pytest.approx(bench["opengate"]["mean_dose_rate_gy_s"], rel=1e-12)
    rate, _ = gate_from_images(isotope, target, data_dir)
    assert rate == pytest.approx(bench["cross_dose"]["opengate"]["mean_dose_rate_gy_s"], rel=1e-12)


@pytest.mark.parametrize("kind", ["self", "cross"])
@pytest.mark.parametrize("isotope", list(BENCHMARK_FILES))
def test_new_simulation_matches_stored(isotope, kind, data_dir, gate_dir):
    # A new Monte Carlo run agrees with the stored one within 3 combined standard uncertainties
    if gate_dir is None:
        pytest.skip("needs --opengate=simulate")
    bench = load_benchmark(isotope)
    source, target = organ_masks(bench, data_dir)
    stored = bench["opengate"] if kind == "self" else bench["cross_dose"]["opengate"]
    if stored["uncertainty_estimate_fraction"] >= 0.5:
        # e.g. Y90 cross-dose: almost no photons, the stored value is noise
        pytest.skip(f"{isotope} {kind}: stored uncertainty {stored['uncertainty_estimate_fraction']:.0%}, not meaningful")
    rate, unc = gate_from_images(isotope, source if kind == "self" else target, data_dir, gate_dir)
    ref = stored["mean_dose_rate_gy_s"]
    sigma = (unc**2 + stored["uncertainty_estimate_fraction"] ** 2) ** 0.5
    diff = abs(rate - ref) / ref
    assert diff < 3 * sigma, f"{isotope} {kind}: new {rate:.4e} vs stored {ref:.4e} Gy/s, {diff:.2%} > 3 x {sigma:.2%}"


@pytest.mark.parametrize("kind", ["self", "cross"])
def test_uncertainty_estimate_matches_batches(kind):
    # The organ uncertainty estimated from the voxel maps (benchmark JSON) agrees with the one
    # measured with 20 independent simulations (batch_uncertainty.py), within 2 sigma of the measurement
    batches = json.loads((BENCH_DIR / "p142_lu177_batch_uncertainty.json").read_text())[kind]
    bench = load_benchmark("Lu177")
    stored = bench["opengate"] if kind == "self" else bench["cross_dose"]["opengate"]
    estimate = stored["uncertainty_estimate_fraction"]
    measured = batches["measured_uncertainty_fraction"]
    precision = batches["measured_uncertainty_precision_fraction"]
    assert abs(estimate - measured) < 2 * precision, (
        f"Lu177 {kind}: estimated {estimate:.3%}, measured {measured:.3%} +/- {precision:.3%}"
    )


def test_validation_page_is_up_to_date():
    # docs/7_validation.md is generated: python validation/opengate/report.py
    report = Path(__file__).resolve().parent / "report.py"
    result = subprocess.run(
        [sys.executable, str(report), "--check"], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
