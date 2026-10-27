import json
import subprocess
import sys
from pathlib import Path

import pytest

from benchmarks import (
    BENCH_DIR,
    BENCHMARK_FILES,
    CASES,
    CROSS_TARGETS,
    Case,
    cross_dose,
    download_data,
    evaluate,
    gate_from_images,
    load_benchmark,
    organ_mask,
    run_opengate,
    SIM_DIR,
)

CHECKED_CASES = [c for c in CASES if c.accept is not None]
# The source organ ("": self-dose), then the cross-dose targets
ORGANS = [""] + list(CROSS_TARGETS)


def stored_gate(bench: dict, organ: str) -> dict:
    return cross_dose(bench, organ)["opengate"] if organ else bench["opengate"]


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
    for organ in ORGANS:
        rate, _ = gate_from_images(isotope, organ_mask(data_dir, organ), data_dir)
        assert rate == pytest.approx(stored_gate(bench, organ)["mean_dose_rate_gy_s"], rel=1e-12), organ


@pytest.mark.parametrize("organ", ORGANS, ids=lambda o: o or "Liver")
@pytest.mark.parametrize("isotope", list(BENCHMARK_FILES))
def test_new_simulation_matches_stored(isotope, organ, data_dir, gate_dir):
    # A new Monte Carlo run agrees with the stored one within 3 combined standard uncertainties
    if gate_dir is None:
        pytest.skip("needs --opengate=simulate")
    stored = stored_gate(load_benchmark(isotope), organ)
    name = organ or "Liver"
    if stored["uncertainty_estimate_fraction"] >= 0.5:
        # e.g. Y90 cross-doses: almost no photons, the stored value is noise
        pytest.skip(f"{isotope} {name}: stored uncertainty {stored['uncertainty_estimate_fraction']:.0%}, not meaningful")
    rate, unc = gate_from_images(isotope, organ_mask(data_dir, organ), data_dir, gate_dir)
    ref = stored["mean_dose_rate_gy_s"]
    sigma = (unc**2 + stored["uncertainty_estimate_fraction"] ** 2) ** 0.5
    diff = abs(rate - ref) / ref
    assert diff < 3 * sigma, f"{isotope} {name}: new {rate:.4e} vs stored {ref:.4e} Gy/s, {diff:.2%} > 3 x {sigma:.2%}"


@pytest.mark.parametrize("organ", ORGANS, ids=lambda o: o or "Liver")
def test_uncertainty_estimate_matches_batches(organ):
    # The organ uncertainty estimated from the voxel maps (benchmark JSON) agrees with the one
    # measured with 20 independent simulations (batch_uncertainty.py), within 2 sigma of the measurement
    bench = load_benchmark("Lu177")
    batches = json.loads((BENCH_DIR / "p142_lu177_batch_uncertainty.json").read_text())
    measured = batches["organs"][organ or bench["roi_name"]]
    estimate = stored_gate(bench, organ)["uncertainty_estimate_fraction"]
    m, precision = measured["measured_uncertainty_fraction"], measured["measured_uncertainty_precision_fraction"]
    assert abs(estimate - m) < 2 * precision, (
        f"Lu177 {organ or 'Liver'}: estimated {estimate:.3%}, measured {m:.3%} +/- {precision:.3%}"
    )


def test_validation_page_is_up_to_date():
    # docs/7_validation.md is generated: python validation/opengate/report.py
    report = Path(__file__).resolve().parent / "report.py"
    result = subprocess.run(
        [sys.executable, str(report), "--check"], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
