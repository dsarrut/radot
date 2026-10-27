import subprocess
import sys
from pathlib import Path

import pytest

from benchmarks import (
    BENCHMARK_FILES,
    CASES,
    Case,
    download_data,
    evaluate,
    gate_from_images,
    load_benchmark,
    organ_masks,
)

CHECKED_CASES = [c for c in CASES if c.accept is not None]


@pytest.fixture(scope="session")
def data_dir(request):
    """None for the default level; the downloaded images for --opengate=images."""
    if request.config.getoption("--opengate") == "json":
        return None
    return download_data()


@pytest.mark.parametrize("case", CHECKED_CASES, ids=str)
def test_radot_vs_opengate(case, data_dir):
    row = evaluate(case, data_dir=data_dir)
    lo, hi = case.accept
    assert row["radot"] > 0.0 and row["gate"] > 0.0
    assert row["passed"], (
        f"{case}: radot {row['radot']:.4g} vs GATE {row['gate']:.4g} mGy/MBq, "
        f"difference {row['diff']:+.2f} % (expected |diff| in [{lo}, {hi}) %)"
    )


def test_ac225_parent_only_and_cascade(data_dir):
    # Parent-only Madsen and LDM agree (both ignore the daughters) ...
    madsen = evaluate(Case("Ac225", "madsen2018", None), data_dir=data_dir)
    ldm = evaluate(Case("Ac225", "ldm", None), data_dir=data_dir)
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


def test_validation_page_is_up_to_date():
    # docs/7_validation.md is generated: python validation/opengate/report.py
    report = Path(__file__).resolve().parent / "report.py"
    result = subprocess.run(
        [sys.executable, str(report), "--check"], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
