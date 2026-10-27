import subprocess
import sys
from pathlib import Path

import pytest

from benchmarks import CASES, Case, evaluate

CHECKED_CASES = [c for c in CASES if c.accept is not None]


@pytest.mark.parametrize("case", CHECKED_CASES, ids=str)
def test_radot_vs_opengate(case):
    row = evaluate(case)
    lo, hi = case.accept
    assert row["radot"] > 0.0 and row["gate"] > 0.0
    assert row["passed"], (
        f"{case}: radot {row['radot']:.4g} vs GATE {row['gate']:.4g} mGy/MBq, "
        f"difference {row['diff']:+.2f} % (expected |diff| in [{lo}, {hi}) %)"
    )


def test_ac225_parent_only_and_cascade():
    # Parent-only Madsen and LDM agree (both ignore the daughters) ...
    madsen = evaluate(Case("Ac225", "madsen2018", None))
    ldm = evaluate(Case("Ac225", "ldm", None))
    assert abs(madsen["radot"] - ldm["radot"]) / ldm["radot"] < 0.01

    # ... and GATE (full decay chain) is ~4.7 times higher than the parent alone
    ratio = ldm["gate"] / ldm["radot"]
    assert 4.5 < ratio < 5.0, f"GATE / LDM parent-only ratio: {ratio:.3f}"


def test_validation_page_is_up_to_date():
    # docs/7_validation.md is generated: python validation/opengate/report.py
    report = Path(__file__).resolve().parent / "report.py"
    result = subprocess.run(
        [sys.executable, str(report), "--check"], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
