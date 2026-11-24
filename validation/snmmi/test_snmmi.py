"""
SNMMI 177Lu Dosimetry Challenge: radot vs the participants' results.

test_challenge_stats runs everywhere (CI). The T5 test needs the Deep Blue data
(downloaded by hand, then prepare_data.py) and is skipped otherwise.
"""

import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cases import PATIENTS, challenge_stats, has_nifti  # noqa: E402

# |radot - participants' median| / median accepted for T5 (the participants' QCD is 3-7 %)
T5_TOLERANCE_PCT = 8.0
# T4, f4a: the participants' spread is larger (QCD up to ~25 %, fit choices); kidney R of A is -12 %
T4_TOLERANCE_PCT = 15.0
T5_STORED = Path(__file__).resolve().parent / "results" / "t5_radot.json"


def test_challenge_stats():
    """The parsed Supplemental Table 1 is consistent with Table 3 of Brosch-Lenz 2023."""
    for task in ("T1", "T4", "T5"):
        for name in PATIENTS:
            for voi, s in challenge_stats(task, name).items():
                assert s["min"] <= s["q1"] <= s["median"] <= s["q3"] <= s["max"], (task, name, voi)
    # Table 3, % difference of the medians, T5 vs T1
    a, a1 = challenge_stats("T5", "patient_4"), challenge_stats("T1", "patient_4")
    b, b1 = challenge_stats("T5", "patient_6"), challenge_stats("T1", "patient_6")
    assert 100 * (a["Liver"]["median"] / a1["Liver"]["median"] - 1) == pytest.approx(38, abs=0.5)
    assert 100 * (a["Spleen"]["median"] / a1["Spleen"]["median"] - 1) == pytest.approx(-12, abs=0.5)
    assert 100 * (b["R Kidney"]["median"] / b1["R Kidney"]["median"] - 1) == pytest.approx(-5, abs=0.5)
    assert "Spleen" not in b  # patient B is splenectomised


@pytest.fixture(scope="module")
def t5_results():
    if not all(has_nifti(name) for name in PATIENTS):
        pytest.skip("SNMMI challenge data not prepared (see validation/snmmi/README.md)")
    from t5 import compare, run_patient

    return compare({name: run_patient(name) for name in PATIENTS})


@pytest.mark.parametrize("name", list(PATIENTS))
def test_t5(t5_results, name):
    for voi, r in t5_results[name].items():
        if voi.startswith("Liver") and PATIENTS[name].challenge_id == "A":
            continue  # liver definition differs between participants, see test_t5_liver_a
        if voi == "Liver (no lesions)":
            continue
        for method in ("svalue", "ldm", "gate"):
            if r.get(method) is not None:
                assert abs(r[f"{method}_diff_pct"]) < T5_TOLERANCE_PCT, (name, voi, method, r)


def test_t5_liver_a(t5_results):
    """Patient A: the liver VOI contains lesion 1 (116 mL). The participants' spread (QCD 32 %)
    matches the two definitions: Q1 ~ liver without the lesions, median ~ liver VOI as given."""
    r = t5_results["patient_4"]
    c = r["Liver"]["challenge"]
    assert r["Liver"]["svalue"] == pytest.approx(c["median"], rel=0.08)
    assert r["Liver (no lesions)"]["svalue"] == pytest.approx(c["q1"], rel=0.08)


@pytest.fixture(scope="module")
def t4_results():
    from cases import has_t4_nifti

    if not all(has_t4_nifti(name) for name in PATIENTS) or not T5_STORED.exists():
        pytest.skip("SNMMI challenge T4 data not prepared (see validation/snmmi/README.md)")
    import json

    from t4 import compare, run_patient

    t5 = json.loads(T5_STORED.read_text())
    return compare({name: run_patient(name, t5[name]) for name in PATIENTS})


def test_t4_reproduces_the_stored_results(t4_results):
    """results/t4_radot.json (the page) is what t4.py computes now."""
    import json

    from t4 import METHODS

    stored = json.loads((T5_STORED.parent / "t4_radot.json").read_text())
    for name, results in t4_results.items():
        for voi, per in results.items():
            for label in METHODS:
                a, b = per[label], stored[name][voi][label]
                assert ("failed" in a) == ("failed" in b), (name, voi, label)
                if "failed" not in a:
                    assert a["tia_mbq_h"] == pytest.approx(b["tia_mbq_h"], rel=1e-6), (name, voi, label)
                    assert (a["ldm_fit_after"] is None) == (b["ldm_fit_after"] is None), (name, voi, label)
                    if a["ldm_fit_after"] is not None:
                        assert a["ldm_fit_after"] == pytest.approx(b["ldm_fit_after"], rel=1e-5), (name, voi, label)


def test_t4_f4a_with_the_participants(t4_results):
    """f4a (uptake then washout) is inside the participants' interquartile range for most VOIs (the
    liver of A, whose definition differs between participants, excluded)."""
    inside, total = 0, 0
    for name, results in t4_results.items():
        for voi, per in results.items():
            if voi == "Liver (no lesions)" or (voi == "Liver" and PATIENTS[name].challenge_id == "A"):
                continue
            r = per["f4a"]
            k = "ldm_scaled" if voi.startswith("Lesion") else "svalue"
            inside += r[f"{k}_in_iqr"]
            total += 1
            assert abs(r[f"{k}_diff_pct"]) < T4_TOLERANCE_PCT, (name, voi, r[f"{k}_diff_pct"])
            # fit after (LDM dose rate per time point with its CT): the same tolerance
            assert abs(r["ldm_fit_after_diff_pct"]) < T4_TOLERANCE_PCT, (name, voi, r["ldm_fit_after_diff_pct"])
    assert inside >= 10 and total == 14


def test_validation_page_is_up_to_date():
    # docs/8_validation_snmmi.md is generated from results/t5_radot.json: python validation/snmmi/report.py
    report = Path(__file__).resolve().parent / "report.py"
    result = subprocess.run([sys.executable, str(report), "--check"], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
