import json

import numpy as np
import pytest
import SimpleITK as sitk
from click.testing import CliRunner

from radot.check import CHECKS, Level, load_thresholds, run_checks, worst_level
from radot.cli import main

GOOD_META = {
    "modality": "spect",
    "unit": "bqml",
    "radionuclide": "Lu177",
    "injection_date": "2026-01-10 10:00:00",
    "acquisition_date": "2026-01-11 10:00:00",
    "activity": 7400.0,
    "weight": 72.0,
    "height": 175.0,
    "sex": "male",
}


def write_image(folder, name, meta, value=1.0):
    path = folder / f"{name}.nii.gz"
    img = sitk.GetImageFromArray(np.full((4, 4, 4), value, dtype=np.float32))
    img.SetSpacing((4.0, 4.0, 4.0))
    sitk.WriteImage(img, str(path))
    if meta is not None:
        (folder / f"{name}.json").write_text(json.dumps(meta))
    return path


def spect(tmp_path, name="spect", **changes):
    meta = {**GOOD_META, **changes}
    meta = {k: v for k, v in meta.items() if v is not None}
    return write_image(tmp_path, name, meta)


def by_code(findings, code, level=None):
    return [f for f in findings if f.code == code and (level is None or f.level == level)]


def test_registry_and_thresholds():
    assert "meta.required" in CHECKS and "meta.dates_order" in CHECKS
    th = load_thresholds()
    assert th["injected_activity_mbq"]["lu177"][0] < 7400 < th["injected_activity_mbq"]["lu177"][1]


def test_good_case_clean(tmp_path):
    findings = run_checks(spect(tmp_path), method="madsen2018")
    assert findings
    assert worst_level(findings) == Level.OK


def test_no_sidecar(tmp_path):
    p = write_image(tmp_path, "nosidecar", None)
    assert by_code(run_checks(p), "meta.required", Level.FAIL)


def test_required_from_sidecar_only(tmp_path):
    # the dates are in the sidecar: OK; missing: FAIL with a method, WARN without
    assert not by_code(run_checks(spect(tmp_path), method="hanscheid2018"), "meta.required", Level.FAIL)
    p = spect(tmp_path, injection_date=None)
    fail = by_code(run_checks(p, method="hanscheid2018"), "meta.required", Level.FAIL)
    assert any("injection_date" in f.message for f in fail)
    warn = by_code(run_checks(p), "meta.required", Level.WARN)
    assert any("injection_date" in f.message for f in warn)


def test_required_suv_and_tia(tmp_path):
    p = spect(tmp_path, unit="suv", weight=None)
    assert any("weight" in f.message for f in by_code(run_checks(p), "meta.required", Level.FAIL))
    # TIA image: no dates needed, but only the TIA methods accept it
    p = spect(tmp_path, name="tia", unit="bqs", injection_date=None, acquisition_date=None)
    assert not by_code(run_checks(p, method="svalue"), "meta.required", Level.FAIL)
    assert by_code(run_checks(p, method="madsen2018"), "meta.required", Level.FAIL)


def test_dates(tmp_path):
    p = spect(tmp_path, acquisition_date="2026-01-09 10:00:00")
    assert by_code(run_checks(p), "meta.dates_order", Level.FAIL)
    p = spect(tmp_path, injection_date="10/01/2026")
    assert by_code(run_checks(p), "meta.dates_order", Level.FAIL)
    # 60 days for Lu177 (T½ 6.6 d) = 9 half-lives: ok; 80 days: warn
    p = spect(tmp_path, acquisition_date="2026-03-11 10:00:00")
    assert not by_code(run_checks(p), "meta.time_from_injection", Level.WARN)
    p = spect(tmp_path, acquisition_date="2026-03-31 10:00:00")
    assert by_code(run_checks(p), "meta.time_from_injection", Level.WARN)


def test_patient_and_activity(tmp_path):
    findings = run_checks(spect(tmp_path, weight=720.0))
    assert by_code(findings, "meta.patient", Level.WARN)
    # GBq typed instead of MBq
    findings = run_checks(spect(tmp_path, activity=7.4))
    assert by_code(findings, "meta.activity", Level.WARN)
    # user thresholds (merged with the defaults)
    findings = run_checks(spect(tmp_path, activity=7.4), thresholds={"injected_activity_mbq": {"lu177": [1, 10]}})
    assert not by_code(findings, "meta.activity", Level.WARN)


def test_radionuclide_and_modality(tmp_path):
    assert by_code(run_checks(spect(tmp_path, radionuclide="Xx999")), "meta.radionuclide", Level.FAIL)
    ct = write_image(tmp_path, "ct", {"modality": "ct", "unit": "bqml"})
    assert by_code(run_checks(ct=ct), "meta.modality", Level.FAIL)
    ct = write_image(tmp_path, "ct", {"modality": "ct", "unit": "hu"})
    assert worst_level(run_checks(ct=ct)) == Level.OK


def test_multi_time_points(tmp_path):
    p1 = spect(tmp_path, "d1")
    p2 = spect(tmp_path, "d4", acquisition_date="2026-01-14 10:00:00")
    assert worst_level(run_checks([p1, p2], method="jackson2013")) == Level.OK
    p3 = spect(tmp_path, "d7", injection_date="2026-01-11 10:00:00", acquisition_date="2026-01-17 10:00:00")
    assert by_code(run_checks([p1, p3]), "meta.consistency", Level.FAIL)
    p4 = spect(tmp_path, "d1bis")
    assert by_code(run_checks([p1, p4]), "meta.time_from_injection", Level.WARN)


def test_only_skip_and_errors(tmp_path):
    p = spect(tmp_path)
    assert {f.code for f in run_checks(p, only=["meta.patient"])} == {"meta.patient"}
    assert not [f for f in run_checks(p, skip=["meta"]) if f.category == "meta"]
    with pytest.raises(ValueError):
        run_checks(p, method="nope")


def test_cli(tmp_path):
    runner = CliRunner()
    good = spect(tmp_path)
    r = runner.invoke(main, ["check", str(good), "-m", "madsen2018"])
    assert r.exit_code == 0, r.output
    assert "0 failures" in r.output

    bad = spect(tmp_path, "bad", acquisition_date="2026-01-09 10:00:00")
    report = tmp_path / "qc.json"
    r = runner.invoke(main, ["check", str(bad), "-o", str(report)])
    assert r.exit_code == 1
    assert "meta.dates_order" in r.output
    assert any(f["code"] == "meta.dates_order" and f["level"] == "FAIL" for f in json.loads(report.read_text()))

    warn = spect(tmp_path, "warn", weight=720.0)
    assert runner.invoke(main, ["check", str(warn)]).exit_code == 0
    assert runner.invoke(main, ["check", str(warn), "--strict"]).exit_code == 1

    r = runner.invoke(main, ["check", str(good), "-a", "--only", "meta.patient"])
    assert "OK" in r.output and "72 kg" in r.output

    r = runner.invoke(main, ["check", "--list"])
    assert r.exit_code == 0 and "meta.required" in r.output
