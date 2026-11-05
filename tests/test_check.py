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


LU177_HL_H = 6.647 * 24
VOXEL_ML = 0.064  # 4 mm voxels
N = 4


def write_image(folder, name, meta, values=1.0):
    path = folder / f"{name}.nii.gz"
    arr = np.broadcast_to(np.asarray(values, dtype=np.float32), (N, N, N))
    img = sitk.GetImageFromArray(np.ascontiguousarray(arr))
    img.SetSpacing((4.0, 4.0, 4.0))
    sitk.WriteImage(img, str(path))
    if meta is not None:
        (folder / f"{name}.json").write_text(json.dumps(meta))
    return path


def activity_values(fraction, meta):
    """Bq/mL values (not integers) holding a fraction of the injected activity, decayed to the acquisition."""
    t = 24.0
    if meta.get("injection_date") and meta.get("acquisition_date"):
        from datetime import datetime

        f = "%Y-%m-%d %H:%M:%S"
        try:
            t = (datetime.strptime(meta["acquisition_date"], f) - datetime.strptime(meta["injection_date"], f)).total_seconds() / 3600
        except ValueError:
            pass
    total_bq = fraction * meta.get("activity", 7400.0) * 1e6 * 0.5 ** (t / LU177_HL_H)
    noise = 1 + 0.01 * np.random.default_rng(0).standard_normal((N, N, N))
    noise /= noise.mean()
    return total_bq / (N**3 * VOXEL_ML) * noise


def spect(tmp_path, name="spect", fraction=0.5, values=None, **changes):
    meta = {**GOOD_META, **changes}
    meta = {k: v for k, v in meta.items() if v is not None}
    if values is None:
        values = activity_values(fraction, meta)
    return write_image(tmp_path, name, meta, values)


def ct_values(offset=0.0, scale=1.0):
    """Half air (-1000), half soft tissue (40 HU), one bone voxel."""
    arr = np.full((N, N, N), 40.0)
    arr[: N // 2] = -1000.0
    arr[-1, -1, -1] = 1200.0
    return arr * scale + offset


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
    ct = write_image(tmp_path, "ct", {"modality": "ct", "unit": "bqml"}, ct_values())
    assert by_code(run_checks(ct=ct), "meta.modality", Level.FAIL)
    ct = write_image(tmp_path, "ct", {"modality": "ct", "unit": "hu"}, ct_values())
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


def test_ct_hu(tmp_path):
    meta = {"modality": "ct", "unit": "hu"}
    assert worst_level(run_checks(ct=write_image(tmp_path, "ct", meta, ct_values()))) == Level.OK
    # padding outside the field of view is ignored
    v = ct_values()
    v[0, 0, 0] = -3024
    assert worst_level(run_checks(ct=write_image(tmp_path, "ct", meta, v))) == Level.OK
    assert by_code(run_checks(ct=write_image(tmp_path, "ct", meta, ct_values(offset=1024))), "units.ct_hu", Level.FAIL)
    mu = (ct_values() + 1000) / 1000 * 0.0096  # attenuation map, 1/mm
    assert by_code(run_checks(ct=write_image(tmp_path, "ct", meta, mu)), "units.ct_hu", Level.FAIL)
    no_air = ct_values()
    no_air[no_air < 0] = 40
    assert by_code(run_checks(ct=write_image(tmp_path, "ct", meta, no_air)), "units.ct_hu", Level.WARN)


def test_activity_values(tmp_path):
    counts = spect(tmp_path, "counts", values=np.arange(N**3).reshape(N, N, N) + 1.0)
    assert by_code(run_checks(counts), "units.activity_values", Level.WARN)
    assert by_code(run_checks(spect(tmp_path, "zero", values=0.0)), "units.activity_values", Level.FAIL)
    neg = activity_values(0.5, GOOD_META)
    neg[0] = -neg[0]
    assert by_code(run_checks(spect(tmp_path, "neg", values=neg)), "units.activity_values", Level.WARN)
    assert by_code(run_checks(spect(tmp_path, "suv", unit="suv", values=500.5)), "units.suv_range", Level.WARN)


def test_dose_range(tmp_path):
    ok = write_image(tmp_path, "dose", {"modality": "dose", "unit": "gy"}, 12.5)
    assert worst_level(run_checks(dose=ok)) == Level.OK
    high = write_image(tmp_path, "dose", {"modality": "dose", "unit": "gy"}, 5000.5)
    assert by_code(run_checks(dose=high), "units.dose_range", Level.WARN)


def test_fov_fraction(tmp_path):
    f = by_code(run_checks(spect(tmp_path, fraction=0.5)), "activity.fov_fraction")
    assert f[0].level == Level.OK and f[0].value == pytest.approx(0.5, rel=1e-3)
    assert by_code(run_checks(spect(tmp_path, fraction=1.5)), "activity.fov_fraction", Level.FAIL)
    assert by_code(run_checks(spect(tmp_path, fraction=0.001)), "activity.fov_fraction", Level.FAIL)
    assert by_code(run_checks(spect(tmp_path, fraction=0.05)), "activity.fov_fraction", Level.WARN)
    # image decay-corrected to the injection time: compared with the injected activity, no decay
    meta = {**GOOD_META, "decay_correction": "admin"}
    values = activity_values(0.5, meta) / 0.5 ** (24 / LU177_HL_H)
    f = by_code(run_checks(write_image(tmp_path, "admin", meta, values)), "activity.fov_fraction")
    assert f[0].value == pytest.approx(0.5, rel=1e-3)
    # no injected activity: skipped (meta.required reports it)
    f = by_code(run_checks(spect(tmp_path, activity=None, values=1.5)), "activity.fov_fraction")
    assert f[0].level == Level.INFO
    # SUV: sum(SUV x mL) / weight (g)
    suv = 0.5 * 72000 / (N**3 * VOXEL_ML)
    f = by_code(run_checks(spect(tmp_path, "suv", unit="suv", values=suv)), "activity.fov_fraction")
    assert f[0].value == pytest.approx(0.5, rel=1e-3)
    # TIA (Bq.s/mL) below the bound (injected activity x physical mean life)
    tia = 0.3 * 7400e6 * LU177_HL_H * 3600 / np.log(2) / (N**3 * VOXEL_ML)
    p = spect(tmp_path, "tia", unit="bqsml", injection_date=None, acquisition_date=None, values=tia)
    f = by_code(run_checks(p), "activity.fov_fraction")
    assert f[0].level == Level.OK and f[0].value == pytest.approx(0.3, rel=1e-3)


def test_mtp_trend(tmp_path):
    # biological clearance: faster than the physical decay -> ok
    d1 = spect(tmp_path, "d1", fraction=0.5)
    d4 = spect(tmp_path, "d4", fraction=0.2, acquisition_date="2026-01-14 10:00:00")
    assert by_code(run_checks([d1, d4]), "activity.mtp_trend", Level.OK)
    # same activity at day 4 as at day 1 = images decay-corrected (or recalibrated)
    v = activity_values(0.5, GOOD_META)
    flat = spect(tmp_path, "d4flat", values=v, acquisition_date="2026-01-14 10:00:00")
    assert by_code(run_checks([d1, flat]), "activity.mtp_trend", Level.WARN)
