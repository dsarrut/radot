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


ROI_NAMES = {"1": "Liver", "2": "Spleen", "3": "Kidney_R", "4": "Kidney_L", "5": "Lung_L", "6": "Thyroid"}


def body(tmp_path, name="labels", mirror=False, names=ROI_NAMES, liver_slices=5):
    """Labels and CT, 20 mm voxels (8 mL), shape z, y, x; +x = patient left (LPS)."""
    lab = np.zeros((12, 12, 22), dtype=np.uint8)  # one empty voxel at the border (not truncated)
    inner = lab[1:-1, 1:-1, 1:-1]
    inner[:liver_slices, :, 0:5] = 1  # liver, patient right, 250 voxels = 2000 mL
    inner[0, 0:4, 15:20] = 2  # spleen, 160 mL
    inner[5:7, 0:5, 5:7] = 3  # right kidney, 160 mL
    inner[5:7, 0:5, 13:15] = 4  # left kidney
    inner[5:10, 5:10, 10:20] = 5  # left lung
    ct = np.full(lab.shape, 40.0)
    ct[lab == 5] = -800.0
    if mirror:
        lab = lab[:, :, ::-1]
    paths = []
    for suffix, arr, meta in (("", lab, names), ("_ct", ct, {"modality": "ct", "unit": "hu"})):
        img = sitk.GetImageFromArray(np.ascontiguousarray(arr))
        img.SetSpacing((20.0, 20.0, 20.0))
        path = tmp_path / f"{name}{suffix}.nii.gz"
        sitk.WriteImage(img, str(path))
        if meta is not None:
            (tmp_path / f"{name}{suffix}.json").write_text(json.dumps(meta))
        paths.append(path)
    return paths


def test_roi_clean(tmp_path):
    labels, ct = body(tmp_path)
    findings = run_checks(ct=ct, labels=labels, only=["roi"])
    assert worst_level(findings) == Level.INFO  # Thyroid empty: reported, not a problem
    assert any("Thyroid" in f.message for f in by_code(findings, "roi.empty", Level.INFO))
    mass = by_code(findings, "roi.mass", Level.OK)
    assert mass and "4 organs" in mass[0].message and "sex unknown" in mass[0].message


def test_roi_laterality(tmp_path):
    labels, ct = body(tmp_path, mirror=True)
    fail = by_code(run_checks(labels=labels, only=["roi.laterality"]), "roi.laterality", Level.FAIL)
    assert len(fail) == 2  # kidney pair swapped, liver and spleen on the wrong side
    names = {**ROI_NAMES, "3": "Kidney_L", "4": "Kidney_R"}
    labels, ct = body(tmp_path, "swapped", names=names)
    fail = by_code(run_checks(labels=labels, only=["roi.laterality"]), "roi.laterality", Level.FAIL)
    assert len(fail) == 1 and "Kidney_L" in fail[0].message


def test_roi_mass_density_names(tmp_path):
    labels, ct = body(tmp_path, liver_slices=1)  # 400 mL liver
    warn = by_code(run_checks(ct=ct, labels=labels), "roi.mass", Level.WARN)
    assert [f.subject for f in warn] == ["Liver"]
    # labels misaligned with the CT: the lung label on soft tissue
    labels, ct = body(tmp_path, "shifted")
    img = sitk.ReadImage(str(ct))
    img = sitk.GetImageFromArray(np.full(sitk.GetArrayFromImage(img).shape, 40.0))
    img.SetSpacing((20.0, 20.0, 20.0))
    sitk.WriteImage(img, str(ct))
    warn = by_code(run_checks(ct=ct, labels=labels), "roi.density", Level.WARN)
    assert warn and "Lung_L" in warn[0].value
    # no sidecar
    (tmp_path / "labels.json").unlink()
    # no labels sidecar: radot dose stops on it, so a FAIL with the reason
    fail = by_code(run_checks(labels=tmp_path / "labels.nii.gz"), "roi.names", Level.FAIL)
    assert fail and "has no sidecar" in fail[0].message


def test_roi_mass_female(tmp_path):
    labels, ct = body(tmp_path)
    act = spect(tmp_path, sex="female")
    ok = by_code(run_checks(act, ct=ct, labels=labels, only=["roi.mass"]), "roi.mass", Level.OK)
    assert "(female)" in ok[0].message


def geom_case(tmp_path, name="g", label_shift=(0, 0, 0), act_shift=(0, 0, 0), left_uptake=1.0):
    """Body (CT) with two hot kidneys, blurred activity, labels; 4 mm voxels, arrays z, y, x."""
    shape = (40, 40, 40)
    ct = np.full(shape, -1000.0)
    ct[5:35, 8:32, 6:34] = 40.0
    lab = np.zeros(shape, dtype=np.uint8)
    lab[15:25, 18:26, 10:15] = 1  # Kidney_R, patient right (low x)
    lab[15:25, 18:26, 25:30] = 2  # Kidney_L
    act = np.where(ct > 0, 1.0, 0.0) + 20.0 * (lab == 1) + 20.0 * left_uptake * (lab == 2)
    act = np.roll(act, act_shift, axis=(0, 1, 2))
    lab = np.roll(lab, label_shift, axis=(0, 1, 2))
    paths = {}
    for kind, arr, meta in (
        ("ct", ct, {"modality": "ct", "unit": "hu"}),
        ("labels", lab, {"1": "Kidney_R", "2": "Kidney_L"}),
        ("spect", act * 1e4, {"modality": "spect", "unit": "bqml"}),
    ):
        img = sitk.GetImageFromArray(np.ascontiguousarray(arr))
        img.SetSpacing((4.0, 4.0, 4.0))
        if kind == "spect":
            img = sitk.SmoothingRecursiveGaussian(img, 5.0)
        path = tmp_path / f"{name}_{kind}.nii.gz"
        sitk.WriteImage(img, str(path))
        (tmp_path / f"{name}_{kind}.json").write_text(json.dumps(meta))
        paths[kind] = path
    return paths


def run_geom(p):
    return run_checks(p["spect"], ct=p["ct"], labels=p["labels"], only=["geom", "activity.lr_symmetry"])


def test_geom_aligned(tmp_path):
    findings = run_geom(geom_case(tmp_path))
    assert worst_level(findings) == Level.OK
    ok = by_code(findings, "geom.organ_shift", Level.OK)
    assert ok and "Kidney_L 0 mm" in ok[0].message


def test_geom_organ_shift(tmp_path):
    warn = by_code(run_geom(geom_case(tmp_path, label_shift=(0, 6, 0))), "geom.organ_shift", Level.WARN)
    # labels moved 24 mm posterior (+y): the activity is 24 mm anterior of them
    assert warn and "24 mm anterior" in warn[0].message


def test_geom_activity_in_body_and_overlap(tmp_path):
    # activity moved 15 voxels (60 mm) to the left: partly outside the body
    findings = run_geom(geom_case(tmp_path, act_shift=(0, 0, 15)))
    assert by_code(findings, "geom.activity_in_body", Level.WARN) or by_code(findings, "geom.activity_in_body", Level.FAIL)
    # kidneys partly outside the activity image
    p = geom_case(tmp_path, "crop")
    img = sitk.ReadImage(str(p["spect"]))
    sitk.WriteImage(img[:, :, :20], str(p["spect"]))  # z 0-19: kidneys z 15-24 half outside
    assert by_code(run_geom(p), "geom.overlap", Level.WARN)


def test_lr_symmetry_and_direction(tmp_path):
    warn = by_code(run_geom(geom_case(tmp_path, left_uptake=0.1)), "activity.lr_symmetry", Level.WARN)
    assert warn and "Kidney_L" in warn[0].message
    p = geom_case(tmp_path, "obl")
    img = sitk.ReadImage(str(p["ct"]))
    c, s = np.cos(0.3), np.sin(0.3)
    img.SetDirection((c, -s, 0, s, c, 0, 0, 0, 1))
    sitk.WriteImage(img, str(p["ct"]))
    assert by_code(run_checks(ct=p["ct"], only=["geom.direction"]), "geom.direction", Level.WARN)


# The time from injection follows the dose rule (meta.get_time_from_injection_h)

def test_time_key_without_dates_accepted(tmp_path):
    """A time_from_injection key without dates: radot dose uses it, so check does too."""
    p = spect(tmp_path, injection_date=None, acquisition_date=None, time_from_injection=24.0)
    findings = run_checks(p, method="hanscheid2018")
    assert not by_code(findings, "meta.required", Level.FAIL)
    ok = by_code(findings, "meta.time_from_injection", Level.OK)
    assert ok and ok[0].value == pytest.approx(24.0)


def test_time_key_disagrees_with_dates_fails(tmp_path):
    p = spect(tmp_path, time_from_injection=48.0)  # the dates give 24 h
    fail = by_code(run_checks(p), "meta.time_from_injection", Level.FAIL)
    assert fail and "differs from the dates" in fail[0].message


def test_time_key_agrees_with_dates_ok(tmp_path):
    p = spect(tmp_path, time_from_injection=24.0)
    assert by_code(run_checks(p), "meta.time_from_injection", Level.OK)


def test_unreadable_sidecar_one_fail_other_images_checked(tmp_path):
    """Two time points, one sidecar broken: one FAIL for it, the other image is still checked."""
    good = spect(tmp_path, name="tp1")
    bad = spect(tmp_path, name="tp2")
    (tmp_path / "tp2.json").write_text("{ broken")
    findings = run_checks([good, bad])
    assert [f.subject for f in by_code(findings, "meta.sidecar", Level.FAIL)] == ["tp2.nii.gz"]
    assert by_code(findings, "meta.required", Level.OK)[0].subject == "tp1.nii.gz"
    assert not [f for f in findings if f.subject == "tp2.nii.gz" and f.category == "meta" and f.code != "meta.sidecar"]
