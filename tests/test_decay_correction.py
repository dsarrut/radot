"""Decay correction of the activity images (sidecar key decay_correction: start, none, admin).

An image decay-corrected to the injection time (admin) holds A(t) exp(lambda t): every computation that
uses the activity at its time point brings it back to the acquisition time. The key has no default:
missing is an error wherever the activity at a time point is used, not for TIA images."""
import json
import math
from datetime import datetime

import numpy as np
import pytest
import SimpleITK as sitk
from click.testing import CliRunner

from radot.check import Level, run_checks
from radot.cli import main
from radot.convert import convert_image_units
from radot.db import get_physical_half_life_h
from radot.dicom import convert_dicom
from radot.dose import compute_roi_dose, compute_roi_dose_mtp
from radot.meta import activity_at_acquisition_factor, get_decay_correction, write_metadata

from test_check import GOOD_META, activity_values, by_code, spect
from test_dicom import _write_dcm
from test_dose_input_errors import DATES, _labels_arr, _write
from test_dose_tia import make_image, make_labels

T_H = 96.0
LAMBDA_H = math.log(2) / get_physical_half_life_h("Lu177")
ROI = "Left kidney"
COMMON = dict(labels_map={1: ROI}, isotope="Lu177")


def _image(value, decay_correction, unit="bq"):
    image = make_image(value, unit, time_from_injection=T_H)
    if decay_correction is None:
        image.EraseMetaData("decay_correction")
    else:
        image.SetMetaData("decay_correction", decay_correction)
    return image


def _dose(image, method):
    params = {ROI: {"teff": 60.0}} if method == "madsen2018" else None
    res = compute_roi_dose(image=image, labels=make_labels(), method=method, method_params=params, **COMMON)
    if method == "ldm":
        res = res[0]
    return res[ROI]["mean_dose_gy"]


# ----------------------------------------------------------------------------------------------
# meta: write, read, validate
# ----------------------------------------------------------------------------------------------


def test_meta_cli_writes_decay_correction(tmp_path):
    path = tmp_path / "spect.nii.gz"
    sitk.WriteImage(make_image(1.0, "bq"), str(path))
    result = CliRunner().invoke(main, ["meta", str(path), "--decay-correction", "ADMIN"])
    assert result.exit_code == 0, result.output
    assert json.loads((tmp_path / "spect.json").read_text())["decay_correction"] == "admin"
    assert get_decay_correction(None, path) == "admin"
    result = CliRunner().invoke(main, ["meta", str(path), "--decay-correction", "later"])
    assert result.exit_code != 0


def test_write_metadata_invalid_decay_correction_errors(tmp_path):
    path = tmp_path / "spect.nii.gz"
    sitk.WriteImage(make_image(1.0, "bq"), str(path))
    with pytest.raises(ValueError, match="Invalid decay_correction 'later'"):
        write_metadata(path, decay_correction="later")


def test_get_decay_correction_missing_or_unknown_errors(tmp_path):
    path = tmp_path / "spect.nii.gz"
    sitk.WriteImage(make_image(1.0, "bq"), str(path))
    (tmp_path / "spect.json").write_text(json.dumps({"unit": "bq"}))
    with pytest.raises(ValueError, match="radot meta IMAGE --decay-correction start\\|none\\|admin"):
        get_decay_correction(None, path)
    (tmp_path / "spect.json").write_text(json.dumps({"unit": "bq", "decay_correction": "later"}))
    with pytest.raises(ValueError, match="Invalid decay_correction 'later'"):
        get_decay_correction(None, path)


def test_activity_at_acquisition_factor():
    assert activity_at_acquisition_factor("start", "Lu177", T_H) == 1.0
    assert activity_at_acquisition_factor("NONE", "Lu177", T_H) == 1.0
    assert activity_at_acquisition_factor("admin", "Lu177", T_H) == pytest.approx(math.exp(-LAMBDA_H * T_H))


# ----------------------------------------------------------------------------------------------
# radot dose, one time point
# ----------------------------------------------------------------------------------------------


@pytest.mark.parametrize("method", ["hanscheid2017", "hanscheid2018", "madsen2018", "ldm"])
def test_dose_admin_image_brought_back_to_acquisition(method):
    a_acq = 1e9  # Bq at the acquisition time
    d_start = _dose(_image(a_acq, "start"), method)
    d_none = _dose(_image(a_acq, "none"), method)
    d_admin = _dose(_image(a_acq * math.exp(LAMBDA_H * T_H), "admin"), method)
    assert d_none == pytest.approx(d_start, rel=1e-9)
    assert d_admin == pytest.approx(d_start, rel=1e-5)
    # the same values read as admin are a lower activity at the acquisition (≈ x0.66 at 96 h for Lu177)
    assert _dose(_image(a_acq, "admin"), method) == pytest.approx(d_start * math.exp(-LAMBDA_H * T_H), rel=1e-5)


def test_dose_missing_decay_correction_errors():
    with pytest.raises(ValueError, match="Decay correction of the image unknown"):
        _dose(_image(1e9, None), "hanscheid2018")


def test_dose_tia_image_needs_no_decay_correction():
    image = make_image(1e12, "bqs")
    assert not image.HasMetaDataKey("decay_correction")
    res = compute_roi_dose(image=image, labels=make_labels(), method="svalue", **COMMON)
    assert res[ROI]["mean_dose_gy"] > 0


def _cli_case(tmp_path, decay_correction):
    labels = _write(tmp_path / "labels.nii.gz", _labels_arr(), {"1": "Liver", "2": "Spleen"})
    sc = {"modality": "spect", "unit": "bqml", "time_from_injection": 24.0, "radionuclide": "Lu177", "sex": "male"}
    if decay_correction:
        sc["decay_correction"] = decay_correction
    arr = np.where(_labels_arr() > 0, 1e4, 0).astype(np.float32)
    spect_path = _write(tmp_path / "spect.nii.gz", arr, sc)
    return ["dose", str(spect_path), "-l", str(labels), "-m", "hanscheid2018"]


def test_dose_cli_decay_correction_from_sidecar_only(tmp_path):
    result = CliRunner().invoke(main, _cli_case(tmp_path, "start"))
    assert result.exit_code == 0, result.output


def test_dose_cli_decay_correction_missing_errors(tmp_path):
    result = CliRunner().invoke(main, _cli_case(tmp_path, None))
    assert result.exit_code != 0
    assert "--decay-correction" in result.output


# ----------------------------------------------------------------------------------------------
# radot dose, multiple time points
# ----------------------------------------------------------------------------------------------


def _mtp(tmp_path, decay_correction):
    folder = tmp_path / str(decay_correction)
    folder.mkdir()
    labels = _write(folder / "labels.nii.gz", _labels_arr(), {"1": "Liver", "2": "Spleen"})
    injection = "2024-01-01 10:00:00"
    spects = []
    for i, date in enumerate(DATES):
        t = (datetime.fromisoformat(date) - datetime.fromisoformat(injection)).total_seconds() / 3600
        a = 1e4 * math.exp(-(LAMBDA_H + 0.01) * t)  # activity at the acquisition time
        if decay_correction == "admin":
            a *= math.exp(LAMBDA_H * t)
        sc = {"modality": "spect", "unit": "bqml", "injection_date": injection, "acquisition_date": date,
              "radionuclide": "Lu177", "sex": "male"}
        if decay_correction:
            sc["decay_correction"] = decay_correction
        arr = np.where(_labels_arr() > 0, a, 0).astype(np.float32)
        spects.append(_write(folder / f"spect_tp{i+1}.nii.gz", arr, sc))
    return spects, labels


def test_mtp_admin_images_same_dose_as_start(tmp_path):
    def run(dc):
        spects, labels = _mtp(tmp_path, dc)
        return compute_roi_dose_mtp(spects, labels, method="jackson2013", labels_path=labels)["Liver"]

    start, admin = run("start"), run("admin")
    assert admin["activities_mbq"] == pytest.approx(start["activities_mbq"], rel=1e-5)
    assert admin["mean_dose_gy"] == pytest.approx(start["mean_dose_gy"], rel=1e-5)


def test_mtp_missing_decay_correction_errors(tmp_path):
    spects, labels = _mtp(tmp_path, None)
    with pytest.raises(ValueError, match="Time point 1: Decay correction of spect_tp1"):
        compute_roi_dose_mtp(spects, labels, method="jackson2013", labels_path=labels)


# ----------------------------------------------------------------------------------------------
# convert-units
# ----------------------------------------------------------------------------------------------


def _pet(tmp_path, **meta):
    path = tmp_path / "pet.nii.gz"
    sitk.WriteImage(sitk.GetImageFromArray(np.full((5, 5, 5), 1000.0, dtype=np.float32)), str(path))
    write_metadata(path, modality="pet", unit="bqml", weight=75.0, activity=250.0, radionuclide="F18",
                   injection_date="2021-04-28 10:00:00", acquisition_date="2021-04-28 11:49:46", **meta)
    return path


def test_convert_suv_missing_decay_correction_errors(tmp_path):
    path = _pet(tmp_path)
    with pytest.raises(ValueError, match="SUV/SUL conversion: Decay correction of pet.nii.gz unknown"):
        convert_image_units(sitk.ReadImage(str(path)), to_unit="suv", image_path=path)
    result = CliRunner().invoke(main, ["convert-units", str(path), "-t", "suv", "-o", str(tmp_path / "suv.nii.gz")])
    assert result.exit_code != 0 and "--decay-correction" in result.output


def test_convert_keeps_decay_correction(tmp_path):
    path = _pet(tmp_path, decay_correction="admin")
    out = tmp_path / "suv.nii.gz"
    result = CliRunner().invoke(main, ["convert-units", str(path), "-t", "suv", "-o", str(out)])
    assert result.exit_code == 0, result.output
    assert json.loads((tmp_path / "suv.json").read_text())["decay_correction"] == "admin"
    # in memory too
    converted = convert_image_units(sitk.ReadImage(str(path)), to_unit="bq", image_path=path)
    assert converted.GetMetaData("decay_correction") == "admin"


# ----------------------------------------------------------------------------------------------
# radot check
# ----------------------------------------------------------------------------------------------


def test_check_decay_correction(tmp_path):
    ok = spect(tmp_path, "ok")
    assert by_code(run_checks(ok), "meta.decay_correction", Level.OK)
    missing = spect(tmp_path, "missing", decay_correction=None)
    assert by_code(run_checks(missing), "meta.decay_correction", Level.FAIL)
    unknown = spect(tmp_path, "unknown", decay_correction="later")
    assert by_code(run_checks(unknown), "meta.decay_correction", Level.FAIL)
    tia = spect(tmp_path, "tia", unit="bqsml", decay_correction=None, values=1e9)
    assert not by_code(run_checks(tia), "meta.decay_correction")


def test_check_mixed_decay_corrections_warn(tmp_path):
    p1 = spect(tmp_path, "d1")
    p2 = spect(tmp_path, "d4", acquisition_date="2026-01-14 10:00:00", decay_correction="admin")
    assert by_code(run_checks([p1, p2]), "meta.decay_correction", Level.WARN)


def test_check_mtp_trend_on_admin_images(tmp_path):
    def admin(name, acquisition_date, t_h):
        meta = {**GOOD_META, "acquisition_date": acquisition_date, "decay_correction": "admin"}
        return spect(tmp_path, name, values=activity_values(0.5, meta) * math.exp(LAMBDA_H * t_h),
                     acquisition_date=acquisition_date, decay_correction="admin")

    findings = run_checks([admin("d1", "2026-01-11 10:00:00", 24.0), admin("d4", "2026-01-14 10:00:00", 96.0)])
    trend = by_code(findings, "activity.mtp_trend")
    assert trend and all(f.level == Level.OK for f in trend), trend


# ----------------------------------------------------------------------------------------------
# radot dicom
# ----------------------------------------------------------------------------------------------


def test_dicom_decay_correction_tag(tmp_path, caplog):
    tags = {"0008|0060": "NM", "0008|002a": "20260622100000", "0054|1001": "BQML"}
    _write_dcm(tmp_path / "a.dcm", tags)
    convert_dicom(tmp_path / "a.dcm", tmp_path / "a.nii.gz")
    assert "decay_correction" not in json.loads((tmp_path / "a.json").read_text())
    assert "no DICOM Decay Correction tag (0054,1102)" in caplog.text
    _write_dcm(tmp_path / "b.dcm", {**tags, "0054|1102": "START"})
    convert_dicom(tmp_path / "b.dcm", tmp_path / "b.nii.gz")
    assert json.loads((tmp_path / "b.json").read_text())["decay_correction"] == "start"
