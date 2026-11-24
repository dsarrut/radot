"""radot integrate on dose-rate images (fit after, ROI level): the sidecar of `radot dose --dose-rate`,
the mass-weighted ROI mean dose rate, the equivalence with fit first under a fixed geometry, errors."""

import json
import logging
import math

import numpy as np
import pytest
import SimpleITK as sitk
from click.testing import CliRunner

from radot.cli import main
from radot.db import get_cascade_energy_constant_mev, get_physical_half_life_h
from radot.dose import compute_roi_dose_tia_table
from radot.integrate import integrate_rois, read_integrate

LP = math.log(2) / get_physical_half_life_h("Lu177")
TIMES = [4.0, 24.0, 96.0, 168.0]
LIVER = (slice(2, 5), slice(2, 5), slice(2, 5))
VOX_ML = 0.008  # 2 mm voxels
MEV_J_G = 1.60218e-10  # MeV -> J, x 1000 (g -> kg)


def _write(path, arr, sidecar):
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing((2.0, 2.0, 2.0))
    sitk.WriteImage(img, str(path))
    if sidecar is not None:
        path.with_name(path.name.split(".")[0] + ".json").write_text(json.dumps(sidecar))
    return path


def _liver_mbq(t):
    return 100.0 * math.exp(-(0.02 + LP) * t)


def _study(tmp_path, decay_correction="start", liver_hu=(50.0, 50.0, 50.0, 50.0)):
    """Activity images (bqml, liver f2: 100 MBq, 0.02/h), one label image, one CT per time point
    (liver HU per time point, 0 HU elsewhere)."""
    lab = np.zeros((8, 8, 8), np.uint8)
    lab[LIVER] = 1
    labels = _write(tmp_path / "labels.nii.gz", lab, {"1": "Liver"})
    images, cts = [], []
    for i, t in enumerate(TIMES):
        a = _liver_mbq(t) * (math.exp(LP * t) if decay_correction == "admin" else 1.0)
        arr = np.zeros((8, 8, 8), np.float32)
        arr[LIVER] = a * 1e6 / (27 * VOX_ML)
        images.append(_write(tmp_path / f"spect_tp{i + 1}.nii.gz", arr, {
            "modality": "spect", "unit": "bqml", "radionuclide": "Lu177", "time_from_injection": t,
            "decay_correction": decay_correction, "sex": "male", "weight": 70}))
        ct = np.zeros((8, 8, 8), np.float32)
        ct[LIVER] = liver_hu[i]
        cts.append(_write(tmp_path / f"ct_tp{i + 1}.nii.gz", ct, {"modality": "ct", "unit": "hu"}))
    return images, labels, cts


def _dose_rates(tmp_path, images, cts):
    """radot dose -m ldm --dose-rate at each time point, with its CT."""
    out = []
    for i, (img, ct) in enumerate(zip(images, cts)):
        dr = tmp_path / f"dr_tp{i + 1}.nii.gz"
        r = CliRunner().invoke(main, ["dose", str(img), "--ct", str(ct), "-m", "ldm", "--dose-rate",
                                      "-o", str(dr)])
        assert r.exit_code == 0, r.output
        out.append(dr)
    return out


def _gy_per_mbq_h(hu):
    """Mean dose (Gy) per MBq.h of the liver (27 voxels, uniform activity and density)."""
    delta = get_cascade_energy_constant_mev("Lu177", cascade=False)[0]
    mass_g = 27 * VOX_ML * (1.0 + hu / 1000.0)
    return 1e6 * 3600.0 * delta * MEV_J_G / mass_g


@pytest.mark.parametrize("dc", ["start", "admin"])
def test_dose_rate_sidecar(tmp_path, dc):
    images, labels, cts = _study(tmp_path, decay_correction=dc)
    dr = _dose_rates(tmp_path, images[:1], cts[:1])[0]
    sc = json.loads(dr.with_name("dr_tp1.json").read_text())
    assert sc["unit"] == "gy/h" and sc["modality"] == "dose" and sc["decay_correction"] == "none"
    assert sc["radionuclide"] == "Lu177" and sc["sex"] == "male" and sc["weight"] == 70
    assert sc["time_from_injection"] == pytest.approx(TIMES[0])
    # at the acquisition time, whatever the decay correction of the activity image
    arr = sitk.GetArrayFromImage(sitk.ReadImage(str(dr)))
    expected = _liver_mbq(TIMES[0]) * _gy_per_mbq_h(50.0)
    assert float(arr[LIVER].mean()) == pytest.approx(expected, rel=1e-5)


def test_dose_rate_unreadable_activity_sidecar_is_an_error(tmp_path):
    images, labels, cts = _study(tmp_path)
    images[0].with_name("spect_tp1.json").write_text("{not json")
    r = CliRunner().invoke(main, ["dose", str(images[0]), "--ct", str(cts[0]), "-m", "ldm", "--dose-rate",
                                  "--isotope", "Lu177", "-o", str(tmp_path / "dr.nii.gz")])
    assert r.exit_code != 0 and "cannot be read" in r.output


@pytest.mark.parametrize("opts", [
    {"model": "f2", "error_model": "proportional"},
    {"model": "f2", "error_model": "fsd:0.1"},
    {"model": "trapezoid", "head": "ramp", "tail": "phys"},
    {"model": "trapezoid", "head": "extrap", "tail": "fit", "tail_points": 2},
])
def test_fit_after_equals_fit_first_with_a_fixed_geometry(tmp_path, opts):
    # same CT at every time point: the ROI dose rate is proportional to the ROI activity, the fit is
    # scale invariant (relative or fractional error model, trapezoid): fit after = fit first x Gy/MBq.h
    images, labels, cts = _study(tmp_path)
    first = integrate_rois(images, labels, **opts)["rois"]["Liver"]
    after_doc = integrate_rois(_dose_rates(tmp_path, images, cts), labels, ct=cts, **opts)
    assert after_doc["quantity"] == "dose_rate"
    after = after_doc["rois"]["Liver"]
    k = _gy_per_mbq_h(50.0)
    assert after["mean_dose_gy"] == pytest.approx(first["tia_mbq_h"] * k, rel=1e-5)
    if opts.get("error_model", "").startswith("fsd"):  # (proportional: the SE of exact data is rounding)
        assert after["mean_dose_se_gy"] == pytest.approx(first["tia_se_mbq_h"] * k, rel=1e-4)
    assert after["masses_g"] == pytest.approx([27 * VOX_ML * 1.05] * 4)
    fit = after["fits"][after["model"]]
    assert not any(key.endswith("mbq_h") for key in fit)
    assert ("dose_gy" in fit) if opts["model"] == "f2" else ("parts_gy" in fit)
    # the ROI dose rate at the first time point: LDM of the ROI activity, the same mass
    assert after["dose_rates_gy_h"][0] == pytest.approx(_liver_mbq(TIMES[0]) * k, rel=1e-5)


def test_fit_after_uses_the_geometry_of_each_time_point(tmp_path):
    hu = (0.0, 100.0, 200.0, 300.0)
    images, labels, cts = _study(tmp_path, liver_hu=hu)
    doc = integrate_rois(_dose_rates(tmp_path, images, cts), labels, ct=cts, model="trapezoid",
                         head="ramp", tail="phys")
    liver = doc["rois"]["Liver"]
    assert liver["masses_g"] == pytest.approx([27 * VOX_ML * (1 + h / 1000) for h in hu])
    rates = [_liver_mbq(t) * _gy_per_mbq_h(h) for t, h in zip(TIMES, hu)]
    assert liver["dose_rates_gy_h"] == pytest.approx(rates, rel=1e-5)


def test_fit_after_without_ct_warns_water(tmp_path, caplog):
    images, labels, cts = _study(tmp_path, liver_hu=(0.0,) * 4)
    drs = _dose_rates(tmp_path, images, cts)
    with caplog.at_level(logging.WARNING, logger="radot"):
        doc = integrate_rois(drs, labels, model="f2", error_model="fsd:0.1")
    assert any("water density" in r.message for r in caplog.records)
    assert doc["rois"]["Liver"]["masses_g"] == pytest.approx([27 * VOX_ML] * 4)


def test_dose_rate_in_gy_per_second(tmp_path):
    images, labels, cts = _study(tmp_path)
    drs = _dose_rates(tmp_path, images, cts)
    ref = integrate_rois(drs, labels, ct=cts, model="trapezoid", head="ramp", tail="phys")
    for p in drs:  # the same images in Gy/s
        img = sitk.ReadImage(str(p))
        sitk.WriteImage(img / 3600.0, str(p))
        sc = p.with_name(p.name.split(".")[0] + ".json")
        sc.write_text(json.dumps({**json.loads(sc.read_text()), "unit": "gy/s"}))
    doc = integrate_rois(drs, labels, ct=cts, model="trapezoid", head="ramp", tail="phys")
    assert doc["rois"]["Liver"]["mean_dose_gy"] == pytest.approx(ref["rois"]["Liver"]["mean_dose_gy"], rel=1e-5)


def test_mixed_quantities_and_ct_with_activity_are_errors(tmp_path):
    images, labels, cts = _study(tmp_path)
    drs = _dose_rates(tmp_path, images[:2], cts[:2])
    with pytest.raises(ValueError, match="cannot be mixed"):
        integrate_rois(drs + images[2:], labels, model="trapezoid", head="ramp", tail="phys")
    with pytest.raises(ValueError, match="CT is used only with dose-rate images"):
        integrate_rois(images, labels, ct=cts, model="trapezoid", head="ramp", tail="phys")
    with pytest.raises(ValueError, match="2 CT images for 4 images"):
        integrate_rois(_dose_rates(tmp_path, images, cts), labels, ct=cts[:2], model="trapezoid",
                       head="ramp", tail="phys")


def test_dose_rate_without_decay_correction_is_an_error(tmp_path):
    images, labels, cts = _study(tmp_path)
    drs = _dose_rates(tmp_path, images, cts)
    sc = drs[0].with_name("dr_tp1.json")
    sc.write_text(json.dumps({k: v for k, v in json.loads(sc.read_text()).items() if k != "decay_correction"}))
    with pytest.raises(ValueError, match="(?i)decay"):
        integrate_rois(drs, labels, ct=cts, model="trapezoid", head="ramp", tail="phys")


def test_dose_on_a_dose_rate_table_is_an_error(tmp_path):
    images, labels, cts = _study(tmp_path)
    doc = integrate_rois(_dose_rates(tmp_path, images, cts), labels, ct=cts, model="trapezoid", head="ramp",
                         tail="phys")
    with pytest.raises(ValueError, match="already gives the mean absorbed dose"):
        compute_roi_dose_tia_table(doc, labels, method="svalue")


def test_cli_fit_after(tmp_path):
    images, labels, cts = _study(tmp_path)
    drs = _dose_rates(tmp_path, images, cts)
    out = tmp_path / "dose_rois.json"
    args = ["integrate", *map(str, drs), "-l", str(labels), "--model", "f2", "--error-model", "fsd:0.1",
            "-o", str(out)]
    for c in cts:
        args += ["--ct", str(c)]
    r = CliRunner().invoke(main, args)
    assert r.exit_code == 0, r.output
    assert "Dose (Gy)" in r.output
    doc = read_integrate(out)
    assert doc["quantity"] == "dose_rate"
    assert doc["rois"]["Liver"]["mean_dose_gy"] == pytest.approx(
        100.0 / (0.02 + LP) * _gy_per_mbq_h(50.0), rel=1e-5)
    # radot dose on it: an error, it is already a dose
    r = CliRunner().invoke(main, ["dose", str(out), "-m", "svalue", "-l", str(labels)])
    assert r.exit_code != 0 and "already gives the mean absorbed dose" in r.output
