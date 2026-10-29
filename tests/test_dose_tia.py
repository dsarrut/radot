import json

import numpy as np
import pytest
import SimpleITK as sitk
from click.testing import CliRunner

from radot.cli import main
from radot.db import get_physical_half_life_h
from radot.dose import compute_roi_dose

# Mocked OpenDose values (tests/conftest.py), Lu177 male phantom
S_KIDNEY_L = 1.5662e-4  # mGy/MBq.s
M_KIDNEY_L = 152.983  # g


def make_image(value, unit, time_from_injection=None):
    """5x5x5 image, 10 mm voxels (1 mL), one non-zero voxel at the centre."""
    image = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    image.SetSpacing((10.0, 10.0, 10.0))
    image.SetMetaData("unit", unit)
    if time_from_injection is not None:
        image.SetMetaData("time_from_injection", str(time_from_injection))
    image[2, 2, 2] = value
    return image


def make_labels():
    labels = sitk.Image([5, 5, 5], sitk.sitkUInt16)
    labels.SetSpacing((10.0, 10.0, 10.0))
    labels[2, 2, 2] = 1
    return labels


def test_svalue_tia_bqs():
    tia_bqs = 1e12
    res = compute_roi_dose(
        image=make_image(tia_bqs, "bqs"),
        labels=make_labels(),
        method="svalue",
        isotope="Lu177",
        labels_map={1: "Left kidney"},
    )
    r = res["Left kidney"]
    # ROI mass: 1 voxel of 1 mL, no CT (water)
    expected = (tia_bqs / 1e6) * S_KIDNEY_L / (1.0 / M_KIDNEY_L) / 1000.0
    assert r["mean_dose_gy"] == pytest.approx(expected, rel=1e-6)
    assert r["tia_mbq_s"] == pytest.approx(tia_bqs / 1e6)
    assert r["tia_mbq_h"] == pytest.approx(tia_bqs / 1e6 / 3600.0)
    assert "roi_teff_h" not in r


def test_svalue_tia_bqsml_same_as_bqs():
    # 1 mL voxel: Bq.s/mL and Bq.s per voxel are the same numbers
    common = dict(labels=make_labels(), method="svalue", isotope="Lu177", labels_map={1: "Left kidney"})
    d_bqs = compute_roi_dose(image=make_image(1e12, "bqs"), **common)
    d_bqsml = compute_roi_dose(image=make_image(1e12, "bqsml"), **common)
    assert d_bqsml["Left kidney"]["mean_dose_gy"] == pytest.approx(d_bqs["Left kidney"]["mean_dose_gy"])

    # 8 mm voxels (0.512 mL): bqsml is multiplied by the voxel volume
    img = make_image(1e12, "bqsml")
    img.SetSpacing((8.0, 8.0, 8.0))
    lab = make_labels()
    lab.SetSpacing((8.0, 8.0, 8.0))
    res = compute_roi_dose(image=img, labels=lab, method="svalue", isotope="Lu177", labels_map={1: "Left kidney"})
    assert res["Left kidney"]["tia_mbq_s"] == pytest.approx(1e12 * 0.512 / 1e6)


def test_svalue_tia_equals_madsen2018():
    # Madsen 2018 integrates A(t) mono-exponentially: TIA = A(t) exp(k t) / k
    a_bq, t_h, teff_h = 1e8, 24.0, 45.0
    k = np.log(2) / (teff_h * 3600.0)
    tia_bqs = a_bq * np.exp(k * t_h * 3600.0) / k

    madsen = compute_roi_dose(
        image=make_image(a_bq, "bq", time_from_injection=t_h),
        labels=make_labels(),
        method="madsen2018",
        isotope="Lu177",
        labels_map={1: "Left kidney"},
        method_params={"Left kidney": {"teff": teff_h}},
    )
    svalue = compute_roi_dose(
        image=make_image(tia_bqs, "bqs"),
        labels=make_labels(),
        method="svalue",
        isotope="Lu177",
        labels_map={1: "Left kidney"},
    )
    assert svalue["Left kidney"]["mean_dose_gy"] == pytest.approx(
        madsen["Left kidney"]["mean_dose_gy"], rel=1e-5
    )


def test_ldm_tia_equals_ldm_activity():
    # LDM on activity integrates with physical decay only: TIA = A(t) exp(lambda t) / lambda
    a_bq, t_h = 1e8, 24.0
    lam = np.log(2) / (get_physical_half_life_h("Lu177") * 3600.0)
    tia_bqs = a_bq * np.exp(lam * t_h * 3600.0) / lam

    res_a, dose_a = compute_roi_dose(
        image=make_image(a_bq, "bq", time_from_injection=t_h),
        labels=make_labels(),
        method="ldm",
        isotope="Lu177",
        labels_map={1: "Left kidney"},
    )
    res_t, dose_t = compute_roi_dose(
        image=make_image(tia_bqs, "bqs"),
        labels=make_labels(),
        method="ldm",
        isotope="Lu177",
        labels_map={1: "Left kidney"},
    )
    assert res_t["Left kidney"]["mean_dose_gy"] == pytest.approx(
        res_a["Left kidney"]["mean_dose_gy"], rel=1e-5
    )
    assert res_t["Left kidney"]["tia_mbq_s"] == pytest.approx(tia_bqs / 1e6, rel=1e-6)
    np.testing.assert_allclose(
        sitk.GetArrayFromImage(dose_t), sitk.GetArrayFromImage(dose_a), rtol=1e-5
    )


def test_tia_input_errors():
    common = dict(labels=make_labels(), isotope="Lu177", labels_map={1: "Left kidney"})

    # Single-time-point methods need an activity image
    with pytest.raises(ValueError, match="cannot be used with a time-integrated"):
        compute_roi_dose(
            image=make_image(1e12, "bqs"),
            method="madsen2018",
            method_params={"Left kidney": {"teff": 45.0}},
            **common,
        )
    # svalue needs a TIA image
    with pytest.raises(ValueError, match="requires a time-integrated"):
        compute_roi_dose(image=make_image(1e8, "bq", time_from_injection=24.0), method="svalue", **common)
    # No dose rate, no activity normalisation from a TIA image
    with pytest.raises(ValueError, match="Dose rate"):
        compute_roi_dose(image=make_image(1e12, "bqs"), method="ldm", dose_rate=True, **common)
    with pytest.raises(ValueError, match="total-activity"):
        compute_roi_dose(image=make_image(1e12, "bqs"), method="svalue", total_activity_mbq=7000.0, **common)


def test_svalue_cli(tmp_path):
    tia_path = tmp_path / "tia.nii.gz"
    labels_path = tmp_path / "labels.nii.gz"
    out_path = tmp_path / "dose.json"
    sitk.WriteImage(make_image(1e12, "bqsml"), str(tia_path))
    sitk.WriteImage(make_labels(), str(labels_path))
    # Sidecars: the TIA unit (no time point needed) and the label names
    (tmp_path / "tia.json").write_text(json.dumps({"modality": "spect", "unit": "bqsml"}))
    (tmp_path / "labels.json").write_text(json.dumps({"1": "Left kidney"}))

    result = CliRunner().invoke(
        main,
        ["dose", str(tia_path), "-l", str(labels_path), "-m", "svalue", "--isotope", "Lu177", "-o", str(out_path)],
    )
    assert result.exit_code == 0, result.output
    res = json.loads(out_path.read_text())
    expected = 1e6 * S_KIDNEY_L * M_KIDNEY_L / 1000.0
    assert res["Left kidney"]["mean_dose_gy"] == pytest.approx(expected, rel=1e-6)
