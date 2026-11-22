"""radot dose on a ROI TIA table (radot integrate output): fit first, ROI level (roadmap E3)."""
import json
import logging

import numpy as np
import pytest
import SimpleITK as sitk
from click.testing import CliRunner

from radot.cli import main
from radot.db import OrganNotFoundError
from radot.dose import compute_roi_dose, compute_roi_dose_file, compute_roi_dose_tia_table, read_results
from radot.integrate import integrate_rois, write_integrate

from test_integrate import _study, _write


@pytest.fixture
def table(tmp_path):
    """TIA table (f2 per ROI) of the synthetic study, patient sex male; and its labels."""
    images, labels = _study(tmp_path, sex="male")
    doc = integrate_rois(images, labels, model="f2", error_model="proportional")
    path = tmp_path / "tia_rois.json"
    write_integrate(doc, path)
    return path, labels, doc


def _tia_image(doc, labels):
    """TIA image (bqs) with the ROI totals of the table, uniform in each ROI."""
    lab = sitk.GetArrayFromImage(sitk.ReadImage(str(labels)))
    arr = np.zeros(lab.shape, np.float64)
    for lid, name in ((1, "Liver"), (2, "Spleen")):
        arr[lab == lid] = doc["rois"][name]["tia_mbq_s"] * 1e6 / np.sum(lab == lid)
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing((2.0, 2.0, 2.0))
    img.SetMetaData("unit", "bqs")
    return img


def test_same_dose_as_svalue_on_a_tia_image(table):
    path, labels, doc = table
    res = compute_roi_dose_tia_table(path, labels)
    ref = compute_roi_dose(image=_tia_image(doc, labels), method="svalue", isotope="Lu177",
                           labels=sitk.ReadImage(str(labels)), labels_map={1: "Liver", 2: "Spleen"},
                           method_params={"phantom": "male"})
    for roi in ("Liver", "Spleen"):
        assert res[roi]["mean_dose_gy"] == pytest.approx(ref[roi]["mean_dose_gy"], rel=1e-9)
        assert res[roi]["roi_mass_g"] == pytest.approx(ref[roi]["roi_mass_g"])
        assert res[roi]["tia_mbq_s"] == pytest.approx(doc["rois"][roi]["tia_mbq_s"])
        assert res[roi]["tia_model"] == "f2"
        se_ratio = doc["rois"][roi]["tia_se_mbq_h"] / doc["rois"][roi]["tia_mbq_h"]
        assert res[roi]["mean_dose_fit_se_gy"] == pytest.approx(res[roi]["mean_dose_gy"] * se_ratio)


def test_ct_mass(table, tmp_path):
    path, labels, _ = table
    water = _write(tmp_path / "ct0.nii.gz", np.zeros((10, 10, 10), np.float32), {"modality": "ct", "unit": "hu"})
    dense = _write(tmp_path / "ct1.nii.gz", np.full((10, 10, 10), 1000.0, np.float32), {"modality": "ct", "unit": "hu"})
    base = compute_roi_dose_tia_table(path, labels)
    r0 = compute_roi_dose_tia_table(path, labels, ct_image=water)
    r1 = compute_roi_dose_tia_table(path, labels, ct_image=dense)
    assert r0["Liver"]["mean_dose_gy"] == pytest.approx(base["Liver"]["mean_dose_gy"])
    assert r1["Liver"]["roi_mass_g"] == pytest.approx(2 * base["Liver"]["roi_mass_g"])
    assert r1["Liver"]["mean_dose_gy"] == pytest.approx(base["Liver"]["mean_dose_gy"] / 2)


# AGENTS rule 5: the phantom from the table (patient sex of the sidecars), else an explicit error
def test_phantom_from_table_or_option(tmp_path, caplog):
    images, labels = _study(tmp_path)  # no sex
    doc = integrate_rois(images, labels, model="f2", error_model="proportional")
    with pytest.raises(ValueError, match="OpenDose phantom not provided"):
        compute_roi_dose_tia_table(doc, labels)
    assert "Liver" in compute_roi_dose_tia_table(doc, labels, method_params={"phantom": "female"})
    doc["sex"] = "female"
    with caplog.at_level(logging.INFO, logger="radot"):
        compute_roi_dose_tia_table(doc, labels)
    assert "OpenDose phantom: female (patient sex recorded" in caplog.text


def test_isotope_option_vs_table(table):
    path, labels, _ = table
    assert compute_roi_dose_tia_table(path, labels, isotope="177Lu")
    with pytest.raises(ValueError, match="Radionuclide mismatch"):
        compute_roi_dose_tia_table(path, labels, isotope="Y90")


@pytest.mark.parametrize("kwargs, err", [
    ({"method": "hanscheid2018"}, "takes method svalue"),
    ({"method": "ldm"}, "takes method svalue"),
    ({"method": "svalue", "total_activity_mbq": 100.0}, "total_activity_mbq"),
    ({"method": "svalue", "method_params": {"Liver": {"teff": 50.0}}}, "teff for Liver"),
    ({"method": "svalue", "cascade": True}, "cascade"),
    ({"method": "svalue", "resample_target": "like_ct"}, "resample 'like_ct'"),
])
def test_unused_options_are_errors(table, kwargs, err):
    path, labels, _ = table
    with pytest.raises(ValueError, match=err):
        compute_roi_dose_file(image_path=path, labels_path=labels, **kwargs)


def test_rois_rules(table, tmp_path, caplog):
    path, labels, _ = table
    with pytest.raises(ValueError, match="not in tia_rois.json: Kidney_L"):
        compute_roi_dose_tia_table(path, labels, rois=["Kidney_L"])
    # labels without the spleen: skipped in the loop, an error when requested
    lab = sitk.GetArrayFromImage(sitk.ReadImage(str(labels)))
    lab[lab == 2] = 0
    only_liver = _write(tmp_path / "liver.nii.gz", lab, {"1": "Liver", "2": "Spleen"})
    with caplog.at_level(logging.WARNING, logger="radot"):
        res = compute_roi_dose_tia_table(path, only_liver)
    assert list(res) == ["Liver"] and "ROI 'Spleen' has no voxel in the labels" in caplog.text
    with pytest.raises(ValueError, match="'Spleen' has no voxel"):
        compute_roi_dose_tia_table(path, only_liver, rois=["Spleen"])


def test_organ_without_opendose_data(tmp_path, caplog):
    images, labels = _study(tmp_path, sex="male")
    (tmp_path / "labels.json").write_text(json.dumps({"1": "Liver", "2": "Liver_met1"}))
    doc = integrate_rois(images, labels, model="f2", error_model="proportional")
    with caplog.at_level(logging.WARNING, logger="radot"):
        res = compute_roi_dose_tia_table(doc, labels)
    assert list(res) == ["Liver"] and "Skipping dosimetry for ROI 'Liver_met1'" in caplog.text
    with pytest.raises(OrganNotFoundError):
        compute_roi_dose_tia_table(doc, labels, rois=["Liver_met1"])


def test_cli_dose_on_tia_table(table, tmp_path):
    path, labels, doc = table
    out = tmp_path / "dose.json"
    result = CliRunner().invoke(main, ["dose", str(path), "-m", "svalue", "-l", str(labels), "-o", str(out)])
    assert result.exit_code == 0, result.output
    rois = read_results(out)
    assert set(rois) == {"Liver", "Spleen"} and rois["Liver"]["tia_model"] == "f2"
    # a dose config (not an integrate file) is still read as a config
    cfg = tmp_path / "cfg.json"
    cfg.write_text(json.dumps({"image": str(path), "labels": str(labels), "method": "svalue",
                               "output": str(tmp_path / "dose2.json")}))
    result = CliRunner().invoke(main, ["dose", str(cfg)])
    assert result.exit_code == 0, result.output
    assert read_results(tmp_path / "dose2.json")["Liver"]["mean_dose_gy"] == pytest.approx(
        rois["Liver"]["mean_dose_gy"])
