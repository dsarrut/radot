"""Audit B9: an option the chosen computation does not use is an error, never silently ignored."""
import json

import numpy as np
import pytest
import SimpleITK as sitk
from click.testing import CliRunner

from radot.cli import main
from radot.dose import compute_roi_dose, compute_roi_dose_file
from radot.dose.compute import check_unused_options

TEFF = {"Liver": {"teff": 45.0}}


@pytest.mark.parametrize("method, kw, match", [
    # S-value methods: no voxel options
    ("hanscheid2018", {"cascade": True}, "cascade"),
    ("hanscheid2018", {"hu_threshold": -800.0}, "hu_threshold"),
    ("madsen2018", {"dose_rate": True}, "dose_rate"),
    ("svalue", {"method_params": TEFF}, "teff for Liver"),
    # ldm: no Teff (physical decay)
    ("ldm", {"method_params": TEFF}, "teff for Liver"),
    # dose-rate map: the map replaces the local deposition
    ("ldm", {"dose_rate_map": True, "cascade": True}, "cascade"),
    ("hanscheid2018", {"dose_rate_map": True, "hu_threshold": -500.0}, "hu_threshold"),
    # multiple time points
    ("jackson2013", {"multi_tp": True, "total_activity_mbq": 7400.0}, "total_activity_mbq"),
    ("jackson2013", {"multi_tp": True, "resample": "like_ct"}, "resample 'like_ct'"),
    ("jackson2013", {"multi_tp": True, "method_params": TEFF}, "teff"),
    ("jackson2013", {"multi_tp": True, "dose_rate": True}, "dose_rate"),
])
def test_unused_option_errors(method, kw, match):
    with pytest.raises(ValueError, match=f"Not used with .*{match}"):
        check_unused_options(method, **kw)


@pytest.mark.parametrize("method, kw", [
    ("ldm", {"cascade": True, "hu_threshold": -500.0, "dose_rate": True, "total_activity_mbq": 1.0}),
    ("hanscheid2018", {"method_params": {**TEFF, "phantom": "male"}, "total_activity_mbq": 1.0}),
    ("madsen2018", {"dose_rate_map": True, "dose_rate": True, "method_params": TEFF}),
    ("hanscheid2017", {"method_params": TEFF, "cascade": False, "dose_rate": False}),
    ("jackson2013", {"multi_tp": True, "method_params": {"phantom": "male", "times_h": [1, 2, 3]}}),
    ("jackson2013", {"multi_tp": True, "resample": "like_labels"}),
])
def test_used_options_accepted(method, kw):
    check_unused_options(method, **kw)


def _study(tmp_path):
    image = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    image[2, 2, 2] = 1000.0
    image.SetMetaData("unit", "bq")
    image.SetMetaData("time_from_injection", "24.0")
    labels = sitk.Image([5, 5, 5], sitk.sitkUInt8)
    labels[2, 2, 2] = 1
    return image, labels


def test_api_hanscheid_with_cascade_errors(tmp_path):
    image, labels = _study(tmp_path)
    with pytest.raises(ValueError, match="method 'hanscheid2018': cascade"):
        compute_roi_dose(image, "hanscheid2018", "Lu177", labels=labels, labels_map={1: "Liver"},
                         cascade=True, method_params={"phantom": "male"})


def test_api_ldm_with_teff_params_errors(tmp_path):
    image, labels = _study(tmp_path)
    with pytest.raises(ValueError, match="method 'ldm': teff for Liver"):
        compute_roi_dose(image, "ldm", "Lu177", labels=labels, labels_map={1: "Liver"},
                         teff_params={"Liver": 45.0})


def test_api_ldm_default_hu_threshold_still_applied(tmp_path):
    """hu_threshold not given: ldm still zeroes the air voxels (-800 HU)."""
    image, labels = _study(tmp_path)
    ct = sitk.Image([5, 5, 5], sitk.sitkFloat32) - 1000.0
    res, dose = compute_roi_dose(image, "ldm", "Lu177", labels=labels, labels_map={1: "Liver"}, ct_image=ct)
    assert res["Liver"]["mean_dose_gy"] == 0.0


def test_api_mtp_with_total_activity_errors(tmp_path):
    image, labels = _study(tmp_path)
    with pytest.raises(ValueError, match="multiple time points: total_activity_mbq"):
        compute_roi_dose([image, image, image], "jackson2013", "Lu177", labels=labels,
                         total_activity_mbq=7400.0)


def test_file_mtp_config_with_hu_threshold_errors(tmp_path):
    cfg = {"method": "jackson2013", "radionuclide": "Lu177", "labels": "labels.nii.gz",
           "timepoints": ["a.nii.gz", "b.nii.gz", "c.nii.gz"], "hu_threshold": -500}
    with pytest.raises(ValueError, match="multiple time points: hu_threshold"):
        compute_roi_dose_file(config=cfg)


def test_cli_hanscheid_with_hu_threshold_errors(tmp_path):
    image, labels = _study(tmp_path)
    sitk.WriteImage(image, str(tmp_path / "spect.nii.gz"))
    (tmp_path / "spect.json").write_text(json.dumps({"unit": "bq", "time_from_injection": 24.0, "sex": "male"}))
    sitk.WriteImage(labels, str(tmp_path / "labels.nii.gz"))
    (tmp_path / "labels.json").write_text(json.dumps({"1": "Liver"}))
    result = CliRunner().invoke(main, ["dose", str(tmp_path / "spect.nii.gz"), "-l", str(tmp_path / "labels.nii.gz"),
                                       "-m", "hanscheid2018", "--isotope", "Lu177", "--hu-threshold", "-500"])
    assert result.exit_code != 0
    assert "hu_threshold (only ldm" in result.output
