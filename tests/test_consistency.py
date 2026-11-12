"""Audit 100, part D: one rule per behaviour (density, unit case, empty ROI list, grid mismatch,
radionuclide from the sidecar)."""
import json

import numpy as np
import pytest
import SimpleITK as sitk

from radot.dose import compute_roi_dose, compute_roi_dose_mtp
from radot.stats import DENSITY_FLOOR_G_ML, compute_stats, ct_density_g_ml


# D1: one CT density rule, voxel-wise, floored

def test_ct_density_rule():
    assert ct_density_g_ml(0.0) == pytest.approx(1.0)
    assert ct_density_g_ml(1000.0) == pytest.approx(2.0)
    assert np.allclose(ct_density_g_ml([-1000.0, -1024.0, -3000.0]), DENSITY_FLOOR_G_ML)


def _roi_below_air():
    """One ROI of 64 voxels (1 mL each): 32 at 0 HU, 32 at -1024 HU (CT padding)."""
    ct = np.zeros((4, 4, 4), dtype=np.float32)
    ct[:2] = -1024.0
    ct_img = sitk.GetImageFromArray(ct)
    ct_img.SetSpacing((10.0, 10.0, 10.0))
    labels = sitk.GetImageFromArray(np.ones((4, 4, 4), dtype=np.uint8))
    labels.SetSpacing((10.0, 10.0, 10.0))
    return ct_img, labels


EXPECTED_MASS_G = 32 * 1.0 + 32 * DENSITY_FLOOR_G_ML  # sum of voxel masses


def test_stats_ct_mass_is_sum_of_voxel_masses():
    """Was: volume x max(0, 1 + mean HU / 1000) = 64 x 0.488."""
    ct, labels = _roi_below_air()
    res = compute_stats(labels, ct, image_type="ct", labels_map={1: "Liver"})
    assert res["1"]["mass_g"] == pytest.approx(EXPECTED_MASS_G)


def test_dose_roi_mass_same_rule_svalue_and_ldm():
    """S-value methods (was: clamped at 0) and ldm give the same ROI mass."""
    ct, labels = _roi_below_air()
    image = sitk.GetImageFromArray(np.ones((4, 4, 4), dtype=np.float32))
    image.SetSpacing((10.0, 10.0, 10.0))
    image.SetMetaData("unit", "bq")
    image.SetMetaData("time_from_injection", "24.0")
    h = compute_roi_dose(image, "hanscheid2017", "Lu177", labels=labels, labels_map={1: "Liver"}, ct_image=ct)
    ldm, _ = compute_roi_dose(image, "ldm", "Lu177", labels=labels, labels_map={1: "Liver"}, ct_image=ct)
    assert h["Liver"]["roi_mass_g"] == pytest.approx(EXPECTED_MASS_G)
    assert ldm["Liver"]["roi_mass_g"] == pytest.approx(EXPECTED_MASS_G)


# D2: sidecar unit not lower case

def test_stats_sidecar_unit_upper_case(tmp_path):
    image = sitk.GetImageFromArray(np.full((4, 4, 4), 2.0, dtype=np.float32))
    path = tmp_path / "s.nii.gz"
    sitk.WriteImage(image, str(path))
    (tmp_path / "s.json").write_text(json.dumps({"modality": "spect", "unit": "BQ"}))
    labels = sitk.GetImageFromArray(np.ones((4, 4, 4), dtype=np.uint8))
    res = compute_stats(labels, image, unit="bq", image_path=path, labels_map={1: "Liver"})
    assert res["1"]["sum"] == pytest.approx(128.0)  # no conversion attempted from 'BQ' to 'bq'


# D3: an empty ROI list is no filter (all the ROIs), on every path

def _single():
    image = sitk.GetImageFromArray(np.ones((4, 4, 4), dtype=np.float32))
    image.SetMetaData("unit", "bq")
    image.SetMetaData("time_from_injection", "24.0")
    arr = np.zeros((4, 4, 4), dtype=np.uint8)
    arr[:2] = 1
    arr[2:] = 2
    return image, sitk.GetImageFromArray(arr), {1: "Liver", 2: "Spleen"}


def test_empty_rois_all_rois_svalue_and_ldm():
    image, labels, names = _single()
    h = compute_roi_dose(image, "hanscheid2017", "Lu177", labels=labels, labels_map=names, rois=[])
    ldm, _ = compute_roi_dose(image, "ldm", "Lu177", labels=labels, labels_map=names, rois=[])
    assert set(h) == set(ldm) == {"Liver", "Spleen"}


def test_empty_rois_all_rois_stats():
    image, labels, names = _single()
    res = compute_stats(labels, image, image_type="spect", unit="bq", labels_map=names, rois=[])
    assert {r["roi"] for r in res.values()} >= {"Liver", "Spleen"}


# D4: multiple time points, grid mismatch: error unless resample like_labels

def _mtp(tmp_path, spect_spacing):
    lab = np.zeros((10, 10, 10), dtype=np.uint8)
    lab[2:6, 2:6, 2:6] = 1
    labels = sitk.GetImageFromArray(lab)
    labels.SetSpacing((2.0, 2.0, 2.0))
    lp = tmp_path / "labels.nii.gz"
    sitk.WriteImage(labels, str(lp))
    (tmp_path / "labels.json").write_text(json.dumps({"1": "Liver"}))
    paths = []
    n = int(round(20 / spect_spacing))
    for i, (t, acq) in enumerate(zip([5.0, 24.0, 120.0],
                                     ["2024-01-01 15:00:00", "2024-01-02 10:00:00", "2024-01-06 10:00:00"])):
        img = sitk.GetImageFromArray(np.full((n, n, n), 1e4 * np.exp(-0.03 * t), dtype=np.float32))
        img.SetSpacing((spect_spacing,) * 3)
        img.SetOrigin(((spect_spacing - 2.0) / 2,) * 3)
        p = tmp_path / f"s{i}.nii.gz"
        sitk.WriteImage(img, str(p))
        (tmp_path / f"s{i}.json").write_text(json.dumps({
            "unit": "bqml", "modality": "spect", "sex": "male", "radionuclide": "Lu177",
            "injection_date": "2024-01-01 10:00:00", "acquisition_date": acq}))
        paths.append(p)
    return paths, lp


def test_mtp_grid_mismatch_errors_without_resample(tmp_path):
    paths, labels = _mtp(tmp_path, 4.0)
    with pytest.raises(ValueError, match="Time point 1: the image grid does not match its labels grid"):
        compute_roi_dose_mtp(images=paths, labels=labels)


def test_mtp_grid_mismatch_resample_like_labels(tmp_path):
    paths, labels = _mtp(tmp_path, 4.0)
    res = compute_roi_dose_mtp(images=paths, labels=labels, resample="like_labels")
    assert res["Liver"]["mean_dose_gy"] > 0


def test_mtp_resample_other_target_errors(tmp_path):
    paths, labels = _mtp(tmp_path, 2.0)
    with pytest.raises(ValueError, match="only 'like_labels'"):
        compute_roi_dose_mtp(images=paths, labels=labels, resample="like_ct")


# D5: radionuclide from the sidecar

def test_isotope_from_sidecar_single(tmp_path):
    image, labels, names = _single()
    image.SetMetaData("radionuclide", "Lu177")
    res = compute_roi_dose(image, "hanscheid2017", None, labels=labels, labels_map=names)
    assert set(res) == {"Liver", "Spleen"}


def test_isotope_option_contradicts_sidecar_errors():
    image, labels, names = _single()
    image.SetMetaData("radionuclide", "Lu177")
    with pytest.raises(ValueError, match="Radionuclide mismatch: the option says 'Y90'"):
        compute_roi_dose(image, "hanscheid2018", "Y90", labels=labels, labels_map=names)


def test_isotope_same_with_other_spelling_ok():
    image, labels, names = _single()
    image.SetMetaData("radionuclide", "Lu177")
    res = compute_roi_dose(image, "hanscheid2017", "lu-177", labels=labels, labels_map=names)
    assert "Liver" in res


def test_isotope_missing_everywhere_errors():
    image, labels, names = _single()
    with pytest.raises(ValueError, match="Radionuclide not given"):
        compute_roi_dose(image, "hanscheid2017", None, labels=labels, labels_map=names)


def test_mtp_isotope_from_sidecars_and_disagreement(tmp_path):
    paths, labels = _mtp(tmp_path, 2.0)
    res = compute_roi_dose_mtp(images=paths, labels=labels)  # Lu177 from the 3 sidecars
    assert res["Liver"]["mean_dose_gy"] > 0
    sc = json.loads((tmp_path / "s2.json").read_text())
    sc["radionuclide"] = "Y90"
    (tmp_path / "s2.json").write_text(json.dumps(sc))
    with pytest.raises(ValueError, match="different radionuclides"):
        compute_roi_dose_mtp(images=paths, labels=labels)


def test_cli_dose_isotope_from_sidecar(tmp_path):
    from click.testing import CliRunner
    from radot.cli import main

    image, labels, _ = _single()
    sitk.WriteImage(image, str(tmp_path / "spect.nii.gz"))
    (tmp_path / "spect.json").write_text(json.dumps({"unit": "bq", "time_from_injection": 24.0,
                                                     "radionuclide": "Lu177"}))
    sitk.WriteImage(labels, str(tmp_path / "labels.nii.gz"))
    (tmp_path / "labels.json").write_text(json.dumps({"1": "Liver", "2": "Spleen"}))
    result = CliRunner().invoke(main, ["dose", str(tmp_path / "spect.nii.gz"), "-l", str(tmp_path / "labels.nii.gz"),
                                       "-m", "hanscheid2017", "-o", str(tmp_path / "out.json")])
    assert result.exit_code == 0, result.output
    assert set(json.loads((tmp_path / "out.json").read_text())) >= {"Liver", "Spleen"}


def test_mtp_config_resample_like_labels(tmp_path):
    from radot.dose import compute_roi_dose_file

    paths, labels = _mtp(tmp_path, 4.0)
    cfg = {"method": "jackson2013", "labels": str(labels), "timepoints": [str(p) for p in paths]}
    with pytest.raises(ValueError, match="grid does not match"):
        compute_roi_dose_file(config=cfg)
    res = compute_roi_dose_file(config={**cfg, "resample": "like_labels"})
    assert res["Liver"]["mean_dose_gy"] > 0
