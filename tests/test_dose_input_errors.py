"""Inputs that were silently replaced by a default (audit A9-A12): now explicit errors."""
import json

import numpy as np
import pytest
import SimpleITK as sitk

from radot.dose import compute_roi_dose, compute_roi_dose_mtp
from radot.dose.config import load_dose_config

DATES = ["2024-01-01 15:00:00", "2024-01-02 10:00:00", "2024-01-06 10:00:00"]


def _write(path, arr, sidecar=None):
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing((2.0, 2.0, 2.0))
    sitk.WriteImage(img, str(path))
    if sidecar is not None:
        path.with_name(path.name.split(".")[0] + ".json").write_text(json.dumps(sidecar))
    return path


def _labels_arr(spleen=True):
    a = np.zeros((10, 10, 10), dtype=np.uint8)
    a[2:5, 2:5, 2:5] = 1
    if spleen:
        a[6:8, 6:8, 6:8] = 2
    return a


def _mtp_study(tmp_path, unit="bqml", labels_sidecar={"1": "Liver", "2": "Spleen"}):
    labels = _write(tmp_path / "labels.nii.gz", _labels_arr(), labels_sidecar)
    spects = []
    for i, date in enumerate(DATES):
        sc = {"modality": "spect", "injection_date": "2024-01-01 10:00:00",
              "acquisition_date": date, "activity": 5000.0, "sex": "male"}
        if unit:
            sc["unit"] = unit
        arr = np.where(_labels_arr() > 0, 1e4 * np.exp(-0.03 * (i + 1) * 24), 0).astype(np.float32)
        spects.append(_write(tmp_path / f"spect_tp{i+1}.nii.gz", arr, sc))
    return spects, labels


# A9: no unit at a time point
def test_mtp_unit_missing_errors(tmp_path):
    spects, labels = _mtp_study(tmp_path, unit=None)
    with pytest.raises(ValueError, match="Time point 1: image unit"):
        compute_roi_dose_mtp(images=spects, labels=labels, isotope="Lu177",
                             method_params={"phantom": "male"})


# A10: ROI absent at a time point
def _per_tp_labels(tmp_path):
    """Spleen absent from the labels of time point 2."""
    out = []
    for i in range(3):
        out.append(_write(tmp_path / f"labels_tp{i+1}.nii.gz", _labels_arr(spleen=(i != 1)),
                          {"1": "Liver", "2": "Spleen"}))
    return out


def test_mtp_roi_absent_at_a_time_point_requested_errors(tmp_path):
    spects, _ = _mtp_study(tmp_path)
    with pytest.raises(ValueError, match="'Spleen' is absent .* time point\\(s\\) 2"):
        compute_roi_dose_mtp(images=spects, labels=_per_tp_labels(tmp_path), isotope="Lu177",
                             rois=["Liver", "Spleen"], method_params={"phantom": "male"})


def test_mtp_roi_absent_at_a_time_point_loop_skips_with_warning(tmp_path, caplog):
    spects, _ = _mtp_study(tmp_path)
    res = compute_roi_dose_mtp(images=spects, labels=_per_tp_labels(tmp_path), isotope="Lu177",
                               method_params={"phantom": "male"})
    assert "Liver" in res and "Spleen" not in res
    assert "'Spleen' is absent" in caplog.text


# A11: config Teff not a number, or not ROI:hours
@pytest.mark.parametrize("cfg, match", [
    ({"teff": {"Liver": "abc"}}, "Teff of ROI 'Liver' must be a number"),
    ({"rois": {"Liver": "abc"}}, "Teff of ROI 'Liver' must be a number"),
    ({"teff": ["Liver=45"]}, "'Liver=45' must be in ROI:hours"),
    ({"teff": ["Liver:4:5"]}, "'Liver:4:5' must be in ROI:hours"),
    ({"teff": ["Liver:abc"]}, "Teff of ROI 'Liver' must be a number"),
    ({"teff": 45.0}, "must be a dict"),
    ({"teff": [{"Liver": 45}]}, "must be a 'ROI:hours' string"),
])
def test_config_bad_teff_errors(cfg, match):
    with pytest.raises(ValueError, match=match):
        load_dose_config({"method": "hanscheid2018", "isotope": "Lu177", **cfg})


def test_config_teff_forms_accepted():
    for cfg in ({"teff": {"Liver": 45}}, {"teff": ["Liver:45"]}, {"teff": "Liver:45"},
                {"rois": {"Liver": 45}}):
        mp = load_dose_config({"method": "hanscheid2018", "isotope": "Lu177", **cfg})["method_params"]
        assert mp["Liver"]["teff"] == 45.0


# A12: labels sidecar missing, unreadable, empty
def _single(tmp_path, labels_sidecar_text):
    labels = _write(tmp_path / "labels.nii.gz", _labels_arr())
    if labels_sidecar_text is not None:
        (tmp_path / "labels.json").write_text(labels_sidecar_text)
    image = sitk.GetImageFromArray(np.ones((10, 10, 10), dtype=np.float32))
    image.SetSpacing((2.0, 2.0, 2.0))
    image.SetMetaData("time_from_injection", "24.0")
    image.SetMetaData("unit", "bq")
    return image, labels


@pytest.mark.parametrize("text, match", [
    (None, "has no sidecar 'labels.json'"),
    ("{ broken", "cannot be read"),
    ('{"Liver": 1}', "cannot be read"),
    ("{}", "has no ROI"),
])
@pytest.mark.parametrize("method", ["hanscheid2017", "ldm"])
def test_labels_sidecar_errors(tmp_path, text, match, method):
    image, labels = _single(tmp_path, text)
    with pytest.raises(ValueError, match=match):
        compute_roi_dose(image=image, labels=sitk.ReadImage(str(labels)), labels_path=labels,
                         method=method, isotope="Lu177")


def test_mtp_labels_sidecar_missing_errors(tmp_path):
    spects, labels = _mtp_study(tmp_path, labels_sidecar=None)
    with pytest.raises(ValueError, match="has no sidecar"):
        compute_roi_dose_mtp(images=spects, labels=labels, isotope="Lu177",
                             method_params={"phantom": "male"})


# B11: a ROI requested with rois but with no voxel is an error (was: absent from the results)
def _with_empty_spleen(tmp_path):
    """Labels sidecar names Liver (1) and Spleen (2); only Liver has voxels."""
    labels = _write(tmp_path / "labels.nii.gz", _labels_arr(spleen=False), {"1": "Liver", "2": "Spleen"})
    image = sitk.GetImageFromArray(np.ones((10, 10, 10), dtype=np.float32))
    image.SetSpacing((2.0, 2.0, 2.0))
    image.SetMetaData("time_from_injection", "24.0")
    image.SetMetaData("unit", "bq")
    return image, labels


@pytest.mark.parametrize("method, extra", [
    ("hanscheid2017", {}),
    ("hanscheid2018", {"method_params": {"phantom": "male"}}),
    ("ldm", {}),
])
def test_requested_roi_without_voxel_errors(tmp_path, method, extra):
    image, labels = _with_empty_spleen(tmp_path)
    with pytest.raises(ValueError, match="no voxel in the label image: Spleen"):
        compute_roi_dose(image=image, labels=sitk.ReadImage(str(labels)), labels_path=labels,
                         method=method, isotope="Lu177", rois=["Liver", "Spleen"], **extra)


def test_requested_roi_without_voxel_rate_map_errors(tmp_path):
    image, labels = _with_empty_spleen(tmp_path)
    rate = sitk.GetImageFromArray(np.full((10, 10, 10), 1e-9, dtype=np.float32))
    rate.SetSpacing((2.0, 2.0, 2.0))
    rate.SetMetaData("unit", "gy/s")
    with pytest.raises(ValueError, match="no voxel in the label image: Spleen"):
        compute_roi_dose(image=image, labels=sitk.ReadImage(str(labels)), labels_path=labels,
                         method="ldm", isotope="Lu177", rois=["Spleen"], dose_rate_map=rate)


def test_requested_roi_not_in_given_labels_map_errors(tmp_path):
    """API labels_map given directly (no sidecar check): an unknown requested ROI is an error too."""
    image, labels = _with_empty_spleen(tmp_path)
    with pytest.raises(ValueError, match="no voxel in the label image: Kidney_L"):
        compute_roi_dose(image=image, labels=sitk.ReadImage(str(labels)), labels_map={1: "Liver"},
                         method="hanscheid2017", isotope="Lu177", rois=["Kidney_L"])


def test_all_rois_loop_skips_empty_label(tmp_path):
    """Without rois, a named label with no voxel is simply not in the results (nothing requested)."""
    image, labels = _with_empty_spleen(tmp_path)
    res = compute_roi_dose(image=image, labels=sitk.ReadImage(str(labels)), labels_path=labels,
                           method="hanscheid2017", isotope="Lu177")
    assert set(res) == {"Liver"}


def test_rois_without_labels_errors(tmp_path):
    image, _ = _with_empty_spleen(tmp_path)
    with pytest.raises(ValueError, match="rois needs a label image"):
        compute_roi_dose(image=image, method="ldm", isotope="Lu177", rois=["Liver"])


def test_mtp_requested_roi_absent_everywhere_errors(tmp_path):
    spects, labels = _mtp_study(tmp_path, labels_sidecar={"1": "Liver", "2": "Spleen", "3": "Kidney_L"})
    with pytest.raises(ValueError, match="'Kidney_L' has no voxel in the labels at any time point"):
        compute_roi_dose_mtp(images=spects, labels=labels, isotope="Lu177", rois=["Liver", "Kidney_L"],
                             method_params={"phantom": "male"})


# B12: no default radionuclide
def test_mtp_no_isotope_errors(tmp_path):
    spects, labels = _mtp_study(tmp_path)
    with pytest.raises(ValueError, match="Radionuclide \\(isotope\\) is required"):
        compute_roi_dose_mtp(images=spects, labels=labels, method_params={"phantom": "male"})
