import json

import numpy as np
import pytest
import SimpleITK as sitk
from click.testing import CliRunner

from radot.cli import main
from radot.dose import compute_roi_dose

from test_dose_tia import M_KIDNEY_L, S_KIDNEY_L, make_image, make_labels

ROI = "Left kidney"
COMMON = dict(isotope="Lu177", labels_map={1: ROI})


def svalue_map(source_mbq=1.0):
    """Dose-rate map (Gy/s) of the S-value self-dose of source_mbq in the 1 mL ROI voxel."""
    m = make_image(source_mbq * S_KIDNEY_L * M_KIDNEY_L / 1.0 / 1000.0, "gy/s")
    m.EraseMetaData("unit")
    return m


@pytest.mark.parametrize("method", ["hanscheid2018", "madsen2018"])
def test_map_reproduces_svalue_methods(method):
    # A map equal to the S-value dose rate gives the same dose as the S-value method,
    # for any simulated activity (the map is scaled to the image activity)
    image = make_image(5e8, "bq", time_from_injection=24.0)
    params = {ROI: {"teff": 45.0}}
    ref = compute_roi_dose(image=image, labels=make_labels(), method=method, method_params=params, **COMMON)
    for source_mbq in (1.0, 250.0):
        res = compute_roi_dose(
            image=image,
            labels=make_labels(),
            method=method,
            method_params=params,
            dose_rate_map=svalue_map(source_mbq),
            dose_rate_map_activity_mbq=source_mbq,
            **COMMON,
        )
        r = res[ROI]
        assert r["mean_dose_gy"] == pytest.approx(ref[ROI]["mean_dose_gy"], rel=1e-6)
        assert r["roi_teff_h"] == pytest.approx(45.0)
        assert r["mean_dose_rate_gys"] == pytest.approx(500.0 * S_KIDNEY_L * M_KIDNEY_L / 1000.0, rel=1e-6)
        assert r["dose_engine"] == "map"


def test_map_hanscheid2018_default_teff():
    image = make_image(5e8, "bq", time_from_injection=24.0)
    res = compute_roi_dose(
        image=image, labels=make_labels(), method="hanscheid2018",
        dose_rate_map=svalue_map(), dose_rate_map_activity_mbq=1.0, **COMMON,
    )
    r = res[ROI]
    assert r["roi_teff_h"] == pytest.approx(48.0)
    assert r["mean_dose_gy"] == pytest.approx(r["mean_dose_rate_gys"] * 2 * 24.0 * 3600.0 / np.log(2))


def test_map_tia_reproduces_svalue():
    tia = make_image(1e12, "bqs")
    ref = compute_roi_dose(image=tia, labels=make_labels(), method="svalue", **COMMON)
    for method in ("svalue", "ldm"):
        res = compute_roi_dose(
            image=tia, labels=make_labels(), method=method,
            dose_rate_map=svalue_map(), dose_rate_map_activity_mbq=1.0, **COMMON,
        )
        if method == "ldm":
            res, dose_image = res
            assert sitk.GetArrayFromImage(dose_image).max() == pytest.approx(ref[ROI]["mean_dose_gy"], rel=1e-6)
        assert res[ROI]["mean_dose_gy"] == pytest.approx(ref[ROI]["mean_dose_gy"], rel=1e-6)
        assert res[ROI]["tia_mbq_s"] == pytest.approx(1e6)


def test_map_reproduces_ldm():
    # A map equal to the LDM dose rate (Gy/h, absolute) gives the LDM dose and dose rate
    image = make_image(5e8, "bq", time_from_injection=24.0)
    ref, _ = compute_roi_dose(image=image, labels=make_labels(), method="ldm", **COMMON)
    _, rate_image = compute_roi_dose(image=image, labels=make_labels(), method="ldm", dose_rate=True, **COMMON)
    res, dose_image = compute_roi_dose(
        image=image, labels=make_labels(), method="ldm", dose_rate_map=rate_image, **COMMON
    )
    assert res[ROI]["mean_dose_gy"] == pytest.approx(ref[ROI]["mean_dose_gy"], rel=1e-5)
    rate, _ = compute_roi_dose(
        image=image, labels=make_labels(), method="ldm", dose_rate=True, dose_rate_map=rate_image, **COMMON
    )
    assert rate[ROI]["mean_dose_rate_gyh"] == pytest.approx(sitk.GetArrayFromImage(rate_image).max(), rel=1e-6)
    # Without labels: the dose image only
    dose_only = compute_roi_dose(image=image, method="ldm", dose_rate_map=rate_image, isotope="Lu177")
    np.testing.assert_allclose(sitk.GetArrayFromImage(dose_only), sitk.GetArrayFromImage(dose_image))


def test_map_errors():
    image = make_image(5e8, "bq", time_from_injection=24.0)
    common = dict(labels=make_labels(), **COMMON)
    with pytest.raises(ValueError, match="cannot be used with a dose-rate map"):
        compute_roi_dose(image=image, method="hanscheid2017", dose_rate_map=svalue_map(), **common)
    with pytest.raises(ValueError, match="activity of the dose-rate map is required"):
        compute_roi_dose(image=make_image(1e12, "bqs"), method="svalue", dose_rate_map=svalue_map(), **common)
    with pytest.raises(ValueError, match="teff"):
        compute_roi_dose(image=image, method="madsen2018", dose_rate_map=svalue_map(), dose_rate_map_activity_mbq=1.0, **common)
    bad_grid = svalue_map()
    bad_grid.SetSpacing((8.0, 8.0, 8.0))
    with pytest.raises(ValueError, match="grid"):
        compute_roi_dose(image=image, method="hanscheid2018", dose_rate_map=bad_grid, **common)
    bad_unit = svalue_map()
    bad_unit.SetMetaData("unit", "gy")
    with pytest.raises(ValueError, match="unit"):
        compute_roi_dose(image=image, method="hanscheid2018", dose_rate_map=bad_unit, **common)


def test_map_cli(tmp_path):
    tia_path, labels_path = tmp_path / "tia.nii.gz", tmp_path / "labels.nii.gz"
    map_path, out_path = tmp_path / "gate_dose.mhd", tmp_path / "dose.json"
    sitk.WriteImage(make_image(1e12, "bqsml"), str(tia_path))
    sitk.WriteImage(make_labels(), str(labels_path))
    sitk.WriteImage(svalue_map(), str(map_path))
    (tmp_path / "tia.json").write_text(json.dumps({"modality": "spect", "unit": "bqsml"}))
    (tmp_path / "labels.json").write_text(json.dumps({"1": ROI}))

    args = ["dose", str(tia_path), "-l", str(labels_path), "-m", "svalue", "--isotope", "Lu177",
            "--dose-rate-map", str(map_path), "-o", str(out_path)]
    # The map activity is required with a TIA image: from the option, or from the map sidecar
    result = CliRunner().invoke(main, args)
    assert result.exit_code != 0
    result = CliRunner().invoke(main, args + ["--dose-rate-map-activity", "1"])
    assert result.exit_code == 0, result.output
    expected = 1e6 * S_KIDNEY_L * M_KIDNEY_L / 1000.0
    assert json.loads(out_path.read_text())[ROI]["mean_dose_gy"] == pytest.approx(expected, rel=1e-6)

    (tmp_path / "gate_dose.json").write_text(json.dumps({"unit": "gy/s", "source_activity_mbq": 1.0}))
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, result.output
    assert json.loads(out_path.read_text())[ROI]["mean_dose_gy"] == pytest.approx(expected, rel=1e-6)
