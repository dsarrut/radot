"""OpenDose phantom: --phantom (or config), else the patient sex in the sidecar, else an explicit error."""

import json

import pytest
import SimpleITK as sitk
from click.testing import CliRunner

from radot.cli import main


def _study(tmp_path, sex=None):
    tmp_path.mkdir(parents=True, exist_ok=True)
    image = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    image.SetSpacing((10.0, 10.0, 10.0))
    image[2, 2, 2] = 1e8
    labels = sitk.Image([5, 5, 5], sitk.sitkUInt16)
    labels.SetSpacing((10.0, 10.0, 10.0))
    labels[2, 2, 2] = 1
    sitk.WriteImage(image, str(tmp_path / "spect.nii.gz"))
    sitk.WriteImage(labels, str(tmp_path / "labels.nii.gz"))
    meta = {"modality": "spect", "unit": "bq", "injection_date": "2026-01-10 10:00:00",
            "acquisition_date": "2026-01-11 10:00:00"}
    if sex is not None:
        meta["sex"] = sex
    (tmp_path / "spect.json").write_text(json.dumps(meta))
    (tmp_path / "labels.json").write_text(json.dumps({"1": "Liver"}))
    return tmp_path


def _dose(study, *args):
    out = study / "dose.json"
    result = CliRunner().invoke(
        main, ["dose", str(study / "spect.nii.gz"), "-l", str(study / "labels.nii.gz"), "--isotope", "Lu177",
               "-o", str(out), *args]
    )
    return result, (json.loads(out.read_text()) if out.exists() else None)


def test_phantom_from_sidecar_sex(tmp_path):
    # only in the sidecar: female phantom used (mocked S-values differ between phantoms)
    _, female = _dose(_study(tmp_path / "f", sex="female"), "-m", "hanscheid2018")
    _, explicit = _dose(_study(tmp_path / "e"), "-m", "hanscheid2018", "--phantom", "female")
    _, male = _dose(_study(tmp_path / "m"), "-m", "hanscheid2018", "--phantom", "male")
    assert female["Liver"]["mean_dose_gy"] == pytest.approx(explicit["Liver"]["mean_dose_gy"])
    assert female["Liver"]["mean_dose_gy"] != pytest.approx(male["Liver"]["mean_dose_gy"])


def test_phantom_option_wins_over_sidecar(tmp_path):
    _, opt = _dose(_study(tmp_path / "a", sex="female"), "-m", "hanscheid2018", "--phantom", "male")
    _, male = _dose(_study(tmp_path / "b", sex="male"), "-m", "hanscheid2018")
    assert opt["Liver"]["mean_dose_gy"] == pytest.approx(male["Liver"]["mean_dose_gy"])


@pytest.mark.parametrize("method", ["hanscheid2018", "madsen2018"])
def test_phantom_missing_is_an_error(tmp_path, method):
    result, _ = _dose(_study(tmp_path), "-m", method, "--teff", "Liver:50")
    assert result.exit_code != 0
    assert "phantom not provided" in result.output and "--sex" in result.output


@pytest.mark.parametrize("method", ["hanscheid2017", "ldm"])
def test_no_phantom_needed_without_s_values(tmp_path, method):
    result, res = _dose(_study(tmp_path), "-m", method)
    assert result.exit_code == 0, result.output
    assert res["Liver"]["mean_dose_gy"] > 0
