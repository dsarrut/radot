import json
from pathlib import Path
import pytest
import SimpleITK as sitk
from click.testing import CliRunner

from radot.cli import main
from radot.dose import load_dose_config, compute_roi_dose_file


@pytest.fixture
def dummy_study(tmp_path):
    """Create dummy activity, labels, and CT images for testing."""
    # Activity image
    image = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    image.SetSpacing((2.0, 2.0, 2.0))
    image.SetMetaData("time_from_injection", "24.0")
    image.SetMetaData("unit", "Bq")
    image[2, 2, 2] = 1000.0
    image_path = tmp_path / "spect.nii.gz"
    sitk.WriteImage(image, str(image_path))

    # Write a companion spect.json sidecar to define metadata
    with open(tmp_path / "spect.json", "w") as f:
        json.dump({"unit": "Bq", "time_from_injection": "24.0"}, f)

    # Labels image
    labels = sitk.Image([5, 5, 5], sitk.sitkUInt16)
    labels.SetSpacing((2.0, 2.0, 2.0))
    labels[2, 2, 2] = 1
    labels_path = tmp_path / "labels.nii.gz"
    sitk.WriteImage(labels, str(labels_path))

    # Sidecar JSON for labels
    labels_json_path = tmp_path / "labels.json"
    with open(labels_json_path, "w") as f:
        json.dump({"1": "Liver"}, f)

    # Companion CT
    ct = sitk.Image([5, 5, 5], sitk.sitkInt16)
    ct.SetSpacing((2.0, 2.0, 2.0))
    ct_path = tmp_path / "ct.nii.gz"
    sitk.WriteImage(ct, str(ct_path))

    return {
        "dir": tmp_path,
        "image_path": image_path,
        "labels_path": labels_path,
        "ct_path": ct_path,
    }


def test_load_dose_config_relative_paths(dummy_study):
    study_dir = dummy_study["dir"]
    config_file = study_dir / "dose_config.json"
    config_data = {
        "image": "spect.nii.gz",
        "labels": "labels.nii.gz",
        "ct": "ct.nii.gz",
        "method": "hanscheid2017",
        "radionuclide": "Lu177",
        "total_activity": 250.0,
        "rois": ["Liver"],
        "teff": {"Liver": 45.2},
        "output": "dose_out.json",
    }
    with open(config_file, "w") as f:
        json.dump(config_data, f)

    loaded = load_dose_config(config_file)
    assert loaded["image_path"] == (study_dir / "spect.nii.gz").resolve()
    assert loaded["labels_path"] == (study_dir / "labels.nii.gz").resolve()
    assert loaded["ct_path"] == (study_dir / "ct.nii.gz").resolve()
    assert loaded["output_path"] == (study_dir / "dose_out.json").resolve()
    assert loaded["method"] == "hanscheid2017"
    assert loaded["isotope"] == "Lu177"
    assert loaded["total_activity_mbq"] == 250.0
    assert loaded["rois"] == ["Liver"]
    assert "method_params" in loaded
    assert loaded["method_params"]["Liver"]["teff"] == 45.2


def test_compute_roi_dose_file_with_config_dict(dummy_study):
    cfg = {
        "image_path": dummy_study["image_path"],
        "labels_path": dummy_study["labels_path"],
        "ct_path": dummy_study["ct_path"],
        "method": "hanscheid2017",
        "isotope": "Lu177",
        "total_activity_mbq": 250.0,
        "rois": ["Liver"],
        "teff": {"Liver": 45.2},
    }
    results = compute_roi_dose_file(config=cfg)
    assert "Liver" in results
    assert results["Liver"]["roi_teff_h"] == 45.2
    assert results["Liver"]["mean_dose_gy"] > 0


def test_compute_roi_dose_file_with_json_file(dummy_study):
    study_dir = dummy_study["dir"]
    config_file = study_dir / "study.json"
    output_json = study_dir / "output.json"
    cfg = {
        "image": "spect.nii.gz",
        "labels": "labels.nii.gz",
        "ct": "ct.nii.gz",
        "method": "hanscheid2017",
        "radionuclide": "Lu177",
        "total_activity": 250.0,
        "rois": ["Liver"],
        "teff": {"Liver": 45.2},
        "output": str(output_json),
    }
    with open(config_file, "w") as f:
        json.dump(cfg, f)

    # Pass config file as first argument
    results = compute_roi_dose_file(config_file)
    assert "Liver" in results
    assert output_json.exists()


def test_cli_dose_with_config_positional(dummy_study):
    study_dir = dummy_study["dir"]
    config_file = study_dir / "config.json"
    output_json = study_dir / "cli_out.json"
    cfg = {
        "image": "spect.nii.gz",
        "labels": "labels.nii.gz",
        "ct": "ct.nii.gz",
        "method": "hanscheid2017",
        "radionuclide": "Lu177",
        "total_activity": 250.0,
        "rois": ["Liver"],
        "teff": {"Liver": 45.2},
        "output": str(output_json),
    }
    with open(config_file, "w") as f:
        json.dump(cfg, f)

    runner = CliRunner()
    result = runner.invoke(main, ["dose", str(config_file)])
    assert result.exit_code == 0, result.output
    assert output_json.exists()

    with open(output_json, "r") as f:
        data = json.load(f)
    assert "Liver" in data
    assert data["Liver"]["roi_teff_h"] == 45.2


def test_cli_dose_with_config_flag(dummy_study):
    study_dir = dummy_study["dir"]
    config_file = study_dir / "config_flag.json"
    output_json = study_dir / "flag_out.json"
    cfg = {
        "image": "spect.nii.gz",
        "labels": "labels.nii.gz",
        "ct": "ct.nii.gz",
        "method": "hanscheid2017",
        "radionuclide": "Lu177",
        "total_activity": 250.0,
        "rois": ["Liver"],
        "teff": {"Liver": 45.2},
        "output": str(output_json),
    }
    with open(config_file, "w") as f:
        json.dump(cfg, f)

    runner = CliRunner()
    result = runner.invoke(main, ["dose", "-c", str(config_file)])
    assert result.exit_code == 0, result.output
    assert output_json.exists()


def test_cli_dose_override_config(dummy_study):
    study_dir = dummy_study["dir"]
    config_file = study_dir / "config_override.json"
    default_out = study_dir / "default_out.json"
    override_out = study_dir / "override_out.json"

    cfg = {
        "image": "spect.nii.gz",
        "labels": "labels.nii.gz",
        "ct": "ct.nii.gz",
        "method": "hanscheid2017",
        "radionuclide": "Lu177",
        "total_activity": 250.0,
        "rois": ["Liver"],
        "teff": {"Liver": 45.2},
        "output": str(default_out),
    }
    with open(config_file, "w") as f:
        json.dump(cfg, f)

    runner = CliRunner()
    # Override output and teff from CLI
    result = runner.invoke(
        main,
        [
            "dose",
            str(config_file),
            "-o",
            str(override_out),
            "--teff",
            "Liver:50.0",
        ],
    )
    assert result.exit_code == 0, result.output
    assert not default_out.exists()
    assert override_out.exists()

    with open(override_out, "r") as f:
        data = json.load(f)
    assert data["Liver"]["roi_teff_h"] == 50.0


def test_cli_dose_missing_params(tmp_path):
    runner = CliRunner()

    # 1. Neither image nor config
    res1 = runner.invoke(main, ["dose"])
    assert res1.exit_code != 0
    assert "Missing argument 'IMAGE_PATH' or option '-c' / '--config'" in res1.output

    # 2. Incomplete config file (missing method)
    incomplete_cfg = tmp_path / "bad.json"
    with open(incomplete_cfg, "w") as f:
        json.dump({"image": "nonexistent.nii.gz", "radionuclide": "Lu177"}, f)

    res2 = runner.invoke(main, ["dose", str(incomplete_cfg)])
    assert res2.exit_code != 0
    assert "Missing option '-m' / '--method'" in res2.output


def test_load_dose_config_grouped_rois(dummy_study):
    study_dir = dummy_study["dir"]

    # 1. Test dictionary with float values
    cfg_file1 = study_dir / "cfg_grouped_dict.json"
    with open(cfg_file1, "w") as f:
        json.dump({
            "image": "spect.nii.gz",
            "labels": "labels.nii.gz",
            "method": "madsen2018",
            "isotope": "Lu177",
            "rois": {
                "Liver": 45.2,
                "Spleen": {"teff": 38.0}
            }
        }, f)

    loaded1 = load_dose_config(cfg_file1)
    assert set(loaded1["rois"]) == {"Liver", "Spleen"}
    assert loaded1["method_params"]["Liver"]["teff"] == 45.2
    assert loaded1["method_params"]["Spleen"]["teff"] == 38.0

    # 2. Test list of dicts
    cfg_file2 = study_dir / "cfg_grouped_list.json"
    with open(cfg_file2, "w") as f:
        json.dump({
            "image": "spect.nii.gz",
            "labels": "labels.nii.gz",
            "method": "madsen2018",
            "isotope": "Lu177",
            "rois": [
                {"name": "Liver", "teff": 45.2},
                {"name": "Spleen", "teff": 38.0}
            ]
        }, f)

    loaded2 = load_dose_config(cfg_file2)
    assert set(loaded2["rois"]) == {"Liver", "Spleen"}
    assert loaded2["method_params"]["Liver"]["teff"] == 45.2
    assert loaded2["method_params"]["Spleen"]["teff"] == 38.0

