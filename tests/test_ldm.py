import pytest
import numpy as np
import SimpleITK as sitk
import json
from click.testing import CliRunner
from radot.cli import main
from radot.db import get_non_penetrating_energy_constant_mev
from radot.dose import compute_roi_dose

def test_get_non_penetrating_energy_constant_mev():
    # Lu-177 non-penetrating energy
    lu_const = get_non_penetrating_energy_constant_mev("Lu177")
    assert lu_const > 0.1  # Typically around 0.14 MeV

    # Y-90 non-penetrating energy
    y_const = get_non_penetrating_energy_constant_mev("Y90")
    assert y_const > 0.8  # Typically around 0.93 MeV


def test_ldm_core_calculation_without_ct():
    # Setup test SPECT image (2x2x2) with voxel values of 1e6 Bq
    image = sitk.Image([2, 2, 2], sitk.sitkFloat32)
    # Voxel spacing: 2mm isotropic
    image.SetSpacing([2.0, 2.0, 2.0])
    image_arr = np.ones((2, 2, 2), dtype=np.float32) * 1e6
    
    # Write array to SITK image
    sitk_image = sitk.GetImageFromArray(image_arr)
    sitk_image.SetSpacing([2.0, 2.0, 2.0])
    sitk_image.SetMetaData("time_from_injection", "24.0")
    sitk_image.SetMetaData("unit", "Bq")

    # Run LDM without CT (assumes water density = 1.0 g/cm3)
    dose_img = compute_roi_dose(
        image=sitk_image,
        labels=None,
        method="ldm",
        isotope="Lu177",
        ct_image=None,
        verbose=True
    )

    assert isinstance(dose_img, sitk.Image)
    assert dose_img.GetSize() == (2, 2, 2)
    assert dose_img.GetSpacing() == (2.0, 2.0, 2.0)

    # Let's inspect the calculated voxel values
    dose_arr = sitk.GetArrayFromImage(dose_img)
    # Dose value should be > 0.0 Gy since voxel has activity
    assert np.all(dose_arr > 0.0)


def test_ldm_core_calculation_with_ct():
    # Setup activity image
    image = sitk.Image([2, 2, 2], sitk.sitkFloat32)
    image.SetSpacing([2.0, 2.0, 2.0])
    image_arr = np.ones((2, 2, 2), dtype=np.float32) * 1e6
    sitk_image = sitk.GetImageFromArray(image_arr)
    sitk_image.SetSpacing([2.0, 2.0, 2.0])
    sitk_image.SetMetaData("time_from_injection", "24.0")
    sitk_image.SetMetaData("unit", "Bq")

    # Setup CT image with density corresponding to bone (HU = 1000)
    ct_arr = np.ones((2, 2, 2), dtype=np.int16) * 1000
    ct_image = sitk.GetImageFromArray(ct_arr)
    ct_image.SetSpacing([2.0, 2.0, 2.0])

    # Run LDM with CT
    dose_img = compute_roi_dose(
        image=sitk_image,
        labels=None,
        method="ldm",
        isotope="Lu177",
        ct_image=ct_image,
        verbose=True
    )

    dose_arr = sitk.GetArrayFromImage(dose_img)
    assert np.all(dose_arr > 0.0)


def test_ldm_cli_execution(tmp_path):
    runner = CliRunner()

    image_path = tmp_path / "spect.nii.gz"
    ct_path = tmp_path / "ct.nii.gz"
    output_path = tmp_path / "dose.nii.gz"

    # Save SPECT image
    img = sitk.Image([3, 3, 3], sitk.sitkFloat32)
    img.SetSpacing([2.0, 2.0, 2.0])
    sitk.WriteImage(img, str(image_path))

    with open(tmp_path / "spect.json", "w") as f:
        json.dump({"unit": "Bq", "time_from_injection": "24.0"}, f)

    # Save CT image
    ct = sitk.Image([3, 3, 3], sitk.sitkFloat32)
    ct.SetSpacing([2.0, 2.0, 2.0])
    sitk.WriteImage(ct, str(ct_path))

    # Run LDM CLI command
    result = runner.invoke(main, [
        "dose",
        str(image_path),
        "-m", "ldm",
        "--isotope", "Lu177",
        "--ct", str(ct_path),
        "-o", str(output_path)
    ])

    assert result.exit_code == 0
    assert output_path.exists()

    json_path = output_path.with_suffix("").with_suffix(".json")
    assert json_path.exists()
    with open(json_path, "r") as f:
        meta_data = json.load(f)
    assert meta_data["unit"] == "gy"
    assert meta_data["modality"] == "dose"


def test_ldm_core_calculation_with_labels_no_output():
    # Setup test SPECT image (2x2x2) with voxel values of 1e6 Bq
    image = sitk.Image([2, 2, 2], sitk.sitkFloat32)
    image.SetSpacing([2.0, 2.0, 2.0])
    image_arr = np.ones((2, 2, 2), dtype=np.float32) * 1e6
    sitk_image = sitk.GetImageFromArray(image_arr)
    sitk_image.SetSpacing([2.0, 2.0, 2.0])
    sitk_image.SetMetaData("time_from_injection", "24.0")
    sitk_image.SetMetaData("unit", "Bq")

    # Setup labels image
    labels = sitk.Image([2, 2, 2], sitk.sitkUInt16)
    labels.SetSpacing([2.0, 2.0, 2.0])
    labels[1, 1, 1] = 1

    # Run LDM with labels
    results = compute_roi_dose(
        image=sitk_image,
        labels=labels,
        method="ldm",
        isotope="Lu177",
        labels_map={1: "Liver"},
        verbose=True
    )

    assert isinstance(results, tuple)
    results_dict, dose_img = results
    assert "Liver" in results_dict
    assert results_dict["Liver"]["mean_dose_gy"] > 0.0
    assert results_dict["Liver"]["roi_volume_ml"] > 0.0
    assert isinstance(dose_img, sitk.Image)


def test_ldm_cli_with_labels_no_output(tmp_path):
    runner = CliRunner()

    image_path = tmp_path / "spect.nii.gz"
    labels_path = tmp_path / "labels.nii.gz"

    # Save SPECT image
    img = sitk.Image([3, 3, 3], sitk.sitkFloat32)
    img.SetSpacing([2.0, 2.0, 2.0])
    sitk.WriteImage(img, str(image_path))

    with open(tmp_path / "spect.json", "w") as f:
        json.dump({"unit": "Bq", "time_from_injection": "24.0"}, f)

    # Save labels image
    labels = sitk.Image([3, 3, 3], sitk.sitkUInt16)
    labels.SetSpacing([2.0, 2.0, 2.0])
    labels[1, 1, 1] = 1
    sitk.WriteImage(labels, str(labels_path))

    with open(tmp_path / "labels.json", "w") as f:
        json.dump({"1": "Liver"}, f)

    # Run LDM CLI command without output path option -o
    result = runner.invoke(main, [
        "dose",
        str(image_path),
        "-l", str(labels_path),
        "-m", "ldm",
        "--isotope", "Lu177"
    ])

    assert result.exit_code == 0
    assert "Liver" in result.output


def test_ldm_hu_threshold():
    # Activity image (2x2x2) with activity = 1e6 Bq
    image = sitk.Image([2, 2, 2], sitk.sitkFloat32)
    image.SetSpacing([2.0, 2.0, 2.0])
    image_arr = np.ones((2, 2, 2), dtype=np.float32) * 1e6
    sitk_image = sitk.GetImageFromArray(image_arr)
    sitk_image.SetSpacing([2.0, 2.0, 2.0])
    sitk_image.SetMetaData("time_from_injection", "24.0")
    sitk_image.SetMetaData("unit", "Bq")

    # CT image where one voxel has HU = -900, another has HU = 0
    ct_arr = np.zeros((2, 2, 2), dtype=np.int16)
    ct_arr[0, 0, 0] = -900
    ct_image = sitk.GetImageFromArray(ct_arr)
    ct_image.SetSpacing([2.0, 2.0, 2.0])

    # Run LDM with hu_threshold = -800
    dose_img = compute_roi_dose(
        image=sitk_image,
        labels=None,
        method="ldm",
        isotope="Lu177",
        ct_image=ct_image,
        hu_threshold=-800.0,
        verbose=True
    )

    dose_arr = sitk.GetArrayFromImage(dose_img)
    # The voxel at (0, 0, 0) should be zeroed out
    assert dose_arr[0, 0, 0] == 0.0
    # The other voxels (HU = 0 >= -800) should have dose > 0.0
    assert dose_arr[1, 1, 1] > 0.0


def test_ldm_dose_rate(tmp_path):
    runner = CliRunner()

    image_path = tmp_path / "spect.nii.gz"
    labels_path = tmp_path / "labels.nii.gz"
    output_path = tmp_path / "dose_rate.nii.gz"

    # Save SPECT image
    img = sitk.Image([3, 3, 3], sitk.sitkFloat32)
    img.SetSpacing([2.0, 2.0, 2.0])
    sitk.WriteImage(img, str(image_path))

    with open(tmp_path / "spect.json", "w") as f:
        json.dump({"unit": "Bq", "time_from_injection": "24.0"}, f)

    # Save labels image
    labels = sitk.Image([3, 3, 3], sitk.sitkUInt16)
    labels.SetSpacing([2.0, 2.0, 2.0])
    labels[1, 1, 1] = 1
    sitk.WriteImage(labels, str(labels_path))

    with open(tmp_path / "labels.json", "w") as f:
        json.dump({"1": "Liver"}, f)

    # Run LDM CLI command with --dose-rate
    result = runner.invoke(main, [
        "dose",
        str(image_path),
        "-l", str(labels_path),
        "-m", "ldm",
        "--isotope", "Lu177",
        "-o", str(output_path),
        "--dose-rate"
    ])

    assert result.exit_code == 0
    assert "Mean Dose Rate (Gy/h)" in result.output
    assert output_path.exists()

    json_path = output_path.with_suffix("").with_suffix(".json")
    assert json_path.exists()
    with open(json_path, "r") as f:
        meta_data = json.load(f)
    assert meta_data["unit"] == "gy/h"
    assert meta_data["modality"] == "dose"



