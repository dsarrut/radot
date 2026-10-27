import json
import numpy as np
import pytest
import SimpleITK as sitk
from click.testing import CliRunner
from radot.cli import main
from radot.dose import (
    compute_roi_dose,
    roi_dose_uncertainty_bounds,
    roi_dose_uncertainty_estimate,
    roi_mean_chord_mm,
)


def test_compute_roi_dose_core():
    # Create dummy images with appropriate metadata tags
    image = sitk.Image([10, 10, 10], sitk.sitkFloat32)
    image.SetMetaData("time_from_injection", "24.0")
    image.SetMetaData("unit", "Bq")

    labels = sitk.Image([10, 10, 10], sitk.sitkUInt16)
    labels[2, 2, 2] = 1  # Mark one voxel with label 1

    # We call with rois=None so it doesn't filter out our "unknown_1" label
    res = compute_roi_dose(
        image=image,
        labels=labels,
        method="hanscheid2017",
        isotope="Lu177",
        rois=None,
        total_activity_mbq=250.0,
        method_params={"unknown_1": {"teff": 45.2}},
    )

    assert "unknown_1" in res
    # With a dummy empty image, concentration is 0, so mean dose is 0
    assert res["unknown_1"]["mean_dose_gy"] == 0.0
    assert res["unknown_1"]["roi_teff_h"] == 45.2


def test_physical_half_life_fallback():
    # Create dummy images
    image = sitk.Image([10, 10, 10], sitk.sitkFloat32)
    image.SetMetaData("time_from_injection", "24.0")
    image.SetMetaData("unit", "Bq")

    labels = sitk.Image([10, 10, 10], sitk.sitkUInt16)
    labels[2, 2, 2] = 1

    res = compute_roi_dose(
        image=image,
        labels=labels,
        method="hanscheid2017",
        isotope="Lu177",
        rois=None,
        total_activity_mbq=None,
        method_params=None,  # No Teff provided, should fall back to physical half life
    )

    assert "unknown_1" in res
    # Physical half life of Lu177 is 159.528 hours
    assert res["unknown_1"]["roi_teff_h"] == pytest.approx(159.528)


def test_madsen2018():
    image = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    image.SetSpacing((10.0, 10.0, 10.0))
    image.SetMetaData("time_from_injection", "24.0")
    image.SetMetaData("unit", "Bq")
    image[2, 2, 2] = 100.0

    labels = sitk.Image([5, 5, 5], sitk.sitkUInt16)
    labels.SetSpacing((10.0, 10.0, 10.0))
    labels[2, 2, 2] = 1

    labels_map = {1: "Left kidney"}

    # Check that it raises an error if teff is missing
    with pytest.raises(ValueError, match="Effective half-life .* must be provided"):
        compute_roi_dose(
            image=image,
            labels=labels,
            method="madsen2018",
            isotope="Lu177",
            labels_map=labels_map,
            rois=None,
        )

    res = compute_roi_dose(
        image=image,
        labels=labels,
        method="madsen2018",
        isotope="Lu177",
        labels_map=labels_map,
        rois=None,
        method_params={"phantom": "male", "Left kidney": {"teff": 45.2}},
    )
    assert "Left kidney" in res
    assert res["Left kidney"]["mean_dose_gy"] > 0.0
    assert res["Left kidney"]["roi_teff_h"] == 45.2


def test_hanscheid2018_custom_teff():
    image = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    image.SetSpacing((10.0, 10.0, 10.0))
    image.SetMetaData("time_from_injection", "24.0")
    image.SetMetaData("unit", "Bq")
    image[2, 2, 2] = 100.0

    labels = sitk.Image([5, 5, 5], sitk.sitkUInt16)
    labels.SetSpacing((10.0, 10.0, 10.0))
    labels[2, 2, 2] = 1

    labels_map = {1: "Left kidney"}

    # Default Teff = 2 * 24 = 48 h
    res_default = compute_roi_dose(
        image=image,
        labels=labels,
        method="hanscheid2018",
        isotope="Lu177",
        labels_map=labels_map,
        method_params={"phantom": "male"},
    )

    # Custom Teff = 36 h
    res_custom = compute_roi_dose(
        image=image,
        labels=labels,
        method="hanscheid2018",
        isotope="Lu177",
        labels_map=labels_map,
        method_params={"phantom": "male", "Left kidney": {"teff": 36.0}},
    )

    assert res_default["Left kidney"]["roi_teff_h"] == 48.0
    assert res_custom["Left kidney"]["roi_teff_h"] == 36.0
    # Ratio of doses should be exactly 36 / 48 = 0.75
    ratio = res_custom["Left kidney"]["mean_dose_gy"] / res_default["Left kidney"]["mean_dose_gy"]
    assert ratio == pytest.approx(36.0 / 48.0, rel=1e-5)


def test_invalid_isotope_error():

    image = sitk.Image([10, 10, 10], sitk.sitkFloat32)
    image.SetMetaData("time_from_injection", "24.0")
    image.SetMetaData("unit", "Bq")

    labels = sitk.Image([10, 10, 10], sitk.sitkUInt16)

    with pytest.raises(ValueError, match="Hanscheid 2017 method .* specific to Lu177"):
        compute_roi_dose(
            image=image,
            labels=labels,
            method="hanscheid2017",
            isotope="Y90",
        )


def test_dose_cli(tmp_path):
    image_path = tmp_path / "spect.nii.gz"
    labels_path = tmp_path / "labels.nii.gz"
    ct_path = tmp_path / "ct.nii.gz"
    output_path = tmp_path / "dose_out.json"

    # Create dummy images with metadata tags
    img = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    img.SetMetaData("time_from_injection", "24.0")
    img.SetMetaData("unit", "Bq")
    sitk.WriteImage(img, str(image_path))

    # Write a companion spect.json sidecar to define metadata
    with open(tmp_path / "spect.json", "w") as f:
        json.dump({"unit": "Bq", "time_from_injection": "24.0"}, f)

    lbls = sitk.Image([5, 5, 5], sitk.sitkUInt16)
    lbls[2, 2, 2] = 1
    sitk.WriteImage(lbls, str(labels_path))

    # Write a companion labels.json sidecar to define the Liver label
    labels_json_path = tmp_path / "labels.json"
    with open(labels_json_path, "w") as f:
        json.dump({"1": "Liver"}, f)

    ct_img = sitk.Image([5, 5, 5], sitk.sitkInt16)
    sitk.WriteImage(ct_img, str(ct_path))

    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "dose",
            str(image_path),
            "-l",
            str(labels_path),
            "--ct",
            str(ct_path),
            "-m",
            "hanscheid2017",
            "--isotope",
            "Lu177",
            "--total-activity",
            "250.0",
            "--rois",
            "Liver",
            "--teff",
            "Liver:45.2",
            "-o",
            str(output_path),
        ],
    )

    assert result.exit_code == 0, f"CLI output: {result.output}"
    assert output_path.exists()

    with open(output_path, "r") as f:
        data = json.load(f)

    assert "Liver" in data
    assert data["Liver"]["roi_teff_h"] == 45.2


def test_grid_mismatch_error():
    image = sitk.Image([10, 10, 10], sitk.sitkFloat32)
    image.SetMetaData("time_from_injection", "24.0")
    image.SetMetaData("unit", "Bq")

    labels = sitk.Image([5, 5, 5], sitk.sitkUInt16)  # Mismatched size

    with pytest.raises(
        ValueError, match="Activity image grid does not match labels grid"
    ):
        compute_roi_dose(
            image=image,
            labels=labels,
            method="hanscheid2017",
            isotope="Lu177",
        )


def test_missing_roi_in_sidecar_error(tmp_path):
    labels_path = tmp_path / "labels.nii.gz"
    labels = sitk.Image([10, 10, 10], sitk.sitkUInt16)
    sitk.WriteImage(labels, str(labels_path))

    # Write label mapping mapping only "Liver"
    with open(tmp_path / "labels.json", "w") as f:
        json.dump({"1": "Liver"}, f)

    image = sitk.Image([10, 10, 10], sitk.sitkFloat32)
    image.SetMetaData("time_from_injection", "24.0")
    image.SetMetaData("unit", "Bq")

    # Requesting "Kidney" which is not in sidecar
    with pytest.raises(ValueError, match="Requested ROI 'Kidney' not found"):
        compute_roi_dose(
            image=image,
            labels=labels,
            method="hanscheid2017",
            isotope="Lu177",
            rois=["Kidney"],
            labels_path=labels_path,
        )


def test_dose_cli_resampling(tmp_path):
    image_path = tmp_path / "spect.nii.gz"
    labels_path = tmp_path / "labels.nii.gz"
    output_path = tmp_path / "dose_out.json"

    # Mismatched image grid (size 10x10x10) vs labels grid (size 5x5x5)
    img = sitk.Image([10, 10, 10], sitk.sitkFloat32)
    img.SetMetaData("time_from_injection", "24.0")
    img.SetMetaData("unit", "Bq")
    sitk.WriteImage(img, str(image_path))

    with open(tmp_path / "spect.json", "w") as f:
        json.dump({"unit": "Bq", "time_from_injection": "24.0"}, f)

    lbls = sitk.Image([5, 5, 5], sitk.sitkUInt16)
    lbls[2, 2, 2] = 1
    sitk.WriteImage(lbls, str(labels_path))

    with open(tmp_path / "labels.json", "w") as f:
        json.dump({"1": "Liver"}, f)

    runner = CliRunner()

    # 1. Run without resampling: should fail due to mismatched grids
    result_fail = runner.invoke(
        main,
        [
            "dose",
            str(image_path),
            "-l",
            str(labels_path),
            "-m",
            "hanscheid2017",
            "--isotope",
            "Lu177",
            "-o",
            str(output_path),
        ],
    )
    assert result_fail.exit_code != 0
    assert "Activity image grid does not match labels grid" in result_fail.output

    # 2. Run with resampling: should succeed
    result_success = runner.invoke(
        main,
        [
            "dose",
            str(image_path),
            "-l",
            str(labels_path),
            "-m",
            "hanscheid2017",
            "--isotope",
            "Lu177",
            "--resample",
            "like_activity",
            "-o",
            str(output_path),
        ],
    )
    assert result_success.exit_code == 0
    assert output_path.exists()


def test_dose_cli_verbose(tmp_path):
    image_path = tmp_path / "spect.nii.gz"
    labels_path = tmp_path / "labels.nii.gz"

    img = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    img.SetMetaData("time_from_injection", "24.0")
    img.SetMetaData("unit", "Bq")
    sitk.WriteImage(img, str(image_path))

    # Write companion spect.json sidecar to define metadata
    with open(tmp_path / "spect.json", "w") as f:
        json.dump({"unit": "Bq", "time_from_injection": "24.0"}, f)

    lbls = sitk.Image([5, 5, 5], sitk.sitkUInt16)
    lbls[2, 2, 2] = 1
    sitk.WriteImage(lbls, str(labels_path))

    with open(tmp_path / "labels.json", "w") as f:
        json.dump({"1": "Liver"}, f)

    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "dose",
            str(image_path),
            "-l",
            str(labels_path),
            "-m",
            "hanscheid2017",
            "--isotope",
            "Lu177",
            "--teff",
            "Liver:45.2",
            "--verbose",
        ],
    )

    assert result.exit_code == 0
    assert (
        "Dosimetry calculation log for ROI 'Liver' using Hanscheid 2017:"
        in result.output
    )
    assert "Effective half-life (Teff): 45.2" in result.output
    assert "Formula: 0.125 * " in result.output


def test_dose_cli_teff_variants(tmp_path):
    image_path = tmp_path / "spect.nii.gz"
    labels_path = tmp_path / "labels.nii.gz"

    img = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    img.SetMetaData("time_from_injection", "24.0")
    img.SetMetaData("unit", "Bq")
    sitk.WriteImage(img, str(image_path))

    with open(tmp_path / "spect.json", "w") as f:
        json.dump({"unit": "Bq", "time_from_injection": "24.0"}, f)

    lbls = sitk.Image([5, 5, 5], sitk.sitkUInt16)
    lbls[2, 2, 2] = 1
    sitk.WriteImage(lbls, str(labels_path))

    with open(tmp_path / "labels.json", "w") as f:
        json.dump({"1": "Liver"}, f)

    runner = CliRunner()

    # 1. Test that global/single numeric value for teff fails (must be ROI:value format)
    result_global = runner.invoke(
        main,
        [
            "dose",
            str(image_path),
            "-l",
            str(labels_path),
            "-m",
            "hanscheid2017",
            "--isotope",
            "Lu177",
            "--teff",
            "45",
        ],
    )
    assert result_global.exit_code != 0
    assert "must be in ROI:value format" in result_global.output

    # 2. Test invalid teff formatting error
    result_invalid = runner.invoke(
        main,
        [
            "dose",
            str(image_path),
            "-l",
            str(labels_path),
            "-m",
            "hanscheid2017",
            "--isotope",
            "Lu177",
            "--teff",
            "Liver:abc",
        ],
    )
    assert result_invalid.exit_code != 0
    assert "Invalid teff value" in result_invalid.output


def test_dose_cli_multiple_and_case_insensitive(tmp_path):
    image_path = tmp_path / "spect.nii.gz"
    labels_path = tmp_path / "labels.nii.gz"

    img = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    img.SetMetaData("time_from_injection", "24.0")
    img.SetMetaData("unit", "Bq")
    sitk.WriteImage(img, str(image_path))

    with open(tmp_path / "spect.json", "w") as f:
        json.dump({"unit": "Bq", "time_from_injection": "24.0"}, f)

    lbls = sitk.Image([5, 5, 5], sitk.sitkUInt16)
    lbls[1, 1, 1] = 1
    lbls[2, 2, 2] = 2
    sitk.WriteImage(lbls, str(labels_path))

    with open(tmp_path / "labels.json", "w") as f:
        json.dump({"1": "Liver", "2": "Kidney_L"}, f)

    runner = CliRunner()

    # 1. Success with case-insensitive multiple rois and teffs
    result = runner.invoke(
        main,
        [
            "dose",
            str(image_path),
            "-l",
            str(labels_path),
            "-m",
            "hanscheid2017",
            "--isotope",
            "Lu177",
            "--rois",
            "liver",
            "--rois",
            "kidney_l",
            "--teff",
            "LIVER:45",
            "--teff",
            "KIDNEY_L:32",
            "--verbose",
        ],
    )
    assert result.exit_code == 0
    assert (
        "Dosimetry calculation log for ROI 'Liver' using Hanscheid 2017:"
        in result.output
    )
    assert "Effective half-life (Teff): 45.000 h" in result.output
    assert (
        "Dosimetry calculation log for ROI 'Kidney_L' using Hanscheid 2017:"
        in result.output
    )
    assert "Effective half-life (Teff): 32.000 h" in result.output

    # 2. Failure with typo suggestion for kiney
    result_typo = runner.invoke(
        main,
        [
            "dose",
            str(image_path),
            "-l",
            str(labels_path),
            "-m",
            "hanscheid2017",
            "--isotope",
            "Lu177",
            "--teff",
            "kiney:32",
        ],
    )
    assert result_typo.exit_code != 0
    assert "Teff ROI 'kiney' not found. Did you mean" in result_typo.output
    assert "'Kidney_L'" in result_typo.output


def test_compute_roi_dose_string_inputs(tmp_path):
    image_path = tmp_path / "spect.nii.gz"
    labels_path = tmp_path / "labels.nii.gz"
    ct_path = tmp_path / "ct.nii.gz"

    img = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    img.SetMetaData("time_from_injection", "24.0")
    img.SetMetaData("unit", "Bq")
    sitk.WriteImage(img, str(image_path))

    with open(tmp_path / "spect.json", "w") as f:
        json.dump({"unit": "Bq", "time_from_injection": "24.0"}, f)

    lbls = sitk.Image([5, 5, 5], sitk.sitkUInt16)
    lbls[2, 2, 2] = 1
    sitk.WriteImage(lbls, str(labels_path))

    with open(tmp_path / "labels.json", "w") as f:
        json.dump({"1": "Liver"}, f)

    ct_img = sitk.Image([5, 5, 5], sitk.sitkInt16)
    sitk.WriteImage(ct_img, str(ct_path))

    # Verify we can pass file paths directly as strings to compute_roi_dose
    # and pass ct_path/teff_params directly.
    res = compute_roi_dose(
        image=str(image_path),
        labels=str(labels_path),
        method="hanscheid2017",
        isotope="Lu177",
        ct_path=str(ct_path),
        teff_params={"Liver": 45.2},
    )

    assert "Liver" in res
    assert res["Liver"]["roi_teff_h"] == 45.2


def test_compute_roi_dose_hanscheid2018():
    # Create dummy images with appropriate metadata tags
    # Voxel size 10x10x10 mm -> voxel volume = 1 mL
    image = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    image.SetSpacing((10.0, 10.0, 10.0))
    image.SetMetaData("time_from_injection", "24.0")
    image.SetMetaData("unit", "Bq")
    # Put 100 Bq in voxel (2,2,2)
    image[2, 2, 2] = 100.0

    labels = sitk.Image([5, 5, 5], sitk.sitkUInt16)
    labels.SetSpacing((10.0, 10.0, 10.0))
    labels[2, 2, 2] = 1  # Label 1 mapping to Left kidney

    labels_map = {1: "Left kidney"}

    # 1. Test male phantom
    res_male = compute_roi_dose(
        image=image,
        labels=labels,
        method="hanscheid2018",
        isotope="Lu177",
        labels_map=labels_map,
        method_params={"phantom": "male"},
        rois=None,
    )
    assert "Left kidney" in res_male
    # Expected dose: ~0.0005974 Gy
    assert res_male["Left kidney"]["mean_dose_gy"] == pytest.approx(
        0.000597389, rel=1e-3
    )
    assert res_male["Left kidney"]["roi_teff_h"] == 48.0

    # 2. Test female phantom
    res_female = compute_roi_dose(
        image=image,
        labels=labels,
        method="hanscheid2018",
        isotope="Lu177",
        labels_map=labels_map,
        method_params={"phantom": "female"},
        rois=None,
    )
    assert "Left kidney" in res_female
    # Expected dose: ~0.0005975 Gy
    assert res_female["Left kidney"]["mean_dose_gy"] == pytest.approx(
        0.000597525, rel=1e-3
    )


def test_dose_missing_organ_halt():
    import numpy as np
    from radot.dose import compute_roi_dose
    import click

    image = sitk.Image([3, 3, 3], sitk.sitkFloat32)
    image.SetMetaData("time_from_injection", "24.0")
    image.SetMetaData("unit", "Bq")

    labels = sitk.Image([3, 3, 3], sitk.sitkUInt16)
    labels[1, 1, 1] = 1

    # prostate is in male, but NOT in female
    labels_map = {1: "Prostate"}

    # It should raise ValueError with instructions rather than skipping
    with pytest.raises(
        ValueError, match="Dosimetry Halt:.*Prostate.*not found in S-value matrix"
    ):
        compute_roi_dose(
            image=image,
            labels=labels,
            method="hanscheid2018",
            isotope="Lu177",
            labels_map=labels_map,
            method_params={"phantom": "female"},
            rois=["Prostate"],
        )


def test_dose_missing_organ_loop_warning():
    import numpy as np
    from radot.dose import compute_roi_dose

    image = sitk.Image([3, 3, 3], sitk.sitkFloat32)
    image.SetMetaData("time_from_injection", "24.0")
    image.SetMetaData("unit", "Bq")

    labels = sitk.Image([3, 3, 3], sitk.sitkUInt16)
    labels[1, 1, 1] = 1

    # prostate is in male, but NOT in female
    labels_map = {1: "Prostate"}

    res = compute_roi_dose(
        image=image,
        labels=labels,
        method="hanscheid2018",
        isotope="Lu177",
        labels_map=labels_map,
        method_params={"phantom": "female"},
        rois=None,
    )

    # It should skip "Prostate" since it does not exist in the female S-value database,
    # rather than crashing. The results dict for Prostate will be empty/omitted.
    assert "Prostate" not in res


def test_dose_list_command():
    from click.testing import CliRunner
    from radot.cli import main

    runner = CliRunner()

    result = runner.invoke(main, ["dose-list"])
    assert result.exit_code == 0
    assert "hanscheid2017" in result.output
    assert "Hanscheid 2017 method" in result.output
    assert "Equation: Dose = 0.125 * " in result.output
    assert "hanscheid2018" in result.output
    assert "uses OpenDose S-value database" in result.output
    assert "Equation: Dose = TIA * S * " in result.output
    assert "ldm" in result.output
    assert "Local Deposition Model" in result.output
    assert "Equation: Dose = (TIA * Delta_local * " in result.output
    assert "Equation: Dose Rate = (Activity * Delta_local * " in result.output


def test_mutually_exclusive_verbose_quiet(tmp_path):
    from click.testing import CliRunner
    from radot.cli import main

    runner = CliRunner()

    # Test global options
    result = runner.invoke(main, ["--verbose", "--quiet", "dose-list"])
    assert result.exit_code != 0
    assert "mutually exclusive" in result.output

    # Create dummy file so exists=True check passes
    dummy_file = tmp_path / "dummy.nii.gz"
    dummy_file.touch()

    # Test command-level options
    result_cmd = runner.invoke(
        main,
        [
            "dose",
            str(dummy_file),
            "-m",
            "madsen2018",
            "--isotope",
            "Lu177",
            "--verbose",
            "--quiet",
        ],
    )
    assert result_cmd.exit_code != 0
    assert "mutually exclusive" in result_cmd.output


def test_roi_dose_uncertainty_bounds():
    # 1. Single voxel: bounds must equal the voxel uncertainty
    dose = np.array([[[10.0]]], dtype=np.float32)
    unc = np.array([[[0.05]]], dtype=np.float32)
    mask = np.array([[[1]]], dtype=np.uint8)

    unc_min, unc_max = roi_dose_uncertainty_bounds(dose, unc, mask)
    assert unc_min == pytest.approx(0.05)
    assert unc_max == pytest.approx(0.05)

    # 2. 4 identical voxels:
    # lower bound = 0.05 / sqrt(4) = 0.025
    # upper bound = 0.05
    dose_4 = np.ones((2, 2, 1), dtype=np.float32) * 5.0
    unc_4 = np.ones((2, 2, 1), dtype=np.float32) * 0.05
    mask_4 = np.ones((2, 2, 1), dtype=np.uint8)

    unc_min_4, unc_max_4 = roi_dose_uncertainty_bounds(dose_4, unc_4, mask_4)
    assert unc_min_4 == pytest.approx(0.025)
    assert unc_max_4 == pytest.approx(0.05)

    # 3. SimpleITK image input
    sitk_dose = sitk.GetImageFromArray(dose_4)
    sitk_unc = sitk.GetImageFromArray(unc_4)
    sitk_mask = sitk.GetImageFromArray(mask_4)

    unc_min_sitk, unc_max_sitk = roi_dose_uncertainty_bounds(
        sitk_dose, sitk_unc, sitk_mask
    )
    assert unc_min_sitk == pytest.approx(0.025)
    assert unc_max_sitk == pytest.approx(0.05)

    # 4. Empty mask raises error
    empty_mask = np.zeros_like(mask_4)
    with pytest.raises(ValueError, match="ROI mask contains no voxels"):
        roi_dose_uncertainty_bounds(dose_4, unc_4, empty_mask)


def test_roi_dose_uncertainty_estimate():
    # 4 voxels with large spacing 10x10x10 mm, isotope Lu177 (mean beta range = 0.23 mm)
    # nu_cluster = (1 + 0.23/10)^3
    dose = np.ones((2, 2, 1), dtype=np.float32)
    unc = np.ones((2, 2, 1), dtype=np.float32) * 0.05
    mask = np.ones((2, 2, 1), dtype=np.uint8)

    unc_est = roi_dose_uncertainty_estimate(
        dose, unc, mask, spacing_mm=(10.0, 10.0, 10.0), isotope="Lu177"
    )
    # unc_min = 0.025, unc_max = 0.05
    assert unc_est == pytest.approx(0.025 * np.sqrt(1.023**3), rel=1e-3)
    assert unc_est >= 0.025
    assert unc_est <= 0.05

    # Photons: nu_cluster = 1 + chord / spacing, capped by unc_max
    unc_ph = roi_dose_uncertainty_estimate(
        dose, unc, mask, spacing_mm=(10.0, 10.0, 10.0), photons=True, correlation_length_mm=10.0
    )
    assert unc_ph == pytest.approx(0.025 * np.sqrt(2.0), rel=1e-6)
    unc_ph = roi_dose_uncertainty_estimate(
        dose, unc, mask, spacing_mm=(10.0, 10.0, 10.0), photons=True, correlation_length_mm=100.0
    )
    assert unc_ph == pytest.approx(0.05)


def test_roi_mean_chord_mm():
    # Sphere of radius r: mean chord 4r/3 (Cauchy formula)
    r = 20
    z, y, x = np.mgrid[-r - 1 : r + 2, -r - 1 : r + 2, -r - 1 : r + 2]
    sphere = (x**2 + y**2 + z**2 <= r**2).astype(np.uint8)
    assert roi_mean_chord_mm(sphere, (1.0, 1.0, 1.0)) == pytest.approx(4 * r / 3, rel=0.01)
    # Spacing scales the chord; sitk input reads the spacing
    img = sitk.GetImageFromArray(sphere)
    img.SetSpacing((2.0, 2.0, 2.0))
    assert roi_mean_chord_mm(img) == pytest.approx(8 * r / 3, rel=0.01)
    with pytest.raises(ValueError, match="ROI mask contains no voxels"):
        roi_mean_chord_mm(np.zeros((2, 2, 2)))
