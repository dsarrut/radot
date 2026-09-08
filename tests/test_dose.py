import json
from pathlib import Path
import numpy as np
import pytest
import SimpleITK as sitk
from click.testing import CliRunner
from radot.cli import main
from radot.dose import compute_roi_dose


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


def test_madsen2018_opengate_p142_comparison():
    # pytest -s tests/test_dose.py -k "test_madsen2018_opengate_p142_comparison" -v
    bench_path = (
        Path(__file__).resolve().parent / "data" / "p142_opengate_benchmark.json"
    )
    if not bench_path.exists():
        pytest.skip("p142 OpenGate benchmark JSON not found.")

    with open(bench_path) as f:
        bench = json.load(f)

    p142_raw_dir = Path(__file__).resolve().parent.parent / "data" / "p142"
    gate_dose_path = p142_raw_dir / "output_vrt" / "output_dose.mhd"
    liver_path = p142_raw_dir / "liver_4mm_cropped.nii.gz"
    ct_path = p142_raw_dir / "ct_4mm_cropped.nii.gz"

    total_activity_bq = float(bench["total_activity_bq"])
    isotope = bench.get("radionuclide", "Lu177")
    teff_h = float(bench["teff_h"])
    time_from_injection_h = float(bench.get("time_from_injection_h", 0.0))
    roi_name = bench["roi_name"]
    phantom = bench.get("phantom", "male")

    # If full patient 3D images exist locally in data/p142, test with actual images
    if gate_dose_path.exists() and liver_path.exists() and ct_path.exists():
        liver_img = sitk.ReadImage(str(liver_path))
        mask = sitk.GetArrayViewFromImage(liver_img) > 0
        num_voxels = int(mask.sum())

        activity_img = sitk.Cast(liver_img, sitk.sitkFloat32) * (
            total_activity_bq / num_voxels
        )
        activity_img.SetMetaData("unit", "Bq")
        activity_img.SetMetaData("time_from_injection", str(time_from_injection_h))
        labels_img = sitk.Cast(liver_img > 0, sitk.sitkUInt16)
        ct_img = sitk.ReadImage(str(ct_path))

        gate_img = sitk.ReadImage(str(gate_dose_path))
        gate_arr = sitk.GetArrayViewFromImage(gate_img)
        gate_mean_dose_rate_gy_s = float(gate_arr[mask].mean())
    else:
        # Otherwise, test using the benchmark parameters without needing heavy image binaries
        roi_volume_ml = float(bench["roi_volume_ml"])
        roi_mass_g = float(bench["roi_mass_g"])
        gate_mean_dose_rate_gy_s = float(bench["opengate"]["mean_dose_rate_gy_s"])

        s = (roi_volume_ml * 1000.0) ** (1.0 / 3.0)
        activity_img = sitk.Image([1, 1, 1], sitk.sitkFloat32)
        activity_img.SetSpacing((s, s, s))
        activity_img.SetMetaData("unit", "Bq")
        activity_img.SetMetaData("time_from_injection", str(time_from_injection_h))
        activity_img[0, 0, 0] = total_activity_bq

        labels_img = sitk.Image([1, 1, 1], sitk.sitkUInt16)
        labels_img.SetSpacing((s, s, s))
        labels_img[0, 0, 0] = 1

        ct_img = sitk.Image([1, 1, 1], sitk.sitkInt16)
        ct_img.SetSpacing((s, s, s))
        ct_val = int(round((roi_mass_g / roi_volume_ml - 1.0) * 1000.0))
        ct_img[0, 0, 0] = ct_val

    results = compute_roi_dose(
        image=activity_img,
        labels=labels_img,
        method="madsen2018",
        isotope=isotope,
        ct_image=ct_img,
        labels_map={1: roi_name},
        teff_params={roi_name: teff_h},
        method_params={"phantom": phantom},
    )
    dose_madsen_gy = results[roi_name]["mean_dose_gy"]

    # Time-integrate OpenGate dose rate with mono-exponential clearance
    k = np.log(2) / (teff_h * 3600.0)
    gate_integrated_dose_gy = gate_mean_dose_rate_gy_s / k

    rel_diff = abs(dose_madsen_gy - gate_integrated_dose_gy) / gate_integrated_dose_gy
    print(f"Madsen 2018 mean dose: {dose_madsen_gy:.6e} Gy")
    print(f"OpenGate integrated dose: {gate_integrated_dose_gy:.6e} Gy")
    print(f"Relative difference: {rel_diff * 100:.2f}%")

    assert (
        rel_diff < 0.05
    ), f"Madsen vs OpenGate difference too large: {rel_diff * 100:.2f}%"
    assert dose_madsen_gy == pytest.approx(
        bench["madsen2018"]["mean_dose_gy"], rel=1e-3
    )


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
