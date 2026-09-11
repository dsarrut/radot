import json
import pytest
import numpy as np
from pathlib import Path
from click.testing import CliRunner
import SimpleITK as sitk
from scipy.integrate import quad

from radot.dose import (
    fit_exp_linear,
    decay_corrected_tac,
    triexpo_fit_jackson2013,
    triexpo_integrate_jackson2013,
    dose_jackson2013,
    compute_roi_dose_mtp,
    compute_roi_dose_file,
    compute_roi_dose,
)
from radot.cli import main
from radot.db import get_physical_half_life_h


def test_fit_exp_linear():
    x = np.array([1.0, 5.0])
    y = 100.0 * np.exp(-0.2 * x)
    A, k = fit_exp_linear(x, y)
    assert np.isclose(A, 100.0, rtol=1e-5)
    assert np.isclose(k, -0.2, rtol=1e-5)

    # In Jackson stripping, y[0] can be negative while y[1] is positive (0.01 * d_next)
    y_strip = np.array([-50.0, 10.0])
    A_strip, k_strip = fit_exp_linear(x, y_strip)
    assert A_strip < 0
    assert k_strip <= 0


def test_triexpo_fit_against_rpt_dosi_reference():
    """
    Compare triexpo_fit_jackson2013 against rpt_dosi reference data
    for liver, spleen, left_kidney, right_kidney.
    """
    ref_file = Path("/Users/dsarrut/src/py/rpt_dosi/rpt_dosi/tests/data/test003/activities_fit_ref.json")
    if not ref_file.exists():
        pytest.skip(f"Reference file not found: {ref_file}")

    with open(ref_file, "r") as f:
        ref_data = json.load(f)

    for organ, organ_ref in ref_data.items():
        times = np.array(organ_ref["times"], dtype=float)
        activities = np.array(organ_ref["activities"], dtype=float)

        fit = triexpo_fit_jackson2013(times, activities, as_dict=True)

        for param in ["A1", "k1", "A2", "k2", "A3", "k3"]:
            expected = organ_ref[param]
            actual = fit[param]
            assert np.isclose(actual, expected, rtol=1e-4, atol=1e-4), (
                f"Mismatch in organ {organ} for {param}: got {actual}, expected {expected}"
            )


def test_triexpo_integrate_against_quadrature():
    """
    Verify analytical integral triexpo_integrate_jackson2013 matches numerical quadrature.
    """
    params = {
        "A1": -606.9115,
        "k1": -1.3,
        "A2": 475.5603,
        "k2": -0.2572,
        "A3": 131.3512,
        "k3": -0.01904,
    }
    t_p_h = get_physical_half_life_h("Lu177")
    decay_constant_h = np.log(2) / t_p_h

    # Analytical integral in MBq.s
    analytical_tia_s = triexpo_integrate_jackson2013(params, decay_constant_h)

    # Numerical integration: int_0^infty sum(A_i * exp((k_i - lambda_p) * t)) dt in hours
    def physical_tac(t):
        val = 0.0
        for i in (1, 2, 3):
            val += params[f"A{i}"] * np.exp((params[f"k{i}"] - decay_constant_h) * t)
        return val

    numerical_tia_h, _ = quad(physical_tac, 0, np.inf)
    numerical_tia_s = numerical_tia_h * 3600.0

    assert np.isclose(analytical_tia_s, numerical_tia_s, rtol=1e-4)


def test_dose_jackson2013_calculation():
    """
    Test direct dose_jackson2013 calculation.
    """
    times_h = [4.5, 24.0, 120.0]
    activities_mbq = [250.0, 70.0, 10.0]
    s_val = 1.5e-6  # mGy / MBq.s
    mass_scaling = 1.1

    dose_gy, tia_s, teff_h, fit_params = dose_jackson2013(
        times_h=times_h,
        activities_mbq=activities_mbq,
        s_value=s_val,
        mass_scaling=mass_scaling,
        isotope="Lu177",
    )

    assert dose_gy > 0.0
    assert tia_s > 0.0
    assert teff_h > 0.0
    assert "A1" in fit_params
    assert "A3" in fit_params
    expected_dose = (tia_s * s_val) / mass_scaling / 1000.0
    assert np.isclose(dose_gy, expected_dose, rtol=1e-6)


def test_compute_roi_dose_mtp_synthetic(tmp_path):
    """
    Test full MTP pipeline with synthetic SimpleITK images and metadata sidecars.
    """
    # 1. Create synthetic label image with Liver (1) and Spleen (2)
    labels_arr = np.zeros((10, 10, 10), dtype=np.uint8)
    labels_arr[2:5, 2:5, 2:5] = 1  # Liver: 27 voxels
    labels_arr[6:8, 6:8, 6:8] = 2  # Spleen: 8 voxels

    labels_img = sitk.GetImageFromArray(labels_arr)
    labels_img.SetSpacing((2.0, 2.0, 2.0))
    labels_path = tmp_path / "labels.nii.gz"
    sitk.WriteImage(labels_img, str(labels_path))

    # Write labels sidecar JSON
    with open(tmp_path / "labels.json", "w") as f:
        json.dump({"1": "Liver", "2": "Spleen"}, f)

    # 2. Create 3 synthetic SPECT images at 5h, 24h, 120h post injection
    times = [5.0, 24.0, 120.0]
    dates = [
        "2024-01-01 15:00:00",
        "2024-01-02 10:00:00",
        "2024-01-06 10:00:00",
    ]
    spect_paths = []
    for i, (t_h, acq_date) in enumerate(zip(times, dates)):
        spect_arr = np.zeros((10, 10, 10), dtype=np.float32)
        spect_arr[labels_arr == 1] = 10000.0 * np.exp(-0.03 * t_h)
        spect_arr[labels_arr == 2] = 5000.0 * np.exp(-0.04 * t_h)

        spect_img = sitk.GetImageFromArray(spect_arr)
        spect_img.SetSpacing((2.0, 2.0, 2.0))
        img_path = tmp_path / f"spect_tp{i+1}.nii.gz"
        sitk.WriteImage(spect_img, str(img_path))

        with open(tmp_path / f"spect_tp{i+1}.json", "w") as f:
            json.dump(
                {
                    "unit": "bqml",
                    "modality": "spect",
                    "injection_date": "2024-01-01 10:00:00",
                    "acquisition_date": acq_date,
                    "activity": 5000.0,
                },
                f,
            )
        spect_paths.append(img_path)

    # Run compute_roi_dose_mtp
    results = compute_roi_dose_mtp(
        images=spect_paths,
        labels=labels_path,
        method="jackson2013",
        isotope="Lu177",
    )

    assert "Liver" in results
    assert "Spleen" in results
    liver = results["Liver"]
    assert liver["mean_dose_gy"] > 0.0
    assert liver["tia_mbq_s"] > 0.0
    assert "fit_params" in liver
    assert len(liver["times_h"]) == 3
    assert np.allclose(liver["times_h"], [5.0, 24.0, 120.0])

    # Also test dispatch from compute_roi_dose with list of images
    results_dispatch = compute_roi_dose(
        image=spect_paths,
        labels=labels_path,
        method="jackson2013",
        isotope="Lu177",
    )
    assert np.isclose(
        results_dispatch["Liver"]["mean_dose_gy"], liver["mean_dose_gy"]
    )


def test_cli_json_config_mtp(tmp_path):
    """
    Test radot dose CLI command with a multi-timepoint JSON config file.
    """
    labels_arr = np.zeros((8, 8, 8), dtype=np.uint8)
    labels_arr[2:5, 2:5, 2:5] = 1
    labels_img = sitk.GetImageFromArray(labels_arr)
    labels_img.SetSpacing((2.0, 2.0, 2.0))
    labels_path = tmp_path / "labels.nii.gz"
    sitk.WriteImage(labels_img, str(labels_path))
    with open(tmp_path / "labels.json", "w") as f:
        json.dump({"1": "Liver"}, f)

    times = [5.0, 24.0, 120.0]
    spect_names = []
    for i, t_h in enumerate(times):
        spect_arr = np.zeros((8, 8, 8), dtype=np.float32)
        spect_arr[labels_arr == 1] = 8000.0 * np.exp(-0.03 * t_h)
        spect_img = sitk.GetImageFromArray(spect_arr)
        spect_img.SetSpacing((2.0, 2.0, 2.0))
        img_name = f"spect_tp{i+1}.nii.gz"
        sitk.WriteImage(spect_img, str(tmp_path / img_name))
        with open(tmp_path / f"spect_tp{i+1}.json", "w") as f:
            json.dump(
                {
                    "unit": "bqml",
                    "modality": "spect",
                    "injection_date": "2024-01-01 10:00:00",
                    "acquisition_date": f"2024-01-{i+1:02d} 10:00:00",
                },
                f,
            )
        spect_names.append(img_name)

    # Create MTP JSON configuration file with explicit time_h
    config_data = {
        "method": "jackson2013",
        "isotope": "Lu177",
        "labels": "labels.nii.gz",
        "output": "dose_out.json",
        "timepoints": [
            {"image": spect_names[0], "time_h": 5.0},
            {"image": spect_names[1], "time_h": 24.0},
            {"image": spect_names[2], "time_h": 120.0},
        ],
    }
    config_path = tmp_path / "dose_jackson.json"
    with open(config_path, "w") as f:
        json.dump(config_data, f, indent=4)

    # Run via CLI
    runner = CliRunner()
    result = runner.invoke(main, ["dose", str(config_path)])
    assert result.exit_code == 0, result.output

    out_file = tmp_path / "dose_out.json"
    assert out_file.exists()
    with open(out_file, "r") as f:
        data = json.load(f)
    assert "Liver" in data
    assert data["Liver"]["mean_dose_gy"] > 0.0
    assert "tia_mbq_h" in data["Liver"]


def test_cli_json_config_mtp_per_tp_labels_override(tmp_path):
    """
    Test that per-timepoint labels override global labels when both are provided.
    We create a global label where Liver is in a wrong/empty region,
    and per-timepoint labels where Liver is at the correct position.
    """
    # 1. Global labels: Liver (1) placed at (0, 0, 0) where activity is ZERO
    global_labels_arr = np.zeros((8, 8, 8), dtype=np.uint8)
    global_labels_arr[0, 0, 0] = 1
    global_lbl_img = sitk.GetImageFromArray(global_labels_arr)
    global_lbl_img.SetSpacing((2.0, 2.0, 2.0))
    global_lbl_path = tmp_path / "global_labels.nii.gz"
    sitk.WriteImage(global_lbl_img, str(global_lbl_path))
    with open(tmp_path / "global_labels.json", "w") as f:
        json.dump({"1": "Liver"}, f)

    # 2. Per-timepoint labels: Liver placed where activity actually is
    times = [5.0, 24.0, 120.0]
    spect_names = []
    tp_label_names = []

    for i, t_h in enumerate(times):
        # Shift organ position across timepoints to simulate patient repositioning
        shift = i  # 0, 1, 2
        tp_arr = np.zeros((8, 8, 8), dtype=np.uint8)
        tp_arr[2 + shift:5 + shift, 2:5, 2:5] = 1
        tp_img = sitk.GetImageFromArray(tp_arr)
        tp_img.SetSpacing((2.0, 2.0, 2.0))
        lbl_name = f"labels_tp{i+1}.nii.gz"
        sitk.WriteImage(tp_img, str(tmp_path / lbl_name))
        with open(tmp_path / f"labels_tp{i+1}.json", "w") as f:
            json.dump({"1": "Liver"}, f)
        tp_label_names.append(lbl_name)

        # SPECT image with activity in the shifted organ
        spect_arr = np.zeros((8, 8, 8), dtype=np.float32)
        spect_arr[tp_arr == 1] = 8000.0 * np.exp(-0.03 * t_h)
        spect_img = sitk.GetImageFromArray(spect_arr)
        spect_img.SetSpacing((2.0, 2.0, 2.0))
        img_name = f"spect_tp{i+1}.nii.gz"
        sitk.WriteImage(spect_img, str(tmp_path / img_name))
        with open(tmp_path / f"spect_tp{i+1}.json", "w") as f:
            json.dump(
                {
                    "unit": "bqml",
                    "modality": "spect",
                    "injection_date": "2024-01-01 10:00:00",
                    "acquisition_date": f"2024-01-{i+1:02d} 10:00:00",
                },
                f,
            )
        spect_names.append(img_name)

    # JSON config providing BOTH global labels AND per-timepoint labels:
    # Per-timepoint labels should override global_labels!
    config_data = {
        "method": "jackson2013",
        "isotope": "Lu177",
        "labels": "global_labels.nii.gz",
        "output": "dose_override.json",
        "timepoints": [
            {"image": spect_names[0], "labels": tp_label_names[0], "time_h": 5.0},
            {"image": spect_names[1], "labels": tp_label_names[1], "time_h": 24.0},
            {"image": spect_names[2], "labels": tp_label_names[2], "time_h": 120.0},
        ],
    }
    config_path = tmp_path / "dose_override.json"
    with open(config_path, "w") as f:
        json.dump(config_data, f, indent=4)

    runner = CliRunner()
    result = runner.invoke(main, ["dose", str(config_path)])
    assert result.exit_code == 0, result.output

    out_file = tmp_path / "dose_override.json"
    assert out_file.exists()
    with open(out_file, "r") as f:
        data = json.load(f)

    # If global labels had been used, Liver activity would be 0 and dose would be 0
    # With per-timepoint labels, dose is positive and non-zero!
    assert "Liver" in data
    assert data["Liver"]["mean_dose_gy"] > 0.0
    assert data["Liver"]["activities_mbq"][0] > 0.0
    assert data["Liver"]["activities_mbq"][1] > 0.0
    assert data["Liver"]["activities_mbq"][2] > 0.0
