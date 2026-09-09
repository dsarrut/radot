import json
from pathlib import Path
import numpy as np
import pytest
import SimpleITK as sitk
from click.testing import CliRunner
from radot.cli import main
from radot.dose import (
    compute_roi_dose,
    dose_hanscheid2018,
    dose_madsen2018,
    roi_dose_uncertainty_bounds,
    roi_dose_uncertainty_estimate,
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


def test_madsen2018_self_dose():
    # pytest -s tests/test_dose.py -k "test_madsen2018_self_dose" -v
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

        gate_unc_path = p142_raw_dir / "output_vrt" / "output_dose_uncertainty.mhd"
        if gate_unc_path.exists():
            gate_unc_img = sitk.ReadImage(str(gate_unc_path))
            unc_min, unc_max = roi_dose_uncertainty_bounds(gate_img, gate_unc_img, mask)
            unc_est = roi_dose_uncertainty_estimate(
                gate_img, gate_unc_img, mask, isotope=isotope
            )
        else:
            unc_min, unc_max, unc_est = None, None, None
    else:
        # Otherwise, test using the benchmark parameters without needing heavy image binaries
        roi_volume_ml = float(bench["roi_volume_ml"])
        roi_mass_g = float(bench["roi_mass_g"])
        gate_mean_dose_rate_gy_s = float(bench["opengate"]["mean_dose_rate_gy_s"])
        op_data = bench.get("opengate", {})
        unc_min = op_data.get("uncertainty_min_fraction")
        unc_max = op_data.get(
            "uncertainty_max_fraction", op_data.get("mean_voxel_uncertainty_fraction")
        )
        unc_est = op_data.get("uncertainty_estimate_fraction")

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
    if unc_est is not None and unc_min is not None and unc_max is not None:
        print(
            f"OpenGate MC uncertainty (liver): {unc_est * 100:.2f}% (bounds: [{unc_min * 100:.3f}%, {unc_max * 100:.2f}%])"
        )
    print(f"Relative difference: {rel_diff * 100:.2f}%  (<5% expected)")

    assert (
        rel_diff < 0.05
    ), f"Madsen vs OpenGate difference too large: {rel_diff * 100:.2f}%"


def test_hanscheid2017_self_dose():
    # pytest -s tests/test_dose.py -k "test_hanscheid2017_self_dose" -v
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

    total_activity_bq = float(bench["total_activity_bq"])
    isotope = bench.get("radionuclide", "Lu177")
    teff_h = float(bench["teff_h"])
    time_from_injection_h = float(bench.get("time_from_injection_h", 0.0))
    roi_name = bench["roi_name"]

    # If full patient 3D images exist locally in data/p142, test with actual images
    if gate_dose_path.exists() and liver_path.exists():
        liver_img = sitk.ReadImage(str(liver_path))
        mask = sitk.GetArrayViewFromImage(liver_img) > 0
        num_voxels = int(mask.sum())

        activity_img = sitk.Cast(liver_img, sitk.sitkFloat32) * (
            total_activity_bq / num_voxels
        )
        activity_img.SetMetaData("unit", "Bq")
        activity_img.SetMetaData("time_from_injection", str(time_from_injection_h))
        labels_img = sitk.Cast(liver_img > 0, sitk.sitkUInt16)

        gate_img = sitk.ReadImage(str(gate_dose_path))
        gate_arr = sitk.GetArrayViewFromImage(gate_img)
        gate_mean_dose_rate_gy_s = float(gate_arr[mask].mean())

        gate_unc_path = p142_raw_dir / "output_vrt" / "output_dose_uncertainty.mhd"
        if gate_unc_path.exists():
            gate_unc_img = sitk.ReadImage(str(gate_unc_path))
            unc_min, unc_max = roi_dose_uncertainty_bounds(gate_img, gate_unc_img, mask)
            unc_est = roi_dose_uncertainty_estimate(
                gate_img, gate_unc_img, mask, isotope=isotope
            )
        else:
            unc_min, unc_max, unc_est = None, None, None
    else:
        # Otherwise, test using the benchmark parameters without needing heavy image binaries
        roi_volume_ml = float(bench["roi_volume_ml"])
        gate_mean_dose_rate_gy_s = float(bench["opengate"]["mean_dose_rate_gy_s"])
        op_data = bench.get("opengate", {})
        unc_min = op_data.get("uncertainty_min_fraction")
        unc_max = op_data.get(
            "uncertainty_max_fraction", op_data.get("mean_voxel_uncertainty_fraction")
        )
        unc_est = op_data.get("uncertainty_estimate_fraction")

        s = (roi_volume_ml * 1000.0) ** (1.0 / 3.0)
        activity_img = sitk.Image([1, 1, 1], sitk.sitkFloat32)
        activity_img.SetSpacing((s, s, s))
        activity_img.SetMetaData("unit", "Bq")
        activity_img.SetMetaData("time_from_injection", str(time_from_injection_h))
        activity_img[0, 0, 0] = total_activity_bq

        labels_img = sitk.Image([1, 1, 1], sitk.sitkUInt16)
        labels_img.SetSpacing((s, s, s))
        labels_img[0, 0, 0] = 1

    results = compute_roi_dose(
        image=activity_img,
        labels=labels_img,
        method="hanscheid2017",
        isotope=isotope,
        labels_map={1: roi_name},
        teff_params={roi_name: teff_h},
    )
    dose_hanscheid_gy = results[roi_name]["mean_dose_gy"]

    # Time-integrate OpenGate dose rate with mono-exponential clearance
    k = np.log(2) / (teff_h * 3600.0)
    gate_integrated_dose_gy = gate_mean_dose_rate_gy_s / k

    rel_diff = (
        abs(dose_hanscheid_gy - gate_integrated_dose_gy) / gate_integrated_dose_gy
    )
    print(f"\nH\u00e4nscheid 2017 mean dose: {dose_hanscheid_gy:.6e} Gy")
    print(f"OpenGate integrated dose: {gate_integrated_dose_gy:.6e} Gy")
    if unc_est is not None and unc_min is not None and unc_max is not None:
        print(
            f"OpenGate MC uncertainty (liver): {unc_est * 100:.2f}% (bounds: [{unc_min * 100:.3f}%, {unc_max * 100:.2f}%])"
        )
    print(f"Relative difference: {rel_diff * 100:.2f}%  (<15% expected)")

    # Hanscheid 2017 assumes tissue density ~1.0 g/mL (patient liver CT density is ~1.08 g/mL)
    # and an empirical factor of 0.125, leading to expected ~9.5% relative difference
    assert (
        rel_diff < 0.15
    ), f"Hanscheid vs OpenGate difference too large: {rel_diff * 100:.2f}%"


def test_hanscheid2018_self_dose():
    # pytest -s tests/test_dose.py -k "test_hanscheid2018_self_dose" -v
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
    # In Hanscheid 2018, Teff = 2 * t, so acquisition time t = Teff / 2 corresponds to Teff = 45h
    time_from_injection_h = teff_h / 2.0
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

        gate_unc_path = p142_raw_dir / "output_vrt" / "output_dose_uncertainty.mhd"
        if gate_unc_path.exists():
            gate_unc_img = sitk.ReadImage(str(gate_unc_path))
            unc_min, unc_max = roi_dose_uncertainty_bounds(gate_img, gate_unc_img, mask)
            unc_est = roi_dose_uncertainty_estimate(
                gate_img, gate_unc_img, mask, isotope=isotope
            )
        else:
            unc_min, unc_max, unc_est = None, None, None
    else:
        # Otherwise, test using the benchmark parameters without needing heavy image binaries
        roi_volume_ml = float(bench["roi_volume_ml"])
        roi_mass_g = float(bench["roi_mass_g"])
        gate_mean_dose_rate_gy_s = float(bench["opengate"]["mean_dose_rate_gy_s"])
        op_data = bench.get("opengate", {})
        unc_min = op_data.get("uncertainty_min_fraction")
        unc_max = op_data.get(
            "uncertainty_max_fraction", op_data.get("mean_voxel_uncertainty_fraction")
        )
        unc_est = op_data.get("uncertainty_estimate_fraction")

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
        method="hanscheid2018",
        isotope=isotope,
        ct_image=ct_img,
        labels_map={1: roi_name},
        method_params={"phantom": phantom},
    )
    dose_hanscheid_gy = results[roi_name]["mean_dose_gy"]

    # Time-integrate OpenGate dose rate with mono-exponential clearance
    k = np.log(2) / (teff_h * 3600.0)
    gate_integrated_dose_gy = gate_mean_dose_rate_gy_s / k

    rel_diff = (
        abs(dose_hanscheid_gy - gate_integrated_dose_gy) / gate_integrated_dose_gy
    )
    print(f"\nH\u00e4nscheid 2018 mean dose: {dose_hanscheid_gy:.6e} Gy")
    print(f"OpenGate integrated dose: {gate_integrated_dose_gy:.6e} Gy")
    if unc_est is not None and unc_min is not None and unc_max is not None:
        print(
            f"OpenGate MC uncertainty (liver): {unc_est * 100:.2f}% (bounds: [{unc_min * 100:.3f}%, {unc_max * 100:.2f}%])"
        )
    print(f"Relative difference: {rel_diff * 100:.2f}%  (<5% expected)")

    assert (
        rel_diff < 0.05
    ), f"Hanscheid 2018 vs OpenGate difference too large: {rel_diff * 100:.2f}% (<5% expected)"


def test_ldm_self_dose():
    # pytest -s tests/test_dose.py -k "test_ldm_self_dose" -v
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
    roi_name = bench["roi_name"]

    # If full patient 3D images exist locally in data/p142, test with actual images
    if gate_dose_path.exists() and liver_path.exists() and ct_path.exists():
        liver_img = sitk.ReadImage(str(liver_path))
        mask = sitk.GetArrayViewFromImage(liver_img) > 0
        num_voxels = int(mask.sum())

        activity_img = sitk.Cast(liver_img, sitk.sitkFloat32) * (
            total_activity_bq / num_voxels
        )
        activity_img.SetMetaData("unit", "Bq")
        activity_img.SetMetaData("time_from_injection", "0.0")
        labels_img = sitk.Cast(liver_img > 0, sitk.sitkUInt16)
        ct_img = sitk.ReadImage(str(ct_path))

        gate_img = sitk.ReadImage(str(gate_dose_path))
        gate_arr = sitk.GetArrayViewFromImage(gate_img)
        gate_mean_dose_rate_gy_s = float(gate_arr[mask].mean())

        gate_unc_path = p142_raw_dir / "output_vrt" / "output_dose_uncertainty.mhd"
        if gate_unc_path.exists():
            gate_unc_img = sitk.ReadImage(str(gate_unc_path))
            unc_min, unc_max = roi_dose_uncertainty_bounds(gate_img, gate_unc_img, mask)
            unc_est = roi_dose_uncertainty_estimate(
                gate_img, gate_unc_img, mask, isotope=isotope
            )
        else:
            unc_min, unc_max, unc_est = None, None, None
    else:
        # Otherwise, test using the benchmark parameters without needing heavy image binaries
        roi_volume_ml = float(bench["roi_volume_ml"])
        roi_mass_g = float(bench["roi_mass_g"])
        gate_mean_dose_rate_gy_s = float(bench["opengate"]["mean_dose_rate_gy_s"])
        op_data = bench.get("opengate", {})
        unc_min = op_data.get("uncertainty_min_fraction")
        unc_max = op_data.get(
            "uncertainty_max_fraction", op_data.get("mean_voxel_uncertainty_fraction")
        )
        unc_est = op_data.get("uncertainty_estimate_fraction")

        s = (roi_volume_ml * 1000.0) ** (1.0 / 3.0)
        activity_img = sitk.Image([1, 1, 1], sitk.sitkFloat32)
        activity_img.SetSpacing((s, s, s))
        activity_img.SetMetaData("unit", "Bq")
        activity_img.SetMetaData("time_from_injection", "0.0")
        activity_img[0, 0, 0] = total_activity_bq

        labels_img = sitk.Image([1, 1, 1], sitk.sitkUInt16)
        labels_img.SetSpacing((s, s, s))
        labels_img[0, 0, 0] = 1

        ct_img = sitk.Image([1, 1, 1], sitk.sitkInt16)
        ct_img.SetSpacing((s, s, s))
        ct_val = int(round((roi_mass_g / roi_volume_ml - 1.0) * 1000.0))
        ct_img[0, 0, 0] = ct_val

    results, _ = compute_roi_dose(
        image=activity_img,
        labels=labels_img,
        method="ldm",
        isotope=isotope,
        ct_image=ct_img,
        labels_map={1: roi_name},
        dose_rate=True,
    )
    ldm_mean_dose_rate_gy_s = results[roi_name]["mean_dose_rate_gyh"] / 3600.0

    # Time-integrate dose rate with mono-exponential clearance
    k = np.log(2) / (teff_h * 3600.0)
    ldm_integrated_dose_gy = ldm_mean_dose_rate_gy_s / k
    gate_integrated_dose_gy = gate_mean_dose_rate_gy_s / k

    rel_diff = (
        abs(ldm_integrated_dose_gy - gate_integrated_dose_gy) / gate_integrated_dose_gy
    )
    print(f"\nLDM mean dose rate:        {ldm_mean_dose_rate_gy_s:.6e} Gy/s")
    print(f"OpenGate mean dose rate:    {gate_mean_dose_rate_gy_s:.6e} Gy/s")
    print(f"LDM integrated dose:        {ldm_integrated_dose_gy:.6e} Gy")
    print(f"OpenGate integrated dose:   {gate_integrated_dose_gy:.6e} Gy")
    if unc_est is not None and unc_min is not None and unc_max is not None:
        print(
            f"OpenGate MC uncertainty (liver): {unc_est * 100:.2f}% (bounds: [{unc_min * 100:.3f}%, {unc_max * 100:.2f}%])"
        )
    print(f"Relative difference:        {rel_diff * 100:.2f}%  (<1% expected)")

    # LDM non-penetrating deposition matches Monte Carlo electron transport closely (<1%)
    assert (
        rel_diff < 0.01
    ), f"LDM vs OpenGate difference too large: {rel_diff * 100:.2f}% (<1% expected)"


def test_y90_self_dose():
    # pytest -s tests/test_dose.py -k "test_y90_self_dose" -v
    bench_path = (
        Path(__file__).resolve().parent / "data" / "p142_y90_opengate_benchmark.json"
    )
    if not bench_path.exists():
        pytest.skip("p142 Y90 OpenGate benchmark JSON not found.")

    with open(bench_path) as f:
        bench = json.load(f)

    p142_raw_dir = Path(__file__).resolve().parent.parent / "data" / "p142"
    gate_dose_path = p142_raw_dir / "output_y90_vrt" / "output_dose.mhd"
    liver_path = p142_raw_dir / "liver_4mm_cropped.nii.gz"
    ct_path = p142_raw_dir / "ct_4mm_cropped.nii.gz"

    total_activity_bq = float(bench["total_activity_bq"])
    isotope = "Y90"
    teff_h = float(bench.get("teff_h", 64.1))
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
        activity_img.SetMetaData("time_from_injection", "0.0")
        labels_img = sitk.Cast(liver_img > 0, sitk.sitkUInt16)
        ct_img = sitk.ReadImage(str(ct_path))

        gate_img = sitk.ReadImage(str(gate_dose_path))
        gate_arr = sitk.GetArrayViewFromImage(gate_img)
        gate_mean_dose_rate_gy_s = float(gate_arr[mask].mean())

        gate_unc_path = p142_raw_dir / "output_y90_vrt" / "output_dose_uncertainty.mhd"
        if gate_unc_path.exists():
            gate_unc_img = sitk.ReadImage(str(gate_unc_path))
            unc_min, unc_max = roi_dose_uncertainty_bounds(gate_img, gate_unc_img, mask)
            unc_est = roi_dose_uncertainty_estimate(
                gate_img, gate_unc_img, mask, isotope=isotope
            )
        else:
            unc_min, unc_max, unc_est = None, None, None
    else:
        # Otherwise, test using the benchmark parameters without needing heavy image binaries
        roi_volume_ml = float(bench["roi_volume_ml"])
        roi_mass_g = float(bench["roi_mass_g"])
        gate_mean_dose_rate_gy_s = float(bench["opengate"]["mean_dose_rate_gy_s"])
        op_data = bench.get("opengate", {})
        unc_min = op_data.get("uncertainty_min_fraction")
        unc_max = op_data.get(
            "uncertainty_max_fraction", op_data.get("mean_voxel_uncertainty_fraction")
        )
        unc_est = op_data.get("uncertainty_estimate_fraction")

        s = (roi_volume_ml * 1000.0) ** (1.0 / 3.0)
        activity_img = sitk.Image([1, 1, 1], sitk.sitkFloat32)
        activity_img.SetSpacing((s, s, s))
        activity_img.SetMetaData("unit", "Bq")
        activity_img.SetMetaData("time_from_injection", "0.0")
        activity_img[0, 0, 0] = total_activity_bq

        labels_img = sitk.Image([1, 1, 1], sitk.sitkUInt16)
        labels_img.SetSpacing((s, s, s))
        labels_img[0, 0, 0] = 1

        ct_img = sitk.Image([1, 1, 1], sitk.sitkInt16)
        ct_img.SetSpacing((s, s, s))
        ct_val = int(round((roi_mass_g / roi_volume_ml - 1.0) * 1000.0))
        ct_img[0, 0, 0] = ct_val

    # 1. Madsen 2018 calculation
    results_madsen = compute_roi_dose(
        image=activity_img,
        labels=labels_img,
        method="madsen2018",
        isotope=isotope,
        ct_image=ct_img,
        labels_map={1: roi_name},
        teff_params={roi_name: teff_h},
        method_params={"phantom": phantom},
    )
    dose_madsen_gy = results_madsen[roi_name]["mean_dose_gy"]

    # 2. LDM calculation
    results_ldm, _ = compute_roi_dose(
        image=activity_img,
        labels=labels_img,
        method="ldm",
        isotope=isotope,
        ct_image=ct_img,
        labels_map={1: roi_name},
        dose_rate=True,
    )
    ldm_mean_dose_rate_gy_s = results_ldm[roi_name]["mean_dose_rate_gyh"] / 3600.0

    k = np.log(2) / (teff_h * 3600.0)
    ldm_integrated_dose_gy = ldm_mean_dose_rate_gy_s / k
    gate_integrated_dose_gy = gate_mean_dose_rate_gy_s / k

    rel_diff_madsen = (
        abs(dose_madsen_gy - gate_integrated_dose_gy) / gate_integrated_dose_gy
    )
    rel_diff_ldm = (
        abs(ldm_integrated_dose_gy - gate_integrated_dose_gy) / gate_integrated_dose_gy
    )

    print(f"\nY90 OpenGate mean dose rate:    {gate_mean_dose_rate_gy_s:.6e} Gy/s")
    print(f"Y90 OpenGate integrated dose:   {gate_integrated_dose_gy:.6e} Gy")
    print(
        f"Y90 Madsen 2018 mean dose:      {dose_madsen_gy:.6e} Gy (diff: {rel_diff_madsen * 100:.2f}%)"
    )
    print(
        f"Y90 LDM integrated dose:        {ldm_integrated_dose_gy:.6e} Gy (diff: {rel_diff_ldm * 100:.2f}%)"
    )
    if unc_est is not None and unc_min is not None and unc_max is not None:
        print(
            f"OpenGate MC uncertainty (liver): {unc_est * 100:.2f}% (bounds: [{unc_min * 100:.3f}%, {unc_max * 100:.2f}%])"
        )

    # For Y90:
    # LDM matches within 3% (slight edge effect due to ~11mm range)
    # Madsen 2018 matches within 10%
    assert (
        rel_diff_ldm < 0.05
    ), f"Y90 LDM vs OpenGate difference too large: {rel_diff_ldm * 100:.2f}% (<5% expected)"
    assert (
        rel_diff_madsen < 0.10
    ), f"Y90 Madsen vs OpenGate difference too large: {rel_diff_madsen * 100:.2f}% (<10% expected)"


def test_i131_self_dose():
    # pytest -s tests/test_dose.py -k "test_i131_self_dose" -v
    bench_path = (
        Path(__file__).resolve().parent / "data" / "p142_i131_opengate_benchmark.json"
    )
    if not bench_path.exists():
        pytest.skip("p142 I131 OpenGate benchmark JSON not found.")

    with open(bench_path) as f:
        bench = json.load(f)

    p142_raw_dir = Path(__file__).resolve().parent.parent / "data" / "p142"
    gate_dose_path = p142_raw_dir / "output_i131_vrt" / "output_dose.mhd"
    liver_path = p142_raw_dir / "liver_4mm_cropped.nii.gz"
    ct_path = p142_raw_dir / "ct_4mm_cropped.nii.gz"

    total_activity_bq = float(bench["total_activity_bq"])
    isotope = "I131"
    teff_h = float(bench.get("teff_h", 192.4968))
    roi_name = bench["roi_name"]
    phantom = bench.get("phantom", "male")

    if gate_dose_path.exists() and liver_path.exists() and ct_path.exists():
        liver_img = sitk.ReadImage(str(liver_path))
        mask = sitk.GetArrayViewFromImage(liver_img) > 0
        num_voxels = int(mask.sum())

        activity_img = sitk.Cast(liver_img, sitk.sitkFloat32) * (
            total_activity_bq / num_voxels
        )
        activity_img.SetMetaData("unit", "Bq")
        activity_img.SetMetaData("time_from_injection", "0.0")
        labels_img = sitk.Cast(liver_img > 0, sitk.sitkUInt16)
        ct_img = sitk.ReadImage(str(ct_path))

        gate_img = sitk.ReadImage(str(gate_dose_path))
        gate_arr = sitk.GetArrayViewFromImage(gate_img)
        gate_mean_dose_rate_gy_s = float(gate_arr[mask].mean())

        gate_unc_path = p142_raw_dir / "output_i131_vrt" / "output_dose_uncertainty.mhd"
        if gate_unc_path.exists():
            gate_unc_img = sitk.ReadImage(str(gate_unc_path))
            unc_min, unc_max = roi_dose_uncertainty_bounds(gate_img, gate_unc_img, mask)
            unc_est = roi_dose_uncertainty_estimate(
                gate_img, gate_unc_img, mask, isotope=isotope
            )
        else:
            unc_min, unc_max, unc_est = None, None, None
    else:
        roi_volume_ml = float(bench["roi_volume_ml"])
        roi_mass_g = float(bench["roi_mass_g"])
        gate_mean_dose_rate_gy_s = float(bench["opengate"]["mean_dose_rate_gy_s"])
        op_data = bench.get("opengate", {})
        unc_min = op_data.get("uncertainty_min_fraction")
        unc_max = op_data.get(
            "uncertainty_max_fraction", op_data.get("mean_voxel_uncertainty_fraction")
        )
        unc_est = op_data.get("uncertainty_estimate_fraction")

        s = (roi_volume_ml * 1000.0) ** (1.0 / 3.0)
        activity_img = sitk.Image([1, 1, 1], sitk.sitkFloat32)
        activity_img.SetSpacing((s, s, s))
        activity_img.SetMetaData("unit", "Bq")
        activity_img.SetMetaData("time_from_injection", "0.0")
        activity_img[0, 0, 0] = total_activity_bq

        labels_img = sitk.Image([1, 1, 1], sitk.sitkUInt16)
        labels_img.SetSpacing((s, s, s))
        labels_img[0, 0, 0] = 1

        ct_img = sitk.Image([1, 1, 1], sitk.sitkInt16)
        ct_img.SetSpacing((s, s, s))
        ct_val = int(round((roi_mass_g / roi_volume_ml - 1.0) * 1000.0))
        ct_img[0, 0, 0] = ct_val

    # 1. Madsen 2018 calculation
    results_madsen = compute_roi_dose(
        image=activity_img,
        labels=labels_img,
        method="madsen2018",
        isotope=isotope,
        ct_image=ct_img,
        labels_map={1: roi_name},
        teff_params={roi_name: teff_h},
        method_params={"phantom": phantom},
    )
    dose_madsen_gy = results_madsen[roi_name]["mean_dose_gy"]

    # 2. LDM calculation
    results_ldm, _ = compute_roi_dose(
        image=activity_img,
        labels=labels_img,
        method="ldm",
        isotope=isotope,
        ct_image=ct_img,
        labels_map={1: roi_name},
        dose_rate=True,
    )
    ldm_mean_dose_rate_gy_s = results_ldm[roi_name]["mean_dose_rate_gyh"] / 3600.0

    k = np.log(2) / (teff_h * 3600.0)
    ldm_integrated_dose_gy = ldm_mean_dose_rate_gy_s / k
    gate_integrated_dose_gy = gate_mean_dose_rate_gy_s / k

    rel_diff_madsen = (
        abs(dose_madsen_gy - gate_integrated_dose_gy) / gate_integrated_dose_gy
    )
    rel_diff_ldm = (
        abs(ldm_integrated_dose_gy - gate_integrated_dose_gy) / gate_integrated_dose_gy
    )

    print(f"\nI131 OpenGate mean dose rate:    {gate_mean_dose_rate_gy_s:.6e} Gy/s")
    print(f"I131 OpenGate integrated dose:   {gate_integrated_dose_gy:.6e} Gy")
    print(
        f"I131 Madsen 2018 mean dose:      {dose_madsen_gy:.6e} Gy (diff: {rel_diff_madsen * 100:.2f}%)"
    )
    print(
        f"I131 LDM integrated dose:        {ldm_integrated_dose_gy:.6e} Gy (diff: {rel_diff_ldm * 100:.2f}%)"
    )
    if unc_est is not None and unc_min is not None and unc_max is not None:
        print(
            f"OpenGate MC uncertainty (liver): {unc_est * 100:.2f}% (bounds: [{unc_min * 100:.3f}%, {unc_max * 100:.2f}%])"
        )

    # Madsen 2018 accounts for both beta and photon self-dose via OpenDose S-values (diff < 2%)
    assert (
        rel_diff_madsen < 0.05
    ), f"I131 Madsen vs OpenGate difference too large: {rel_diff_madsen * 100:.2f}% (<5% expected)"
    # LDM only accounts for local non-penetrating electrons (Delta_local = 0.192 MeV),
    # ignoring ~25% self-dose contributed by 364 keV photons in the large liver volume.
    assert (
        0.15 < rel_diff_ldm < 0.35
    ), f"I131 LDM vs OpenGate difference outside expected photon self-dose range: {rel_diff_ldm * 100:.2f}%"


def test_ac225_self_dose():
    # pytest -s tests/test_dose.py -k "test_ac225_self_dose" -v
    bench_path = (
        Path(__file__).resolve().parent / "data" / "p142_ac225_opengate_benchmark.json"
    )
    if not bench_path.exists():
        pytest.skip("p142 Ac225 OpenGate benchmark JSON not found.")

    with open(bench_path) as f:
        bench = json.load(f)

    p142_raw_dir = Path(__file__).resolve().parent.parent / "data" / "p142"
    gate_dose_path = p142_raw_dir / "output_ac225_vrt" / "output_dose.mhd"
    liver_path = p142_raw_dir / "liver_4mm_cropped.nii.gz"
    ct_path = p142_raw_dir / "ct_4mm_cropped.nii.gz"

    total_activity_bq = float(bench["total_activity_bq"])
    isotope = "Ac225"
    teff_h = float(bench.get("teff_h", 240.0))
    roi_name = bench["roi_name"]
    phantom = bench.get("phantom", "male")

    if gate_dose_path.exists() and liver_path.exists() and ct_path.exists():
        liver_img = sitk.ReadImage(str(liver_path))
        mask = sitk.GetArrayViewFromImage(liver_img) > 0
        num_voxels = int(mask.sum())

        activity_img = sitk.Cast(liver_img, sitk.sitkFloat32) * (
            total_activity_bq / num_voxels
        )
        activity_img.SetMetaData("unit", "Bq")
        activity_img.SetMetaData("time_from_injection", "0.0")
        labels_img = sitk.Cast(liver_img > 0, sitk.sitkUInt16)
        ct_img = sitk.ReadImage(str(ct_path))

        gate_img = sitk.ReadImage(str(gate_dose_path))
        gate_arr = sitk.GetArrayViewFromImage(gate_img)
        gate_mean_dose_rate_gy_s = float(gate_arr[mask].mean())

        gate_unc_path = (
            p142_raw_dir / "output_ac225_vrt" / "output_dose_uncertainty.mhd"
        )
        if gate_unc_path.exists():
            gate_unc_img = sitk.ReadImage(str(gate_unc_path))
            unc_min, unc_max = roi_dose_uncertainty_bounds(gate_img, gate_unc_img, mask)
            unc_est = roi_dose_uncertainty_estimate(
                gate_img, gate_unc_img, mask, isotope=isotope
            )
        else:
            unc_min, unc_max, unc_est = None, None, None
    else:
        roi_volume_ml = float(bench["roi_volume_ml"])
        roi_mass_g = float(bench["roi_mass_g"])
        gate_mean_dose_rate_gy_s = float(bench["opengate"]["mean_dose_rate_gy_s"])
        op_data = bench.get("opengate", {})
        unc_min = op_data.get("uncertainty_min_fraction")
        unc_max = op_data.get(
            "uncertainty_max_fraction", op_data.get("mean_voxel_uncertainty_fraction")
        )
        unc_est = op_data.get("uncertainty_estimate_fraction")

        s = (roi_volume_ml * 1000.0) ** (1.0 / 3.0)
        activity_img = sitk.Image([1, 1, 1], sitk.sitkFloat32)
        activity_img.SetSpacing((s, s, s))
        activity_img.SetMetaData("unit", "Bq")
        activity_img.SetMetaData("time_from_injection", "0.0")
        activity_img[0, 0, 0] = total_activity_bq

        labels_img = sitk.Image([1, 1, 1], sitk.sitkUInt16)
        labels_img.SetSpacing((s, s, s))
        labels_img[0, 0, 0] = 1

        ct_img = sitk.Image([1, 1, 1], sitk.sitkInt16)
        ct_img.SetSpacing((s, s, s))
        ct_val = int(round((roi_mass_g / roi_volume_ml - 1.0) * 1000.0))
        ct_img[0, 0, 0] = ct_val

    # 1. Madsen 2018 calculation (single Ac-225 parent decay)
    results_madsen = compute_roi_dose(
        image=activity_img,
        labels=labels_img,
        method="madsen2018",
        isotope=isotope,
        ct_image=ct_img,
        labels_map={1: roi_name},
        teff_params={roi_name: teff_h},
        method_params={"phantom": phantom},
    )
    dose_madsen_gy = results_madsen[roi_name]["mean_dose_gy"]

    # 2. LDM calculation (single Ac-225 parent decay)
    results_ldm, _ = compute_roi_dose(
        image=activity_img,
        labels=labels_img,
        method="ldm",
        isotope=isotope,
        ct_image=ct_img,
        labels_map={1: roi_name},
        dose_rate=True,
    )
    ldm_mean_dose_rate_gy_s = results_ldm[roi_name]["mean_dose_rate_gyh"] / 3600.0

    k = np.log(2) / (teff_h * 3600.0)
    ldm_integrated_dose_gy = ldm_mean_dose_rate_gy_s / k
    gate_integrated_dose_gy = gate_mean_dose_rate_gy_s / k

    # Madsen and LDM single-parent calculations match closely (<1%)
    rel_diff_models = (
        abs(dose_madsen_gy - ldm_integrated_dose_gy) / ldm_integrated_dose_gy
    )
    # OpenGate analog simulation tracks the full decay cascade (4 alphas + 2 betas) down to Bi209,
    # 3. LDM calculation with cascade=True (secular equilibrium decay chain: Ac225 -> Bi209)
    results_ldm_cascade, _ = compute_roi_dose(
        image=activity_img,
        labels=labels_img,
        method="ldm",
        isotope=isotope,
        ct_image=ct_img,
        labels_map={1: roi_name},
        dose_rate=True,
        cascade=True,
    )
    ldm_cascade_rate_gy_s = results_ldm_cascade[roi_name]["mean_dose_rate_gyh"] / 3600.0
    ldm_cascade_dose_gy = ldm_cascade_rate_gy_s / k
    rel_diff_cascade = (
        abs(ldm_cascade_dose_gy - gate_integrated_dose_gy) / gate_integrated_dose_gy
    )

    cascade_ratio = gate_integrated_dose_gy / ldm_integrated_dose_gy

    print(f"\nAc225 OpenGate mean dose rate:   {gate_mean_dose_rate_gy_s:.6e} Gy/s")
    print(f"Ac225 OpenGate integrated dose:  {gate_integrated_dose_gy:.6e} Gy")
    print(f"Ac225 Madsen 2018 mean dose:     {dose_madsen_gy:.6e} Gy (parent only)")
    print(
        f"Ac225 LDM integrated dose:       {ldm_integrated_dose_gy:.6e} Gy (parent only)"
    )
    print(
        f"Ac225 LDM cascade dose:          {ldm_cascade_dose_gy:.6e} Gy (diff vs Gate: {rel_diff_cascade * 100:.2f}%)"
    )
    print(
        f"Ac225 Madsen vs LDM diff:        {rel_diff_models * 100:.2f}% (<1% expected)"
    )
    print(f"Ac225 Cascade ratio (Gate/LDM):  {cascade_ratio:.3f}x (~4.72x expected)")
    if unc_est is not None and unc_min is not None and unc_max is not None:
        print(
            f"OpenGate MC uncertainty (liver): {unc_est * 100:.2f}% (bounds: [{unc_min * 100:.3f}%, {unc_max * 100:.2f}%])"
        )

    # Parent-only Madsen and LDM calculations match closely (<1%)
    assert (
        rel_diff_models < 0.01
    ), f"Ac225 Madsen vs LDM difference too large: {rel_diff_models * 100:.2f}% (<1% expected)"
    assert (
        4.5 < cascade_ratio < 5.0
    ), f"Ac225 OpenGate decay cascade ratio unexpected: {cascade_ratio:.3f} (expected ~4.72x)"
    # With cascade=True, LDM matches OpenGate full-chain analog simulation within 5% (actual ~2.68%)
    assert (
        rel_diff_cascade < 0.05
    ), f"Ac225 LDM cascade vs OpenGate difference too large: {rel_diff_cascade * 100:.2f}% (<5% expected)"


def test_madsen2018_cross_dose():
    # pytest -s tests/test_dose.py -k "test_madsen2018_cross_dose" -v
    bench_path = (
        Path(__file__).resolve().parent / "data" / "p142_opengate_benchmark.json"
    )
    if not bench_path.exists():
        pytest.skip("p142 OpenGate benchmark JSON not found.")

    with open(bench_path) as f:
        bench = json.load(f)

    if "cross_dose" not in bench:
        pytest.skip("p142 cross_dose benchmark data not found in JSON.")

    cross_bench = bench["cross_dose"]
    p142_raw_dir = Path(__file__).resolve().parent.parent / "data" / "p142"
    gate_dose_path = p142_raw_dir / "output_vrt" / "output_dose.mhd"
    liver_path = p142_raw_dir / "liver_4mm_cropped.nii.gz"
    ct_path = p142_raw_dir / "ct_4mm_cropped.nii.gz"
    labels_crop_path = p142_raw_dir / "labels_all_4mm_cropped.nii.gz"

    total_activity_bq = float(bench["total_activity_bq"])
    isotope = bench.get("radionuclide", "Lu177")
    teff_h = float(bench["teff_h"])
    time_from_injection_h = float(bench.get("time_from_injection_h", 0.0))
    phantom = bench.get("phantom", "male")
    source_roi = cross_bench.get("source_roi_name", "Liver")
    target_roi = cross_bench["target_roi_name"]

    from radot.db import resolve_organ_name, get_s_value, get_reference_mass

    od_source = resolve_organ_name(source_roi, phantom)
    od_target = resolve_organ_name(target_roi, phantom)
    s_cross = get_s_value(phantom, isotope, od_source, od_target)
    target_ref_mass = get_reference_mass(phantom, od_target)

    if (
        gate_dose_path.exists()
        and liver_path.exists()
        and ct_path.exists()
        and labels_crop_path.exists()
    ):
        liver_img = sitk.ReadImage(str(liver_path))
        mask_source = sitk.GetArrayViewFromImage(liver_img) > 0
        num_source_voxels = int(mask_source.sum())

        act_img = sitk.Cast(liver_img, sitk.sitkFloat32) * (
            total_activity_bq / num_source_voxels
        )
        act_arr = sitk.GetArrayViewFromImage(act_img)

        lbl_crop = sitk.ReadImage(str(labels_crop_path))
        mask_target = sitk.GetArrayViewFromImage(lbl_crop) == 6  # Kidney_R is ID 6

        ct_img = sitk.ReadImage(str(ct_path))
        ct_arr = sitk.GetArrayViewFromImage(ct_img)
        spacing = ct_img.GetSpacing()
        vox_vol = (spacing[0] * spacing[1] * spacing[2]) / 1000.0
        density = np.maximum(0.0, 1.0 + (ct_arr[mask_target] / 1000.0))
        target_mass_g = float(np.sum(density) * vox_vol)

        gate_img = sitk.ReadImage(str(gate_dose_path))
        gate_arr = sitk.GetArrayViewFromImage(gate_img)
        gate_mean_rate_gy_s = float(gate_arr[mask_target].mean())

        gate_unc_path = p142_raw_dir / "output_vrt" / "output_dose_uncertainty.mhd"
        if gate_unc_path.exists():
            gate_unc_img = sitk.ReadImage(str(gate_unc_path))
            unc_min, unc_max = roi_dose_uncertainty_bounds(
                gate_img, gate_unc_img, mask_target
            )
            unc_est = roi_dose_uncertainty_estimate(
                gate_img, gate_unc_img, mask_target, isotope=isotope
            )
        else:
            unc_min, unc_max, unc_est = None, None, None

        mass_scaling = target_mass_g / target_ref_mass
        dose_madsen_cross_gy = dose_madsen2018(
            spect_bq_a=act_arr,
            roi_a=mask_target.astype(np.uint8),
            time_from_injection_h=time_from_injection_h,
            s_value=s_cross,
            mass_scaling=mass_scaling,
            effective_time_h=teff_h,
            source_roi_a=mask_source.astype(np.uint8),
        )
    else:
        target_mass_g = float(cross_bench["target_mass_g"])
        gate_mean_rate_gy_s = float(cross_bench["opengate"]["mean_dose_rate_gy_s"])
        op_data = cross_bench["opengate"]
        unc_min = op_data.get("uncertainty_min_fraction")
        unc_max = op_data.get("uncertainty_max_fraction")
        unc_est = op_data.get("uncertainty_estimate_fraction")

        a0_mbq = total_activity_bq / 1e6
        mass_scaling = target_mass_g / target_ref_mass
        k = np.log(2) / (teff_h * 3600.0)
        integrated_act = a0_mbq * np.exp(k * time_from_injection_h * 3600.0) / k
        dose_madsen_cross_gy = integrated_act * s_cross / mass_scaling / 1000.0

    # Time-integrate OpenGate cross dose rate over source clearance Teff
    k = np.log(2) / (teff_h * 3600.0)
    gate_integrated_cross_dose_gy = gate_mean_rate_gy_s / k

    rel_diff = (
        abs(dose_madsen_cross_gy - gate_integrated_cross_dose_gy)
        / gate_integrated_cross_dose_gy
    )
    print(
        f"\nMadsen 2018 cross dose ({target_roi} <- {source_roi}): {dose_madsen_cross_gy:.6e} Gy"
    )
    print(
        f"OpenGate integrated cross dose ({target_roi}): {gate_integrated_cross_dose_gy:.6e} Gy"
    )
    if unc_est is not None and unc_min is not None and unc_max is not None:
        print(
            f"OpenGate MC uncertainty ({target_roi}): {unc_est * 100:.2f}% (bounds: [{unc_min * 100:.3f}%, {unc_max * 100:.2f}%])"
        )
    print(f"Relative difference: {rel_diff * 100:.2f}%  (<50% expected)")

    # Cross dose is secondary (~100x lower than self dose), driven by photons with greater
    # geometric sensitivity to distance differences between stylized phantom and patient anatomy
    assert (
        rel_diff < 0.50
    ), f"Madsen vs OpenGate cross dose difference too large: {rel_diff * 100:.2f}%"
    assert dose_madsen_cross_gy > 0.0
    assert gate_integrated_cross_dose_gy > 0.0


def test_hanscheid2018_cross_dose():
    # pytest -s tests/test_dose.py -k "test_hanscheid2018_cross_dose" -v
    bench_path = (
        Path(__file__).resolve().parent / "data" / "p142_opengate_benchmark.json"
    )
    if not bench_path.exists():
        pytest.skip("p142 OpenGate benchmark JSON not found.")

    with open(bench_path) as f:
        bench = json.load(f)

    if "cross_dose" not in bench:
        pytest.skip("p142 cross_dose benchmark data not found in JSON.")

    cross_bench = bench["cross_dose"]
    p142_raw_dir = Path(__file__).resolve().parent.parent / "data" / "p142"
    gate_dose_path = p142_raw_dir / "output_vrt" / "output_dose.mhd"
    liver_path = p142_raw_dir / "liver_4mm_cropped.nii.gz"
    ct_path = p142_raw_dir / "ct_4mm_cropped.nii.gz"
    labels_crop_path = p142_raw_dir / "labels_all_4mm_cropped.nii.gz"

    total_activity_bq = float(bench["total_activity_bq"])
    isotope = bench.get("radionuclide", "Lu177")
    teff_h = float(bench["teff_h"])
    # In Hanscheid 2018, Teff = 2 * t, so t = Teff / 2 corresponds to Teff = 45h
    time_from_injection_h = teff_h / 2.0
    phantom = bench.get("phantom", "male")
    source_roi = cross_bench.get("source_roi_name", "Liver")
    target_roi = cross_bench["target_roi_name"]

    from radot.db import resolve_organ_name, get_s_value, get_reference_mass

    od_source = resolve_organ_name(source_roi, phantom)
    od_target = resolve_organ_name(target_roi, phantom)
    s_cross = get_s_value(phantom, isotope, od_source, od_target)
    target_ref_mass = get_reference_mass(phantom, od_target)

    if (
        gate_dose_path.exists()
        and liver_path.exists()
        and ct_path.exists()
        and labels_crop_path.exists()
    ):
        liver_img = sitk.ReadImage(str(liver_path))
        mask_source = sitk.GetArrayViewFromImage(liver_img) > 0
        num_source_voxels = int(mask_source.sum())

        act_img = sitk.Cast(liver_img, sitk.sitkFloat32) * (
            total_activity_bq / num_source_voxels
        )
        act_arr = sitk.GetArrayViewFromImage(act_img)

        lbl_crop = sitk.ReadImage(str(labels_crop_path))
        mask_target = sitk.GetArrayViewFromImage(lbl_crop) == 6  # Kidney_R is ID 6

        ct_img = sitk.ReadImage(str(ct_path))
        ct_arr = sitk.GetArrayViewFromImage(ct_img)
        spacing = ct_img.GetSpacing()
        vox_vol = (spacing[0] * spacing[1] * spacing[2]) / 1000.0
        density = np.maximum(0.0, 1.0 + (ct_arr[mask_target] / 1000.0))
        target_mass_g = float(np.sum(density) * vox_vol)

        gate_img = sitk.ReadImage(str(gate_dose_path))
        gate_arr = sitk.GetArrayViewFromImage(gate_img)
        gate_mean_rate_gy_s = float(gate_arr[mask_target].mean())

        gate_unc_path = p142_raw_dir / "output_vrt" / "output_dose_uncertainty.mhd"
        if gate_unc_path.exists():
            gate_unc_img = sitk.ReadImage(str(gate_unc_path))
            unc_min, unc_max = roi_dose_uncertainty_bounds(
                gate_img, gate_unc_img, mask_target
            )
            unc_est = roi_dose_uncertainty_estimate(
                gate_img, gate_unc_img, mask_target, isotope=isotope
            )
        else:
            unc_min, unc_max, unc_est = None, None, None

        mass_scaling = target_mass_g / target_ref_mass
        dose_h2018_cross_gy = dose_hanscheid2018(
            spect_bq_a=act_arr,
            roi_a=mask_target.astype(np.uint8),
            time_from_injection_h=time_from_injection_h,
            s_value=s_cross,
            mass_scaling=mass_scaling,
            source_roi_a=mask_source.astype(np.uint8),
        )
    else:
        target_mass_g = float(cross_bench["target_mass_g"])
        gate_mean_rate_gy_s = float(cross_bench["opengate"]["mean_dose_rate_gy_s"])
        op_data = cross_bench["opengate"]
        unc_min = op_data.get("uncertainty_min_fraction")
        unc_max = op_data.get("uncertainty_max_fraction")
        unc_est = op_data.get("uncertainty_estimate_fraction")

        a0_mbq = total_activity_bq / 1e6
        mass_scaling = target_mass_g / target_ref_mass
        dose_h2018_cross_gy = (
            a0_mbq
            * (2 * time_from_injection_h * 3600.0)
            / np.log(2)
            * s_cross
            / mass_scaling
            / 1000.0
        )

    # Time-integrate OpenGate cross dose rate over source clearance Teff
    k = np.log(2) / (teff_h * 3600.0)
    gate_integrated_cross_dose_gy = gate_mean_rate_gy_s / k

    rel_diff = (
        abs(dose_h2018_cross_gy - gate_integrated_cross_dose_gy)
        / gate_integrated_cross_dose_gy
    )
    print(
        f"\nH\u00e4nscheid 2018 cross dose ({target_roi} <- {source_roi}): {dose_h2018_cross_gy:.6e} Gy"
    )
    print(
        f"OpenGate integrated cross dose ({target_roi}): {gate_integrated_cross_dose_gy:.6e} Gy"
    )
    if unc_est is not None and unc_min is not None and unc_max is not None:
        print(
            f"OpenGate MC uncertainty ({target_roi}): {unc_est * 100:.2f}% (bounds: [{unc_min * 100:.3f}%, {unc_max * 100:.2f}%])"
        )
    print(f"Relative difference: {rel_diff * 100:.2f}%  (<50% expected)")

    assert (
        rel_diff < 0.50
    ), f"Hanscheid 2018 vs OpenGate cross dose difference too large: {rel_diff * 100:.2f}%"
    assert dose_h2018_cross_gy > 0.0
    assert gate_integrated_cross_dose_gy > 0.0


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
    # 4 voxels with large spacing 10x10x10 mm, isotope Lu177 (range = 1.8 mm)
    # nu_cluster = (1 + 1.8/10)^3 = 1.18^3 = 1.643
    dose = np.ones((2, 2, 1), dtype=np.float32)
    unc = np.ones((2, 2, 1), dtype=np.float32) * 0.05
    mask = np.ones((2, 2, 1), dtype=np.uint8)

    unc_est = roi_dose_uncertainty_estimate(
        dose, unc, mask, spacing_mm=(10.0, 10.0, 10.0), isotope="Lu177"
    )
    # unc_min = 0.025, unc_max = 0.05
    # expected: 0.025 * sqrt(1.643032) approx 0.03204
    assert unc_est == pytest.approx(0.025 * np.sqrt(1.18**3), rel=1e-3)
    assert unc_est >= 0.025
    assert unc_est <= 0.05
