import pytest
import numpy as np
import SimpleITK as sitk
from click.testing import CliRunner
from radot.convert import convert_image_units


def test_convert_bq_to_bqml():
    # Voxel volume: 2mm x 2mm x 2mm = 8 mm^3 = 0.008 mL
    image = sitk.Image(5, 5, 5, sitk.sitkFloat32)
    image.SetSpacing((2.0, 2.0, 2.0))

    # Fill with 8 Bq per voxel
    arr = np.full((5, 5, 5), 8.0, dtype=np.float32)
    image = sitk.GetImageFromArray(arr)
    image.SetSpacing((2.0, 2.0, 2.0))

    # Convert to bqml: should be 8 Bq / 0.008 mL = 1000 Bq/mL
    converted = convert_image_units(image, from_unit="bq", to_unit="bqml")
    res_arr = sitk.GetArrayFromImage(converted)

    assert res_arr[0, 0, 0] == pytest.approx(1000.0)

    # Convert back to bq: should be 8.0 Bq
    back = convert_image_units(converted, from_unit="bqml", to_unit="bq")
    back_arr = sitk.GetArrayFromImage(back)
    assert back_arr[0, 0, 0] == pytest.approx(8.0)


def test_convert_bqml_to_suv():
    # Voxel volume is not relevant for bqml <-> suv
    image = sitk.Image(5, 5, 5, sitk.sitkFloat32)
    arr = np.full((5, 5, 5), 1000.0, dtype=np.float32)  # 1000 Bq/mL
    image = sitk.GetImageFromArray(arr)
    image.SetSpacing((1.0, 1.0, 1.0))

    # Weight = 75 kg = 75000 g
    # Activity = 250 MBq = 250,000,000 Bq
    # Expected SUV factor = 250,000,000 / 75,000 = 3333.333 Bq/g
    # SUV = 1000 / 3333.333 = 0.3
    converted = convert_image_units(
        image, from_unit="bqml", to_unit="suv", weight=75.0, activity=250.0, decay_factor=1.0
    )
    res_arr = sitk.GetArrayFromImage(converted)
    assert res_arr[0, 0, 0] == pytest.approx(0.3)

    # Convert back
    back = convert_image_units(
        converted, from_unit="suv", to_unit="bqml", weight=75.0, activity=250.0, decay_factor=1.0
    )
    back_arr = sitk.GetArrayFromImage(back)
    assert back_arr[0, 0, 0] == pytest.approx(1000.0)


def test_convert_validation():
    image = sitk.Image(5, 5, 5, sitk.sitkFloat32)

    # Missing weight/activity for SUV conversion
    with pytest.raises(ValueError, match="weight.*required"):
        convert_image_units(
            image, from_unit="bqml", to_unit="suv", weight=None, activity=250.0
        )

    with pytest.raises(ValueError, match="activity.*required"):
        convert_image_units(
            image, from_unit="bqml", to_unit="suv", weight=75.0, activity=None
        )


def test_convert_integer_input():
    # Input is an Int16 image
    image = sitk.Image(5, 5, 5, sitk.sitkInt16)
    arr = np.full((5, 5, 5), 10, dtype=np.int16)
    image = sitk.GetImageFromArray(arr)
    image.SetSpacing((2.0, 2.0, 2.0))

    # Convert to bqml: should cast to Float32 and execute correctly
    converted = convert_image_units(image, from_unit="bq", to_unit="bqml")

    # Assert output data type is Float32
    assert converted.GetPixelID() == sitk.sitkFloat32

    res_arr = sitk.GetArrayFromImage(converted)
    # 10 Bq / 0.008 mL = 1250.0 Bq/mL
    assert res_arr[0, 0, 0] == pytest.approx(1250.0)


def test_convert_metadata_fallback(tmp_path):
    from click.testing import CliRunner
    from radot.cli import main
    from radot.meta import manage_metadata

    image_path = tmp_path / "pet_convert.nii.gz"
    image = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    # Fill with 1000.0 Bq/mL
    arr = np.full((5, 5, 5), 1000.0, dtype=np.float32)
    image = sitk.GetImageFromArray(arr)
    image.SetSpacing((1.0, 1.0, 1.0))
    sitk.WriteImage(image, str(image_path))

    # Set metadata: unit=bqml, modality=pet, weight=75.0, activity=250.0, scan at the injection time
    manage_metadata(
        image_path,
        modality="pet",
        unit="bqml",
        weight=75.0,
        activity=250.0,
        injection_date="2026-01-10 10:00:00",
        acquisition_date="2026-01-10 10:00:00",
    )
    _set_radionuclide(image_path, "F18")

    runner = CliRunner()
    # Run conversion without specifying --from-unit, --weight, or --activity (should read from metadata)
    output_path = tmp_path / "out_suv.nii.gz"
    result = runner.invoke(
        main, ["convert-units", str(image_path), "-t", "suv", "-o", str(output_path)]
    )
    assert result.exit_code == 0
    assert output_path.exists()

    res_img = sitk.ReadImage(str(output_path))
    res_arr = sitk.GetArrayFromImage(res_img)
    # Expected SUV = 1000 / (250,000,000 / 75,000) = 0.3
    assert res_arr[0, 0, 0] == pytest.approx(0.3)


def test_convert_missing_metadata_errors(tmp_path):
    from click.testing import CliRunner
    from radot.cli import main

    image_path = tmp_path / "pet_no_meta.nii.gz"
    image = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    sitk.WriteImage(image, str(image_path))

    runner = CliRunner()
    # Fails because source unit is not specified and not in metadata
    result = runner.invoke(main, ["convert-units", str(image_path), "-t", "suv"])
    assert result.exit_code != 0
    assert "Source unit is not specified" in result.output


def test_calculate_lbm():
    from radot.convert import calculate_lbm

    # James 1976 formula
    # Male: weight = 70.0, height = 175.0
    # LBM = 1.10 * 70 - 128 * (70 / 175)^2 = 77.0 - 128 * 0.16 = 77.0 - 20.48 = 56.52
    lbm_m = calculate_lbm(weight=70.0, height=175.0, sex="male", method="james")
    assert lbm_m == pytest.approx(56.52)

    # Female: weight = 60.0, height = 160.0
    # LBM = 1.07 * 60 - 148 * (60 / 160)^2 = 64.2 - 148 * 0.140625 = 64.2 - 20.8125 = 43.3875
    lbm_f = calculate_lbm(weight=60.0, height=160.0, sex="female", method="james")
    assert lbm_f == pytest.approx(43.3875)

    # Janmahasatian 2005 formula
    # Male: weight = 70.0, height = 175.0
    # BMI = 70 / (1.75)^2 = 22.857
    # LBM = (9270 * 70) / (6680 + 216 * 22.857) = 648900 / (6680 + 4937.11) = 648900 / 11617.11 = 55.857
    lbm_jan_m = calculate_lbm(
        weight=70.0, height=175.0, sex="male", method="janmahasatian"
    )
    assert lbm_jan_m == pytest.approx(55.857, rel=1e-3)


def test_calculate_lbm_segmentation():
    from radot.convert import calculate_lbm

    # Generate a dummy label map image
    # Size 5x5x5 = 125 voxels. Spacing 2.0 mm, voxel volume = 0.008 mL
    labels = sitk.Image(5, 5, 5, sitk.sitkUInt8)
    labels.SetSpacing((2.0, 2.0, 2.0))
    # Let's say label ID 1 is subcutaneous_fat, fill it in half the image (60 voxels)
    arr = np.zeros((5, 5, 5), dtype=np.uint8)
    arr[:3, :, :] = 1  # 3 * 5 * 5 = 75 voxels
    labels = sitk.GetImageFromArray(arr)
    labels.SetSpacing((2.0, 2.0, 2.0))

    labels_map = {1: "VAT_tissue", 2: "liver"}
    # VAT stands for Visceral Adipose Tissue, should be matched as fat
    # Fat voxels = 75. Fat volume = 75 * 0.008 mL = 0.6 mL
    # Fat mass = 0.6 * 0.92 g/mL * 1e-3 kg/g = 0.000552 kg
    # LBM = weight - fat_mass = 70.0 - 0.000552 = 69.999448 kg
    lbm = calculate_lbm(
        weight=70.0, method="segmentation", labels=labels, labels_map=labels_map
    )
    assert lbm == pytest.approx(69.999448)


def test_convert_bqml_to_sul_james():
    # Voxel volume is not relevant for bqml <-> sul
    image = sitk.Image(5, 5, 5, sitk.sitkFloat32)
    arr = np.full((5, 5, 5), 1000.0, dtype=np.float32)  # 1000 Bq/mL
    image = sitk.GetImageFromArray(arr)
    image.SetSpacing((1.0, 1.0, 1.0))

    # Weight = 70.0 kg, height = 175.0 cm, sex = male => LBM = 56.52 kg
    # Expected SUL = bqml / (activity_bq / lbm_g)
    # activity_bq = 250 MBq = 250,000,000 Bq
    # lbm_g = 56.52 * 1000 = 56,520 g
    # Expected SUL = 1000 / (250,000,000 / 56,520) = 1000 / 4423.213 = 0.22608
    converted = convert_image_units(
        image,
        from_unit="bqml",
        to_unit="sul",
        weight=70.0,
        activity=250.0,
        height=175.0,
        sex="male",
        lbm_method="james",
        decay_factor=1.0,
    )
    res_arr = sitk.GetArrayFromImage(converted)
    assert res_arr[0, 0, 0] == pytest.approx(0.22608, rel=1e-3)


def test_convert_units_cli_sul(tmp_path):
    from click.testing import CliRunner
    from radot.cli import main
    from radot.meta import manage_metadata

    image_path = tmp_path / "pet_convert_sul.nii.gz"
    image = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    arr = np.full((5, 5, 5), 1000.0, dtype=np.float32)
    image = sitk.GetImageFromArray(arr)
    image.SetSpacing((1.0, 1.0, 1.0))
    sitk.WriteImage(image, str(image_path))

    # Write metadata first
    manage_metadata(
        image_path,
        modality="pet",
        unit="bqml",
        weight=70.0,
        activity=250.0,
        height=175.0,
        sex="male",
        injection_date="2026-01-10 10:00:00",
        acquisition_date="2026-01-10 10:00:00",
    )
    _set_radionuclide(image_path, "F18")

    runner = CliRunner()
    output_path = tmp_path / "out_sul.nii.gz"
    result = runner.invoke(
        main,
        [
            "convert-units",
            str(image_path),
            "-t",
            "sul",
            "-o",
            str(output_path),
            "--lbm-method",
            "james",
        ],
    )
    assert result.exit_code == 0
    assert output_path.exists()

    res_img = sitk.ReadImage(str(output_path))
    res_arr = sitk.GetArrayFromImage(res_img)
    # LBM = 56.52 kg. SUL = 1000 / (250,000,000 / 56,520) = 0.22608
    assert res_arr[0, 0, 0] == pytest.approx(0.22608, rel=1e-3)


def test_compute_decay_correction():
    from radot.convert import compute_decay_correction

    # Ga-68 (half-life = 67.71 min) with exactly 67.71 min delay => factor = 2.0
    factor, delay = compute_decay_correction(
        injection_date="2021-04-28 10:00:00",
        acquisition_date="2021-04-28 11:07:42.600000",
        radionuclide="Ga-68",
    )
    assert delay == pytest.approx(67.71, rel=1e-3)
    assert factor == pytest.approx(2.0, rel=1e-3)

    # F-18 (half-life = 109.77 min), other date format
    factor_f18, delay_f18 = compute_decay_correction(
        injection_date="2021/04/28, 10:00:00",
        acquisition_date="2021/04/28, 11:49:46.2",
        radionuclide="F18",
    )
    assert delay_f18 == pytest.approx(109.77, rel=1e-3)
    assert factor_f18 == pytest.approx(2.0, rel=1e-3)

    # A tracer name is not a radionuclide ('psma' can be Ga-68 or Lu-177)
    for tracer in ("psma", "fdg"):
        with pytest.raises(ValueError):
            compute_decay_correction("2021-04-28 10:00:00", "2021-04-28 11:00:00", tracer)
    # acquisition before injection, unreadable date
    with pytest.raises(ValueError):
        compute_decay_correction("2021-04-28 10:00:00", "2021-04-28 09:00:00", "F18")
    with pytest.raises(ValueError):
        compute_decay_correction("28 April", "2021-04-28 11:00:00", "F18")

    # Explicit decay factor
    factor_exp, _ = compute_decay_correction(decay_factor=1.85)
    assert factor_exp == 1.85

    # Missing inputs => returns (1.0, None)
    factor_none, delay_none = compute_decay_correction()
    assert factor_none == 1.0
    assert delay_none is None


def test_convert_with_decay_correction_api():
    # 1000 Bq/mL, Weight 75 kg, Activity 250 MBq
    # Without decay: factor = 1.0 => SUV = 0.3
    # With decay: delay = 67.71 min, factor = 2.0 => effective_activity = 125 MBq => SUV = 0.6
    image = sitk.Image(5, 5, 5, sitk.sitkFloat32)
    arr = np.full((5, 5, 5), 1000.0, dtype=np.float32)
    image = sitk.GetImageFromArray(arr)
    image.SetSpacing((1.0, 1.0, 1.0))

    converted = convert_image_units(
        image,
        from_unit="bqml",
        to_unit="suv",
        weight=75.0,
        activity=250.0,
        radionuclide="Ga-68",
        injection_date="2021-04-28 10:00:00",
        acquisition_date="2021-04-28 11:07:42.6",
    )
    res_arr = sitk.GetArrayFromImage(converted)
    assert res_arr[0, 0, 0] == pytest.approx(0.6, rel=1e-3)

    # Invert back to bqml: should yield 1000.0 Bq/mL
    back = convert_image_units(
        converted,
        from_unit="suv",
        to_unit="bqml",
        weight=75.0,
        activity=250.0,
        radionuclide="Ga-68",
        injection_date="2021-04-28 10:00:00",
        acquisition_date="2021-04-28 11:07:42.6",
    )
    back_arr = sitk.GetArrayFromImage(back)
    assert back_arr[0, 0, 0] == pytest.approx(1000.0, rel=1e-3)


def test_convert_units_cli_with_decay(tmp_path):
    import json
    from click.testing import CliRunner
    from radot.cli import main

    image_path = tmp_path / "pet_decay.nii.gz"
    image = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    arr = np.full((5, 5, 5), 1000.0, dtype=np.float32)
    image = sitk.GetImageFromArray(arr)
    image.SetSpacing((1.0, 1.0, 1.0))
    sitk.WriteImage(image, str(image_path))

    output_path = tmp_path / "out_decay_suv.nii.gz"
    runner = CliRunner()
    result = runner.invoke(
        main,
        [
            "convert-units",
            str(image_path),
            "-f",
            "bqml",
            "-t",
            "suv",
            "-o",
            str(output_path),
            "--weight",
            "75.0",
            "--activity",
            "250.0",
            "--radionuclide",
            "Ga-68",
            "--injection-date",
            "2021-04-28 10:00:00",
            "--acquisition-date",
            "2021-04-28 11:07:42.6",
        ],
    )
    assert result.exit_code == 0
    assert output_path.exists()

    sidecar_path = tmp_path / "out_decay_suv.json"
    assert sidecar_path.exists()
    with open(sidecar_path, "r") as f:
        sidecar = json.load(f)

    assert sidecar["unit"] == "suv"
    assert sidecar["activity"] == pytest.approx(250.0)
    assert sidecar["decay_correction_factor"] == pytest.approx(2.0, rel=1e-3)
    assert sidecar["decay_delay_minutes"] == pytest.approx(67.71, rel=1e-3)
    assert sidecar["decay_correction_reference"] == "injection_time"
    assert "sha256" in sidecar


def _pet(tmp_path, **meta):
    from radot.meta import manage_metadata

    image_path = tmp_path / "pet.nii.gz"
    image = sitk.GetImageFromArray(np.full((5, 5, 5), 1000.0, dtype=np.float32))
    sitk.WriteImage(image, str(image_path))
    rad = meta.pop("radionuclide", None)
    manage_metadata(image_path, modality="pet", unit="bqml", weight=75.0, activity=250.0, **meta)
    if rad:
        _set_radionuclide(image_path, rad)
    return image_path


def test_convert_suv_decay_required(tmp_path):
    from radot.cli import main

    # no radionuclide nor dates: explicit error, not a silent conversion without decay
    image_path = _pet(tmp_path)
    with pytest.raises(ValueError, match="radionuclide, injection_date, acquisition_date"):
        convert_image_units(sitk.ReadImage(str(image_path)), to_unit="suv", image_path=image_path)
    result = CliRunner().invoke(main, ["convert-units", str(image_path), "-t", "suv"])
    assert result.exit_code != 0 and "decay-factor" in result.output
    # explicitly without decay
    out = tmp_path / "out.nii.gz"
    result = CliRunner().invoke(main, ["convert-units", str(image_path), "-t", "suv", "--decay-factor", "1", "-o", str(out)])
    assert result.exit_code == 0, result.output
    assert sitk.GetArrayFromImage(sitk.ReadImage(str(out)))[0, 0, 0] == pytest.approx(0.3)


def test_convert_suv_decay_from_sidecar(tmp_path):
    import json
    from radot.cli import main

    # F-18, one half-life after injection: SUV doubled; decay keys in the output sidecar
    image_path = _pet(tmp_path, radionuclide="F18", injection_date="2021-04-28 10:00:00", acquisition_date="2021-04-28 11:49:46")
    out = tmp_path / "out_suv.nii.gz"
    result = CliRunner().invoke(main, ["convert-units", str(image_path), "-t", "suv", "-o", str(out)])
    assert result.exit_code == 0, result.output
    assert sitk.GetArrayFromImage(sitk.ReadImage(str(out)))[0, 0, 0] == pytest.approx(0.6, rel=1e-3)
    meta = json.loads((tmp_path / "out_suv.json").read_text())
    assert meta["decay_correction_factor"] == pytest.approx(2.0, rel=1e-3)
    assert meta["decay_correction_reference"] == "injection_time"
    # back to Bq/mL: activity at the acquisition time again, no decay keys
    back = tmp_path / "back.nii.gz"
    result = CliRunner().invoke(main, ["convert-units", str(out), "-t", "bqml", "-o", str(back)])
    assert result.exit_code == 0, result.output
    assert sitk.GetArrayFromImage(sitk.ReadImage(str(back)))[0, 0, 0] == pytest.approx(1000.0, rel=1e-3)
    assert "decay_correction_reference" not in json.loads((tmp_path / "back.json").read_text())


def test_convert_suv_decay_admin(tmp_path):
    import json

    # image decay-corrected to the injection time (DICOM ADMIN): the injected activity is not decayed
    image_path = _pet(tmp_path, radionuclide="F18", injection_date="2021-04-28 10:00:00", acquisition_date="2021-04-28 11:49:46")
    sidecar = tmp_path / "pet.json"
    meta = json.loads(sidecar.read_text())
    meta["decay_correction"] = "admin"
    sidecar.write_text(json.dumps(meta))
    converted = convert_image_units(sitk.ReadImage(str(image_path)), to_unit="suv", image_path=image_path)
    assert sitk.GetArrayFromImage(converted)[0, 0, 0] == pytest.approx(0.3, rel=1e-3)


def _set_radionuclide(image_path, radionuclide):
    """Radionuclide in the sidecar (radot meta has no option for it on this branch)."""
    import json
    from radot.meta import get_sidecar_path

    sidecar = get_sidecar_path(image_path)
    meta = json.loads(sidecar.read_text())
    meta["radionuclide"] = radionuclide
    sidecar.write_text(json.dumps(meta))


def test_convert_suv_to_sul_no_dates(tmp_path):
    # SUV -> SUL: the decayed activity cancels out, no dates needed
    image = sitk.GetImageFromArray(np.full((5, 5, 5), 2.0, dtype=np.float32))
    sul = convert_image_units(image, from_unit="suv", to_unit="sul", weight=70.0, activity=250.0, height=175.0, sex="male")
    # SUL = SUV x LBM / W, LBM (James, male 70 kg 175 cm) = 56.52 kg
    assert sitk.GetArrayFromImage(sul)[0, 0, 0] == pytest.approx(2.0 * 56.52 / 70.0, rel=1e-3)
