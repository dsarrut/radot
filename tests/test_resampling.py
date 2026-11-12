import pytest
from pathlib import Path
import SimpleITK as sitk
import numpy as np
from radot.resampling import resample_image, resample_image_file

def test_resample_image_ct():
    # Create 2x2x2 image with 1mm spacing
    ct_data = np.zeros((10, 10, 10), dtype=np.int16)
    ct_data[3:7, 3:7, 3:7] = 100
    image = sitk.GetImageFromArray(ct_data)
    image.SetSpacing((1.0, 1.0, 1.0))

    # Resample to 2mm spacing
    resampled = resample_image(image, (2.0, 2.0, 2.0), image_type="ct", gaussian=False)
    
    assert resampled.GetSpacing() == (2.0, 2.0, 2.0)
    assert resampled.GetSize() == (5, 5, 5)

def test_resample_image_labels():
    # Labels must use NearestNeighbor interpolation
    label_data = np.zeros((10, 10, 10), dtype=np.uint8)
    label_data[4:6, 4:6, 4:6] = 5
    image = sitk.GetImageFromArray(label_data)
    image.SetSpacing((1.0, 1.0, 1.0))

    # Resample to 2.5mm spacing
    resampled = resample_image(image, (2.5, 2.5, 2.5), image_type="labels", gaussian=False)
    
    res_arr = sitk.GetArrayFromImage(resampled)
    # Check that label values are kept as integers (0 and 5 only, no interpolated float values)
    unique_vals = np.unique(res_arr)
    for v in unique_vals:
        assert v in [0, 5]

def test_resample_image_preserve_total():
    # Create absolute activity image (Bq per voxel)
    activity_data = np.zeros((10, 10, 10), dtype=np.float32)
    activity_data[3:7, 3:7, 3:7] = 10.0
    image = sitk.GetImageFromArray(activity_data)
    image.SetSpacing((1.0, 1.0, 1.0))
    
    original_sum = np.sum(activity_data)  # 4 * 4 * 4 * 10 = 640

    # Resample to 2.0mm spacing with preserve_total=True
    resampled = resample_image(image, (2.0, 2.0, 2.0), image_type="pet", preserve_total=True, gaussian=False)
    res_arr = sitk.GetArrayFromImage(resampled)
    resampled_sum = np.sum(res_arr)
    
    # The sum of voxel values should be preserved (allowing small interpolation edge differences)
    assert pytest.approx(resampled_sum, rel=1e-2) == original_sum

def test_resample_image_gaussian():
    # Verify that applying Gaussian filter works and reduces the maximum gradient
    ct_data = np.zeros((15, 15, 15), dtype=np.float32)
    ct_data[5:10, 5:10, 5:10] = 500.0
    image = sitk.GetImageFromArray(ct_data)
    image.SetSpacing((1.0, 1.0, 1.0))

    # Resample with gaussian=True
    resampled_smooth = resample_image(image, (3.0, 3.0, 3.0), image_type="ct", gaussian=True)
    resampled_blocky = resample_image(image, (3.0, 3.0, 3.0), image_type="ct", gaussian=False)
    
    arr_smooth = sitk.GetArrayFromImage(resampled_smooth)
    arr_blocky = sitk.GetArrayFromImage(resampled_blocky)
    
    # Smooth resampled image should have lower maximum value due to blurring
    assert np.max(arr_smooth) < np.max(arr_blocky)

def test_resample_image_file(tmp_path):
    # Setup dummy CT file
    ct_data = np.zeros((6, 6, 6), dtype=np.int16)
    image = sitk.GetImageFromArray(ct_data)
    image.SetSpacing((1.0, 1.0, 1.0))
    
    input_path = tmp_path / "input_ct.nii.gz"
    sitk.WriteImage(image, str(input_path))
    
    output_path = tmp_path / "output_ct_3mm.nii.gz"
    
    # Run the file resampling
    resample_image_file(
        image_path=input_path,
        spacing=(3.0, 3.0, 3.0),
        image_type="ct",
        output_path=output_path,
        unit="hu",
        preserve_total=False,
        gaussian=False
    )
    
    assert output_path.exists()
    out_img = sitk.ReadImage(str(output_path))
    assert out_img.GetSpacing() == (3.0, 3.0, 3.0)
    assert out_img.GetSize() == (2, 2, 2)


def test_resample_metadata_fallback(tmp_path):
    from click.testing import CliRunner
    from radot.cli import main
    from radot.meta import manage_metadata

    image_path = tmp_path / "pet_resample.nii.gz"
    image = sitk.Image([6, 6, 6], sitk.sitkFloat32)
    image.SetSpacing((1.0, 1.0, 1.0))
    sitk.WriteImage(image, str(image_path))

    # Set metadata: modality=pet, unit=bq
    manage_metadata(image_path, modality="pet", unit="bq")

    output_path = tmp_path / "output_resample.nii.gz"

    runner = CliRunner()
    # Run resample CLI without --type or --unit (should resolve modality and unit from metadata)
    # Since unit is 'bq', it should preserve total activity
    result = runner.invoke(main, ["resample", str(image_path), "-s", "2", "2", "2", "-o", str(output_path)])
    assert result.exit_code == 0
    assert output_path.exists()

    res_img = sitk.ReadImage(str(output_path))
    assert res_img.GetSpacing() == (2.0, 2.0, 2.0)


def test_resample_missing_metadata_errors(tmp_path):
    from click.testing import CliRunner
    from radot.cli import main

    image_path = tmp_path / "no_meta_resample.nii.gz"
    image = sitk.Image([6, 6, 6], sitk.sitkFloat32)
    image.SetSpacing((1.0, 1.0, 1.0))
    sitk.WriteImage(image, str(image_path))

    output_path = tmp_path / "output_resample_fail.nii.gz"

    runner = CliRunner()
    # Fails because type/modality is not specified and not in metadata
    result = runner.invoke(main, ["resample", str(image_path), "-s", "2", "2", "2", "-o", str(output_path)])
    assert result.exit_code != 0
    assert "Image modality type is not specified" in result.output

    # Fails because modality is PET but no unit is provided / in metadata
    from radot.meta import manage_metadata
    pet_image_path = tmp_path / "pet_no_unit.nii.gz"
    sitk.WriteImage(image, str(pet_image_path))
    manage_metadata(pet_image_path, modality="pet")
    
    result2 = runner.invoke(main, ["resample", str(pet_image_path), "-s", "2", "2", "2", "-o", str(output_path)])
    assert result2.exit_code != 0
    assert "unit option is required" in result2.output


def test_resample_image_file_like(tmp_path):
    from radot.resampling import resample_image_file_like
    from radot.meta import manage_metadata

    # Create source image (6x6x6, 1mm spacing)
    src_path = tmp_path / "src.nii.gz"
    src_img = sitk.Image([6, 6, 6], sitk.sitkFloat32)
    src_img.SetSpacing((1.0, 1.0, 1.0))
    sitk.WriteImage(src_img, str(src_path))

    # Create like image (3x3x3, 2mm spacing)
    like_path = tmp_path / "like.nii.gz"
    like_img = sitk.Image([3, 3, 3], sitk.sitkUInt16)
    like_img.SetSpacing((2.0, 2.0, 2.0))
    sitk.WriteImage(like_img, str(like_path))

    out_path = tmp_path / "out.nii.gz"

    # 1. Fail because no type parameter and no modality in metadata
    with pytest.raises(ValueError, match="Image modality type is not specified"):
        resample_image_file_like(src_path, like_path, out_path)

    # 2. Succeed when providing explicit type
    resample_image_file_like(src_path, like_path, out_path, type="ct")
    assert out_path.exists()
    res_img = sitk.ReadImage(str(out_path))
    assert res_img.GetSpacing() == (2.0, 2.0, 2.0)
    assert res_img.GetSize() == (3, 3, 3)

    # Clean up output
    out_path.unlink()

    # 3. Set metadata to "spect", "bq"
    manage_metadata(src_path, modality="spect", unit="bq")

    # 4. Succeed if type parameter is ignored/omitted since type is in metadata
    resample_image_file_like(src_path, like_path, out_path)
    assert out_path.exists()
    out_img = sitk.ReadImage(str(out_path))
    assert out_img.GetSpacing() == (2.0, 2.0, 2.0)

    # Clean up output
    out_path.unlink()

    # 5. Succeed if both provided and they match
    resample_image_file_like(src_path, like_path, out_path, type="spect")
    assert out_path.exists()

    # Clean up output
    out_path.unlink()

    # 6. Fail if both provided and they are different (e.g. meta="spect", parameter="ct")
    with pytest.raises(ValueError, match="Image type mismatch"):
        resample_image_file_like(src_path, like_path, out_path, type="ct")


def test_resample_cli_like(tmp_path):
    from click.testing import CliRunner
    from radot.cli import main
    from radot.meta import manage_metadata

    src_path = tmp_path / "src.nii.gz"
    src_img = sitk.Image([6, 6, 6], sitk.sitkFloat32)
    src_img.SetSpacing((1.0, 1.0, 1.0))
    sitk.WriteImage(src_img, str(src_path))
    manage_metadata(src_path, modality="spect", unit="bq")

    like_path = tmp_path / "like.nii.gz"
    like_img = sitk.Image([3, 3, 3], sitk.sitkUInt16)
    like_img.SetSpacing((2.0, 2.0, 2.0))
    sitk.WriteImage(like_img, str(like_path))

    out_path = tmp_path / "out.nii.gz"

    runner = CliRunner()

    # 1. Test resampling with --like option
    result = runner.invoke(main, [
        "resample",
        str(src_path),
        "--like", str(like_path),
        "-o", str(out_path),
    ])
    assert result.exit_code == 0
    assert out_path.exists()

    res_img = sitk.ReadImage(str(out_path))
    assert res_img.GetSpacing() == (2.0, 2.0, 2.0)
    assert res_img.GetSize() == (3, 3, 3)

    # 2. Test failure if neither --spacing nor --like is specified
    result_none = runner.invoke(main, [
        "resample",
        str(src_path),
        "-o", str(out_path),
    ])
    assert result_none.exit_code != 0
    assert "Either --spacing (-s) or --like (-l) must be specified." in result_none.output

    # 3. Test failure if both --spacing and --like are specified
    result_both = runner.invoke(main, [
        "resample",
        str(src_path),
        "-s", "2", "2", "2",
        "--like", str(like_path),
        "-o", str(out_path),
    ])
    assert result_both.exit_code != 0
    assert "Cannot specify both --spacing (-s) and --like (-l)." in result_both.output





def _activity_file(tmp_path, name, sidecar):
    """4 mm image, uniform 10 per voxel, with a sidecar."""
    import json
    from radot.meta import get_sidecar_path

    image = sitk.GetImageFromArray(np.full((8, 8, 8), 10.0, dtype=np.float32))
    image.SetSpacing((4.0, 4.0, 4.0))
    path = tmp_path / name
    sitk.WriteImage(image, str(path))
    get_sidecar_path(path).write_text(json.dumps(sidecar))
    return path


@pytest.mark.parametrize("unit", ["bq", "bqs", "BQS"])
def test_resample_file_per_voxel_unit_keeps_total(tmp_path, unit):
    """A per-voxel unit (activity bq, time-integrated activity bqs) keeps its total, also when the
    sidecar unit is not lower case."""
    path = _activity_file(tmp_path, "a.nii.gz", {"modality": "spect", "unit": unit})
    out = tmp_path / "a_2mm.nii.gz"
    resample_image_file(path, (2.0, 2.0, 2.0), output_path=out, gaussian=False)
    total = sitk.GetArrayFromImage(sitk.ReadImage(str(out))).sum()
    assert total == pytest.approx(8**3 * 10.0, rel=1e-3)


def test_resample_like_per_voxel_unit_keeps_total(tmp_path):
    from radot.resampling import resample_image_file_like

    path = _activity_file(tmp_path, "tia.nii.gz", {"modality": "spect", "unit": "bqs"})
    like = sitk.Image([16, 16, 16], sitk.sitkFloat32)
    like.SetSpacing((2.0, 2.0, 2.0))
    like.SetOrigin((-1.0, -1.0, -1.0))  # same extent as the 4 mm image
    like_path = tmp_path / "like.nii.gz"
    sitk.WriteImage(like, str(like_path))
    out = tmp_path / "tia_like.nii.gz"
    resample_image_file_like(path, like_path, out, gaussian=False)
    total = sitk.GetArrayFromImage(sitk.ReadImage(str(out))).sum()
    assert total == pytest.approx(8**3 * 10.0, rel=1e-3)


def test_resample_output_sidecar(tmp_path):
    """The output sidecar keeps the source keys, writes the resolved modality and unit, and the SHA-256
    of the new file (not the source one)."""
    import json
    from radot.meta import compute_file_sha256, get_sidecar_path

    src = {"unit": "bqml", "radionuclide": "Lu177", "weight": 70.0, "sha256": "0" * 64}
    path = _activity_file(tmp_path, "s.nii.gz", src)
    out = tmp_path / "s_2mm.nii.gz"
    resample_image_file(path, (2.0, 2.0, 2.0), image_type="spect", output_path=out)
    meta = json.loads(get_sidecar_path(out).read_text())
    assert meta["sha256"] == compute_file_sha256(out)
    assert meta["modality"] == "spect"
    assert meta["unit"] == "bqml"
    assert meta["radionuclide"] == "Lu177" and meta["weight"] == 70.0


def test_resample_labels_sidecar_copied(tmp_path):
    """A label map sidecar stays a pure {id: name} map."""
    import json
    from radot.meta import get_sidecar_path

    image = sitk.Image([8, 8, 8], sitk.sitkUInt8)
    image.SetSpacing((4.0, 4.0, 4.0))
    path = tmp_path / "labels.nii.gz"
    sitk.WriteImage(image, str(path))
    get_sidecar_path(path).write_text(json.dumps({"1": "Liver", "2": "Spleen"}))
    out = tmp_path / "labels_2mm.nii.gz"
    resample_image_file(path, (2.0, 2.0, 2.0), image_type="labels", output_path=out)
    assert json.loads(get_sidecar_path(out).read_text()) == {"1": "Liver", "2": "Spleen"}


# Audit B4-B5: the image type and unit are never guessed nor overridden silently

@pytest.mark.parametrize("image_type", ["mask", "seg", "label", None])
def test_resample_image_unknown_or_missing_type_errors(image_type):
    image = sitk.Image([4, 4, 4], sitk.sitkUInt8)
    with pytest.raises(ValueError, match="Unknown image type|Image type is required"):
        resample_image(image, (2.0, 2.0, 2.0), image_type=image_type)


def test_resample_image_type_case_insensitive_labels_nearest():
    arr = np.zeros((4, 4, 4), dtype=np.uint8)
    arr[:, :, 2:] = 7
    image = sitk.GetImageFromArray(arr)
    out = resample_image(image, (0.5, 0.5, 0.5), image_type="Labels")
    assert set(np.unique(sitk.GetArrayFromImage(out))) <= {0, 7}


def test_resample_file_sidecar_unknown_modality_errors(tmp_path):
    path = _activity_file(tmp_path, "m.nii.gz", {"modality": "mask"})
    with pytest.raises(ValueError, match="Unknown image type 'mask'"):
        resample_image_file(path, (2.0, 2.0, 2.0), output_path=tmp_path / "o.nii.gz")


@pytest.mark.parametrize("like", [False, True])
@pytest.mark.parametrize("opt, match", [
    ({"type": "ct"}, "Image type mismatch: the sidecar says 'spect', the option says 'ct'"),
    ({"unit": "bq"}, "Unit mismatch: the sidecar says 'bqml', the option says 'bq'"),
])
def test_resample_option_contradicts_sidecar_errors(tmp_path, like, opt, match):
    """--spacing and --like: an option that contradicts the sidecar is an error (was: silently used with
    --spacing, the unit relabelled without conversion)."""
    from radot.resampling import resample_image_file_like

    path = _activity_file(tmp_path, "s.nii.gz", {"modality": "spect", "unit": "bqml"})
    out = tmp_path / "o.nii.gz"
    with pytest.raises(ValueError, match=match):
        if like:
            resample_image_file_like(path, path, out, type=opt.get("type"), unit=opt.get("unit"))
        else:
            resample_image_file(path, (2.0, 2.0, 2.0), image_type=opt.get("type"), unit=opt.get("unit"),
                                output_path=out)
    assert not out.exists()


def test_resample_option_same_as_sidecar_ok(tmp_path):
    path = _activity_file(tmp_path, "s.nii.gz", {"modality": "spect", "unit": "bqml"})
    out = tmp_path / "o.nii.gz"
    resample_image_file(path, (2.0, 2.0, 2.0), image_type="SPECT", unit="BQML", output_path=out)
    assert out.exists()
