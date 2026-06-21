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

