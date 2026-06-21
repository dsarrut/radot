import pytest
import numpy as np
import SimpleITK as sitk
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
    arr = np.full((5, 5, 5), 1000.0, dtype=np.float32) # 1000 Bq/mL
    image = sitk.GetImageFromArray(arr)
    image.SetSpacing((1.0, 1.0, 1.0))

    # Weight = 75 kg = 75000 g
    # Activity = 250 MBq = 250,000,000 Bq
    # Expected SUV factor = 250,000,000 / 75,000 = 3333.333 Bq/g
    # SUV = 1000 / 3333.333 = 0.3
    converted = convert_image_units(
        image,
        from_unit="bqml",
        to_unit="suv",
        weight=75.0,
        activity=250.0
    )
    res_arr = sitk.GetArrayFromImage(converted)
    assert res_arr[0, 0, 0] == pytest.approx(0.3)

    # Convert back
    back = convert_image_units(
        converted,
        from_unit="suv",
        to_unit="bqml",
        weight=75.0,
        activity=250.0
    )
    back_arr = sitk.GetArrayFromImage(back)
    assert back_arr[0, 0, 0] == pytest.approx(1000.0)


def test_convert_validation():
    image = sitk.Image(5, 5, 5, sitk.sitkFloat32)
    
    # Missing weight/activity for SUV conversion
    with pytest.raises(ValueError, match="weight.*required"):
        convert_image_units(image, from_unit="bqml", to_unit="suv", weight=None, activity=250.0)
        
    with pytest.raises(ValueError, match="activity.*required"):
        convert_image_units(image, from_unit="bqml", to_unit="suv", weight=75.0, activity=None)

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

    # Set metadata: unit=bqml, modality=pet, weight=75.0, activity=250.0
    manage_metadata(image_path, modality="pet", unit="bqml", weight=75.0, activity=250.0)

    runner = CliRunner()
    # Run conversion without specifying --from-unit, --weight, or --activity (should read from metadata)
    output_path = tmp_path / "out_suv.nii.gz"
    result = runner.invoke(main, ["convert-units", str(image_path), "-t", "suv", "-o", str(output_path)])
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
