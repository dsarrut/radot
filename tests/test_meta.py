import pytest
import SimpleITK as sitk
from pathlib import Path
from radot.meta import manage_metadata, get_metadata_case_insensitive

@pytest.mark.parametrize("filename", ["test_image.nii", "test_image.nrrd", "test_image.mhd"])
def test_metadata_write_and_read(tmp_path, filename):
    image_path = tmp_path / filename
    
    # Create dummy image
    image = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    sitk.WriteImage(image, str(image_path))
    
    # Write metadata
    manage_metadata(
        image_path=image_path,
        modality="PET",
        unit="bq",
        weight=70.5,
        activity=180.2,
        injection_date="2026-06-21 10:15:00",
        acquisition_date="2026-06-21 12:30:00"
    )
    
    # Read back image
    img = sitk.ReadImage(str(image_path))
    assert get_metadata_case_insensitive(img, "modality", image_path) == "pet"
    assert get_metadata_case_insensitive(img, "unit", image_path) == "bq"
    assert get_metadata_case_insensitive(img, "weight", image_path) == "70.5"
    assert get_metadata_case_insensitive(img, "activity", image_path) == "180.2"
    assert get_metadata_case_insensitive(img, "injection_date", image_path) == "2026-06-21 10:15:00"
    assert get_metadata_case_insensitive(img, "acquisition_date", image_path) == "2026-06-21 12:30:00"
    
    # Assert JSON sidecar exists
    from radot.meta import get_sidecar_path
    sidecar_path = get_sidecar_path(image_path)
    assert sidecar_path.exists()

def test_metadata_print(tmp_path, capsys):
    image_path = tmp_path / "test_image_print.nii"
    
    # Create dummy image
    image = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    sitk.WriteImage(image, str(image_path))
    
    # Write metadata
    manage_metadata(
        image_path=image_path,
        modality="SPECT",
        unit="bqml",
        weight=65.0,
        activity=120.0,
        injection_date="2026-06-21 10:15:00",
        acquisition_date="2026-06-21 12:30:00"
    )
    
    # Clear stdout buffer
    capsys.readouterr()
    
    # Read/Print metadata
    manage_metadata(image_path=image_path)
    
    captured = capsys.readouterr()
    assert "Image Size: (5, 5, 5)" in captured.out
    assert "Spacing: (1.0, 1.0, 1.0)" in captured.out
    assert "Origin: (0.0, 0.0, 0.0)" in captured.out
    assert "Orientation Matrix:" in captured.out
    assert "Pixel Type: 32-bit float" in captured.out
    assert "Number of Voxels: 125" in captured.out
    assert "Memory Size: 500.00 B" in captured.out
    assert "Modality: SPECT" in captured.out
    assert "Unit: bqml" in captured.out
    assert "Weight: 65.0 kg" in captured.out
    assert "Activity: 120.0 MBq" in captured.out
    assert "Injection Date: 2026-06-21 10:15:00" in captured.out
    assert "Acquisition Date: 2026-06-21 12:30:00" in captured.out
    assert "Time from injection: 2.25 hours" in captured.out

def test_metadata_date_validation(tmp_path):
    image_path = tmp_path / "date_val.nii"
    image = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    sitk.WriteImage(image, str(image_path))
    
    with pytest.raises(ValueError, match="Invalid injection_date format"):
        manage_metadata(image_path, injection_date="2026/06/21 10:15:00")
        
    with pytest.raises(ValueError, match="Invalid acquisition_date format"):
        manage_metadata(image_path, acquisition_date="2026-06-21")

def test_metadata_time_from_injection(tmp_path):
    image_path = tmp_path / "test_time_inj.nii"
    image = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    sitk.WriteImage(image, str(image_path))

    # Case 1: Only injection_date is set, calculate acquisition_date
    manage_metadata(image_path, injection_date="2026-06-21 10:15:00", time_from_injection="1:2:30")
    img = sitk.ReadImage(str(image_path))
    # 2026-06-21 10:15:00 + 1 day, 2 hours, 30 min = 2026-06-22 12:45:00
    assert get_metadata_case_insensitive(img, "injection_date", image_path) == "2026-06-21 10:15:00"
    assert get_metadata_case_insensitive(img, "acquisition_date", image_path) == "2026-06-22 12:45:00"

    # Case 2: Only acquisition_date is set (clear injection_date first), calculate injection_date
    # We clean sidecar first by writing only acquisition_date
    import json
    from radot.meta import get_sidecar_path
    sidecar_path = get_sidecar_path(image_path)
    with open(sidecar_path, "w") as f:
        json.dump({"acquisition_date": "2026-06-22 12:45:00"}, f)

    manage_metadata(image_path, time_from_injection="1:2:30")
    # 2026-06-22 12:45:00 - 1 day, 2 hours, 30 min = 2026-06-21 10:15:00
    assert get_metadata_case_insensitive(img, "injection_date", image_path) == "2026-06-21 10:15:00"
    assert get_metadata_case_insensitive(img, "acquisition_date", image_path) == "2026-06-22 12:45:00"

    # Case 3: None are present, defaults injection_date to 2000-01-01 12:00:00
    with open(sidecar_path, "w") as f:
        json.dump({}, f)
    manage_metadata(image_path, time_from_injection="0:5:15")
    # 2000-01-01 12:00:00 + 5 hours, 15 min = 2000-01-01 17:15:00
    assert get_metadata_case_insensitive(img, "injection_date", image_path) == "2000-01-01 12:00:00"
    assert get_metadata_case_insensitive(img, "acquisition_date", image_path) == "2000-01-01 17:15:00"

def test_metadata_time_from_injection_errors(tmp_path):
    image_path = tmp_path / "test_time_inj_err.nii"
    image = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    sitk.WriteImage(image, str(image_path))

    # Error if both dates are present
    with pytest.raises(ValueError, match="both.*present"):
        manage_metadata(
            image_path, 
            injection_date="2026-06-21 10:15:00", 
            acquisition_date="2026-06-22 12:45:00", 
            time_from_injection="1:2:30"
        )

def test_metadata_ct_default_unit(tmp_path):
    image_path = tmp_path / "test_ct_unit.nii"
    image = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    sitk.WriteImage(image, str(image_path))

    # Write metadata with modality CT and no unit specified
    manage_metadata(image_path, modality="ct")
    img = sitk.ReadImage(str(image_path))
    assert get_metadata_case_insensitive(img, "modality", image_path) == "ct"
    assert get_metadata_case_insensitive(img, "unit", image_path) == "hu"


def test_info_alias(tmp_path):
    from click.testing import CliRunner
    from radot.cli import main

    image_path = tmp_path / "test_info.nii"
    image = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    sitk.WriteImage(image, str(image_path))

    runner = CliRunner()
    # Call with 'info' command (writing unit)
    result_write = runner.invoke(main, ["info", str(image_path), "--unit", "bq"])
    assert result_write.exit_code == 0
    assert "Metadata saved successfully." in result_write.output

    # Call with 'info' command (reading details)
    result_read = runner.invoke(main, ["info", str(image_path)])
    assert result_read.exit_code == 0
    assert "Image Size: (5, 5, 5)" in result_read.output
    assert "Unit: bq" in result_read.output

