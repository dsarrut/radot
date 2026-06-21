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
    assert get_metadata_case_insensitive(img, "modality") == "pet"
    assert get_metadata_case_insensitive(img, "unit") == "bq"
    assert get_metadata_case_insensitive(img, "weight") == "70.5"
    assert get_metadata_case_insensitive(img, "activity") == "180.2"
    assert get_metadata_case_insensitive(img, "injection_date") == "2026-06-21 10:15:00"
    assert get_metadata_case_insensitive(img, "acquisition_date") == "2026-06-21 12:30:00"

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
    assert "Modality: SPECT" in captured.out
    assert "Unit: bqml" in captured.out
    assert "Weight: 65.0 kg" in captured.out
    assert "Activity: 120.0 MBq" in captured.out
    assert "Injection Date: 2026-06-21 10:15:00" in captured.out
    assert "Acquisition Date: 2026-06-21 12:30:00" in captured.out

def test_metadata_descrip_fallback():
    # Test descrip key fallback
    image1 = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    image1.SetMetaData("descrip", "radot:modality=pet;unit=bqml;weight=72.0;activity=150.5;injection_date=2026-06-21 10:15:00")
    
    assert get_metadata_case_insensitive(image1, "modality") == "pet"
    assert get_metadata_case_insensitive(image1, "unit") == "bqml"
    assert get_metadata_case_insensitive(image1, "weight") == "72.0"
    assert get_metadata_case_insensitive(image1, "activity") == "150.5"
    assert get_metadata_case_insensitive(image1, "injection_date") == "2026-06-21 10:15:00"

    # Test ITK_FileNotes key fallback (which ITK uses for NIfTI description when reading)
    image2 = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    image2.SetMetaData("ITK_FileNotes", "radot:modality=spect;unit=bq;weight=65.5;activity=90.0;acquisition_date=2026-06-21 12:30:00")
    
    assert get_metadata_case_insensitive(image2, "modality") == "spect"
    assert get_metadata_case_insensitive(image2, "unit") == "bq"
    assert get_metadata_case_insensitive(image2, "weight") == "65.5"
    assert get_metadata_case_insensitive(image2, "activity") == "90.0"
    assert get_metadata_case_insensitive(image2, "acquisition_date") == "2026-06-21 12:30:00"

def test_metadata_date_validation(tmp_path):
    image_path = tmp_path / "date_val.nii"
    image = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    sitk.WriteImage(image, str(image_path))
    
    with pytest.raises(ValueError, match="Invalid injection_date format"):
        manage_metadata(image_path, injection_date="2026/06/21 10:15:00")
        
    with pytest.raises(ValueError, match="Invalid acquisition_date format"):
        manage_metadata(image_path, acquisition_date="2026-06-21")
