import pytest
from pathlib import Path
import SimpleITK as sitk
from radot.segmentation import get_tg263_name, get_tg263_id, segment_ct

def test_get_tg263_name():
    # Standard organs
    assert get_tg263_name("kidney_left.nii.gz") == "Kidney_L"
    assert get_tg263_name("adrenal_gland_right.nii") == "Adrenal_R"
    
    # Crop suffix
    assert get_tg263_name("kidney_left_crop.nii.gz") == "Kidney_L"
    
    # Ribs
    assert get_tg263_name("rib_left_3.nii.gz") == "Bone_Rib_L_03"
    assert get_tg263_name("rib_right_12_crop.nii.gz") == "Bone_Rib_R_12"
    
    # Vertebrae (different prefix/suffix combinations)
    assert get_tg263_name("vertebra_C1.nii.gz") == "Bone_V_C01"
    assert get_tg263_name("vertebrae_T12.nii.gz") == "Bone_V_T12"
    assert get_tg263_name("vertebra__L5.nii.gz") == "Bone_V_L05"
    
    # Unknown
    assert get_tg263_name("unknown_organ.nii.gz") == "unknown_organ"

def test_get_tg263_id():
    # Standard organs
    assert get_tg263_id("kidney_left.nii.gz") == 5
    assert get_tg263_id("adrenal_gland_right.nii") == 2
    
    # Crop suffix
    assert get_tg263_id("kidney_left_crop.nii.gz") == 5
    
    # Ribs
    assert get_tg263_id("rib_left_3.nii.gz") == 74
    assert get_tg263_id("rib_right_12_crop.nii.gz") == 93
    
    # Vertebrae
    assert get_tg263_id("vertebra_C1.nii.gz") == 94
    assert get_tg263_id("vertebrae_T12.nii.gz") == 112
    assert get_tg263_id("vertebra__L5.nii.gz") == 117
    
    # Unknown
    assert get_tg263_id("unknown_organ.nii.gz") is None

def test_segment_ct(tmp_path):
    from unittest.mock import patch
    import json

    # Create a small dummy CT image
    ct_image = sitk.Image([10, 10, 10], sitk.sitkInt16)
    ct_path = tmp_path / "dummy_ct.nii.gz"
    sitk.WriteImage(ct_image, str(ct_path))
    
    # Define output path
    output_path = tmp_path / "output_labels.nii.gz"
    
    def mock_totalsegmentator(input, output, task, quiet, ml=False, fast=False):
        # Create a single multi-label image
        from totalsegmentator.map_to_binary import class_map
        ml_image = sitk.Image([10, 10, 10], sitk.sitkUInt16)
        
        spleen_id = class_map[task]["spleen"]
        kidney_left_id = class_map[task]["kidney_left"]
        
        ml_image.SetPixel((5, 5, 5), spleen_id)
        ml_image.SetPixel((6, 6, 6), kidney_left_id)
        
        sitk.WriteImage(ml_image, str(output))

    with patch("radot.segmentation.totalsegmentator", side_effect=mock_totalsegmentator) as mock_ts:
        segment_ct(ct_path, output_path, mode="total", fast=True)
        
        assert mock_ts.call_count == 1
        kwargs = mock_ts.call_args[1]
        assert kwargs["input"] == str(ct_path)
        assert kwargs["task"] == "total"
        assert kwargs["fast"] is True
        assert kwargs["ml"] is True
        assert kwargs["quiet"] is False
    
    # Verify the output files exist
    assert output_path.exists()
    
    json_path = output_path.with_suffix("").with_suffix(".json")
    assert json_path.exists()
    
    # Verify JSON content maps correctly
    with open(json_path, "r") as f:
        mapping = json.load(f)
    assert mapping == {"5": "Kidney_L", "12": "Spleen"}
    
    # Verify output properties match input
    output_image = sitk.ReadImage(str(output_path))
    assert output_image.GetSize() == ct_image.GetSize()
    assert output_image.GetSpacing() == ct_image.GetSpacing()

def test_segment_ct_test_mode(tmp_path):
    from unittest.mock import patch
    import json

    # Create a small dummy CT image with 1mm spacing, size 20x20x20
    ct_image = sitk.Image([20, 20, 20], sitk.sitkInt16)
    ct_image.SetSpacing((1.0, 1.0, 1.0))
    ct_path = tmp_path / "dummy_ct_high_res.nii.gz"
    sitk.WriteImage(ct_image, str(ct_path))
    
    # Define output path
    output_path = tmp_path / "output_labels.nii.gz"
    
    def mock_totalsegmentator(input, output, task, quiet, ml=False, fast=False):
        # Verify that the input image to totalsegmentator has 5mm spacing
        input_img = sitk.ReadImage(input)
        assert input_img.GetSpacing() == (5.0, 5.0, 5.0)
        assert input_img.GetSize() == (4, 4, 4)
        
        # Create a single multi-label image matching input size
        from totalsegmentator.map_to_binary import class_map
        ml_image = sitk.Image(input_img.GetSize(), sitk.sitkUInt16)
        ml_image.CopyInformation(input_img)
        
        spleen_id = class_map[task]["spleen"]
        ml_image.SetPixel((2, 2, 2), spleen_id)
        
        sitk.WriteImage(ml_image, str(output))

    with patch("radot.segmentation.totalsegmentator", side_effect=mock_totalsegmentator) as mock_ts:
        segment_ct(ct_path, output_path, mode="total", test_mode=True)
        
        assert mock_ts.call_count == 1
    
    # Verify the output files exist
    assert output_path.exists()
    
    output_image = sitk.ReadImage(str(output_path))
    assert output_image.GetSize() == (20, 20, 20)
    assert output_image.GetSpacing() == (1.0, 1.0, 1.0)

