import importlib.util
import pytest
from pathlib import Path
import SimpleITK as sitk
from radot.segmentation import get_tg263_name, get_tg263_id, segment_ct
from radot.nomenclature import (
    PRIMARY_ORGAN_ID_RANGE,
    DISCRETE_TISSUE_IDS,
    FALLBACK_UNMAPPED_ID_RANGE,
)

requires_ts = pytest.mark.skipif(
    importlib.util.find_spec("totalsegmentator") is None,
    reason="TotalSegmentator not installed (pip install 'radot[seg]')",
)

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

@requires_ts
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
        
        spleen_id = [k for k, v in class_map[task].items() if v == "spleen"][0]
        kidney_left_id = [k for k, v in class_map[task].items() if v == "kidney_left"][0]
        
        ml_image.SetPixel((5, 5, 5), spleen_id)
        ml_image.SetPixel((6, 6, 6), kidney_left_id)
        
        sitk.WriteImage(ml_image, str(output))

    with patch("totalsegmentator.python_api.totalsegmentator", side_effect=mock_totalsegmentator) as mock_ts:
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


@requires_ts
def test_segment_ct_unmapped_class_uses_fallback_range(tmp_path, caplog):
    from unittest.mock import patch
    import json
    import logging

    ct_image = sitk.Image([10, 10, 10], sitk.sitkInt16)
    ct_path = tmp_path / "dummy_ct.nii.gz"
    sitk.WriteImage(ct_image, str(ct_path))

    output_path = tmp_path / "output_labels.nii.gz"

    def mock_totalsegmentator(input, output, task, quiet, ml=False, fast=False):
        # "body" task mixes one mapped class (body_trunc) with one class not
        # yet present in TS_TO_TG263 (body_extremities).
        from totalsegmentator.map_to_binary import class_map
        ml_image = sitk.Image([10, 10, 10], sitk.sitkUInt16)

        trunc_id = [k for k, v in class_map[task].items() if v == "body_trunc"][0]
        extremities_id = [k for k, v in class_map[task].items() if v == "body_extremities"][0]

        ml_image.SetPixel((5, 5, 5), trunc_id)
        ml_image.SetPixel((6, 6, 6), extremities_id)

        sitk.WriteImage(ml_image, str(output))

    with caplog.at_level(logging.WARNING, logger="radot"):
        with patch("totalsegmentator.python_api.totalsegmentator", side_effect=mock_totalsegmentator):
            segment_ct(ct_path, output_path, mode="body", fast=True)

    json_path = output_path.with_suffix("").with_suffix(".json")
    with open(json_path, "r") as f:
        mapping = json.load(f)

    # The mapped class keeps its registered TG-263 id, from PRIMARY_ORGAN_ID_RANGE.
    mapped_ids = {int(k) for k, v in mapping.items() if v == "External"}
    assert mapped_ids
    primary_lo, primary_hi = PRIMARY_ORGAN_ID_RANGE
    assert all(primary_lo <= i <= primary_hi for i in mapped_ids)

    # The unmapped class gets an id from the fallback range, never from
    # DISCRETE_TISSUE_IDS or PRIMARY_ORGAN_ID_RANGE.
    fallback_ids = {int(k) for k, v in mapping.items() if v == "body_extremities"}
    assert fallback_ids
    fallback_lo, fallback_hi = FALLBACK_UNMAPPED_ID_RANGE
    for fid in fallback_ids:
        assert fallback_lo <= fid <= fallback_hi
        assert fid not in DISCRETE_TISSUE_IDS
        assert not (primary_lo <= fid <= primary_hi)

    # A warning is emitted flagging the coverage gap.
    assert any(
        "TS_TO_TG263 coverage gap" in record.message and "body_extremities" in record.message
        for record in caplog.records
    )


@requires_ts
def test_segment_ct_fallback_range_exhausted_raises(tmp_path):
    from unittest.mock import patch

    ct_image = sitk.Image([5, 5, 5], sitk.sitkInt16)
    ct_path = tmp_path / "dummy_ct.nii.gz"
    sitk.WriteImage(ct_image, str(ct_path))

    output_path = tmp_path / "output_labels.nii.gz"

    def mock_totalsegmentator(input, output, task, quiet, ml=False, fast=False):
        # "body" has two classes; force both through the unmapped branch below.
        from totalsegmentator.map_to_binary import class_map
        ml_image = sitk.Image([5, 5, 5], sitk.sitkUInt16)
        ids = [k for k, v in class_map[task].items() if v in ("body_trunc", "body_extremities")]
        ml_image.SetPixel((1, 1, 1), ids[0])
        ml_image.SetPixel((2, 2, 2), ids[1])
        sitk.WriteImage(ml_image, str(output))

    fallback_lo, _fallback_hi = FALLBACK_UNMAPPED_ID_RANGE

    # Shrink the fallback range to a single slot so the second unmapped class
    # overflows it, and force both classes through the unmapped branch.
    with patch("radot.segmentation.FALLBACK_UNMAPPED_ID_RANGE", (fallback_lo, fallback_lo)), \
         patch("totalsegmentator.python_api.totalsegmentator", side_effect=mock_totalsegmentator), \
         patch("radot.segmentation.get_tg263_id", return_value=None):
        with pytest.raises(RuntimeError, match="Exhausted FALLBACK_UNMAPPED_ID_RANGE"):
            segment_ct(ct_path, output_path, mode="body", fast=True)

