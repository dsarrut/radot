import json
import pytest
from pathlib import Path
import SimpleITK as sitk
import numpy as np
from radot.stats import compute_stats

def test_compute_stats_matching_grid(tmp_path):
    # Create matching labels and CT image (10x10x10, 1mm spacing)
    spacing = (1.0, 1.0, 1.0)
    
    labels_arr = np.zeros((10, 10, 10), dtype=np.uint16)
    labels_arr[3:7, 3:7, 3:7] = 5  # Label 5
    labels_image = sitk.GetImageFromArray(labels_arr)
    labels_image.SetSpacing(spacing)
    
    ct_arr = np.zeros((10, 10, 10), dtype=np.int16)
    ct_arr[3:7, 3:7, 3:7] = 100  # Mean intensity 100
    ct_image = sitk.GetImageFromArray(ct_arr)
    ct_image.SetSpacing(spacing)
    
    labels_path = tmp_path / "labels.nii.gz"
    image_path = tmp_path / "ct.nii.gz"
    sitk.WriteImage(labels_image, str(labels_path))
    sitk.WriteImage(ct_image, str(image_path))
    
    # Write a mapping JSON sidecar
    mapping = {"5": "Kidney_L"}
    with open(tmp_path / "labels.json", "w") as f:
        json.dump(mapping, f)
        
    output_path = tmp_path / "stats.json"
    
    # Run compute_stats
    compute_stats(
        labels_path=labels_path,
        image_path=image_path,
        image_type="ct",
        output_path=output_path
    )
    
    assert output_path.exists()
    with open(output_path, "r") as f:
        stats = json.load(f)
        
    assert "5" in stats
    assert stats["5"]["roi"] == "Kidney_L"
    assert stats["5"]["volume_ml"] == pytest.approx(0.064)  # 4 * 4 * 4 = 64 voxels = 0.064 mL
    assert stats["5"]["mean"] == pytest.approx(100.0)
    assert "mass_g" in stats["5"]
    assert stats["5"]["mass_g"] == pytest.approx(0.064 * 1.1)

def test_compute_stats_mismatch_grid(tmp_path):
    # Create labels with 1mm spacing
    labels_arr = np.zeros((10, 10, 10), dtype=np.uint16)
    labels_arr[4:6, 4:6, 4:6] = 12
    labels_image = sitk.GetImageFromArray(labels_arr)
    labels_image.SetSpacing((1.0, 1.0, 1.0))
    
    # Create image with 2mm spacing (size 5x5x5)
    image_arr = np.zeros((5, 5, 5), dtype=np.float32)
    image_arr[2, 2, 2] = 50.0  # Represents a value in the center
    image_image = sitk.GetImageFromArray(image_arr)
    image_image.SetSpacing((2.0, 2.0, 2.0))
    
    labels_path = tmp_path / "labels.nii.gz"
    image_path = tmp_path / "pet.nii.gz"
    sitk.WriteImage(labels_image, str(labels_path))
    sitk.WriteImage(image_image, str(image_path))
    
    output_path = tmp_path / "stats.json"
    
    # Run compute_stats (should trigger resampling PET image to labels grid)
    compute_stats(
        labels_path=labels_path,
        image_path=image_path,
        image_type="pet",
        output_path=output_path
    )
    
    assert output_path.exists()
    with open(output_path, "r") as f:
        stats = json.load(f)
        
    assert "12" in stats
    assert stats["12"]["roi"] == "unknown_12"
    # Should successfully run and compute a positive mean uptake value
    assert stats["12"]["mean"] > 0.0

def test_compute_stats_filtering(tmp_path):
    # Create matching labels containing label 5 and label 12
    labels_arr = np.zeros((10, 10, 10), dtype=np.uint16)
    labels_arr[2:4, 2:4, 2:4] = 5
    labels_arr[6:8, 6:8, 6:8] = 12
    labels_image = sitk.GetImageFromArray(labels_arr)
    labels_image.SetSpacing((1.0, 1.0, 1.0))
    
    ct_arr = np.zeros((10, 10, 10), dtype=np.int16)
    ct_image = sitk.GetImageFromArray(ct_arr)
    ct_image.SetSpacing((1.0, 1.0, 1.0))
    
    labels_path = tmp_path / "labels.nii.gz"
    image_path = tmp_path / "ct.nii.gz"
    sitk.WriteImage(labels_image, str(labels_path))
    sitk.WriteImage(ct_image, str(image_path))
    
    # Write a mapping JSON sidecar
    mapping = {"5": "Kidney_L", "12": "Spleen"}
    with open(tmp_path / "labels.json", "w") as f:
        json.dump(mapping, f)
        
    # Test filtering by labels (only label 5)
    output_path_labels = tmp_path / "stats_labels.json"
    compute_stats(
        labels_path=labels_path,
        image_path=image_path,
        image_type="ct",
        output_path=output_path_labels,
        labels=[5]
    )
    with open(output_path_labels, "r") as f:
        stats_labels = json.load(f)
    assert "5" in stats_labels
    assert "12" not in stats_labels
    
    # Test filtering by rois (only Spleen)
    output_path_rois = tmp_path / "stats_rois.json"
    compute_stats(
        labels_path=labels_path,
        image_path=image_path,
        image_type="ct",
        output_path=output_path_rois,
        rois=["Spleen"]
    )
    with open(output_path_rois, "r") as f:
        stats_rois = json.load(f)
    assert "12" in stats_rois
    assert "5" not in stats_rois

def test_compute_stats_mismatch_grid_resample_to_image(tmp_path):
    # Create labels with 1mm spacing
    labels_arr = np.zeros((10, 10, 10), dtype=np.uint16)
    labels_arr[4:6, 4:6, 4:6] = 12
    labels_image = sitk.GetImageFromArray(labels_arr)
    labels_image.SetSpacing((1.0, 1.0, 1.0))
    
    # Create image with 2mm spacing (size 5x5x5)
    image_arr = np.zeros((5, 5, 5), dtype=np.float32)
    image_arr[2, 2, 2] = 50.0
    image_image = sitk.GetImageFromArray(image_arr)
    image_image.SetSpacing((2.0, 2.0, 2.0))
    
    labels_path = tmp_path / "labels.nii.gz"
    image_path = tmp_path / "pet.nii.gz"
    sitk.WriteImage(labels_image, str(labels_path))
    sitk.WriteImage(image_image, str(image_path))
    
    output_path = tmp_path / "stats_resample_to_image.json"
    
    # Run compute_stats (should trigger resampling labels to image grid)
    compute_stats(
        labels_path=labels_path,
        image_path=image_path,
        image_type="pet",
        output_path=output_path,
        resample_to="image"
    )
    
    assert output_path.exists()
    with open(output_path, "r") as f:
        stats = json.load(f)
        
    assert "12" in stats
    assert stats["12"]["volume_ml"] > 0.0
