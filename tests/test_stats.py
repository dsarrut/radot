import json
import pytest
from pathlib import Path
import SimpleITK as sitk
import numpy as np
from radot.stats import compute_stats, compute_stats_file

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
    
    # Run compute_stats_file
    compute_stats_file(
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
    
    # Run compute_stats_file (should trigger resampling PET image to labels grid)
    compute_stats_file(
        labels_path=labels_path,
        image_path=image_path,
        image_type="pet",
        output_path=output_path,
        unit="suv"
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
    compute_stats_file(
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
    compute_stats_file(
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
    
    # Run compute_stats_file (should trigger resampling labels to image grid)
    compute_stats_file(
        labels_path=labels_path,
        image_path=image_path,
        image_type="pet",
        output_path=output_path,
        resample_to="image",
        unit="suv"
    )
    
    assert output_path.exists()
    with open(output_path, "r") as f:
        stats = json.load(f)
        
    assert "12" in stats
    assert stats["12"]["volume_ml"] > 0.0

def test_compute_stats_image_direct():
    # Test image-based compute_stats directly
    spacing = (1.0, 1.0, 1.0)
    
    labels_arr = np.zeros((10, 10, 10), dtype=np.uint16)
    labels_arr[3:7, 3:7, 3:7] = 5
    labels_image = sitk.GetImageFromArray(labels_arr)
    labels_image.SetSpacing(spacing)
    
    ct_arr = np.zeros((10, 10, 10), dtype=np.int16)
    ct_arr[3:7, 3:7, 3:7] = 200
    ct_image = sitk.GetImageFromArray(ct_arr)
    ct_image.SetSpacing(spacing)
    
    stats_ct = compute_stats(
        labels_image=labels_image,
        image=ct_image,
        image_type="ct",
        labels_map={5: "Liver"}
    )
    
    stats_spect = compute_stats(
        labels_image=labels_image,
        image=ct_image,
        image_type="spect",
        labels_map={5: "Liver"},
        unit="suv"
    )

    stats_custom = compute_stats(
        labels_image=labels_image,
        image=ct_image,
        image_type="ct",
        labels_map={5: "Liver"},
        sigma=1.0
    )
    
    assert "5" in stats_ct
    assert stats_ct["5"]["roi"] == "Liver"
    assert stats_ct["5"]["mean"] == pytest.approx(200.0)
    assert "mass_g" in stats_ct["5"]
    assert "peak" in stats_ct["5"]
    
    assert "5" in stats_spect
    assert "peak" in stats_spect["5"]
    assert stats_spect["5"]["peak"] < stats_ct["5"]["peak"]

    assert "5" in stats_custom
    assert "peak" in stats_custom["5"]
    assert stats_custom["5"]["peak"] > stats_ct["5"]["peak"]

def test_pretty_print_stats(capsys):
    from radot.stats import pretty_print_stats
    
    stats_data = {
        "5": {
            "roi": "Liver",
            "volume_ml": 12.3456,
            "mean": 150.123,
            "std": 12.3,
            "min": 100.0,
            "max": 200.0,
            "sum": 1846.51,
            "mass_g": 14.197,
            "peak": 175.2
        }
    }
    
    # Test CT pretty print
    pretty_print_stats(stats_data, "ct")
    captured = capsys.readouterr()
    assert "Mass (g)" in captured.out
    assert "Liver" in captured.out
    assert "14.197" in captured.out
    assert "175.200" in captured.out

    # Test PET pretty print
    pretty_print_stats(stats_data, "pet")
    captured = capsys.readouterr()
    assert "Mass (g)" not in captured.out
    assert "Liver" in captured.out
    assert "175.200" in captured.out

def test_compute_stats_mismatch_grid_preserve_total(tmp_path):
    # Create labels with 1mm spacing (1x1x1)
    labels_arr = np.zeros((10, 10, 10), dtype=np.uint16)
    labels_arr[4:6, 4:6, 4:6] = 12
    labels_image = sitk.GetImageFromArray(labels_arr)
    labels_image.SetSpacing((1.0, 1.0, 1.0))
    
    # Create image with 2mm spacing (2x2x2)
    image_arr = np.zeros((5, 5, 5), dtype=np.float32)
    image_arr[2, 2, 2] = 100.0
    image_image = sitk.GetImageFromArray(image_arr)
    image_image.SetSpacing((2.0, 2.0, 2.0))
    
    labels_path = tmp_path / "labels.nii.gz"
    image_path = tmp_path / "pet.nii.gz"
    sitk.WriteImage(labels_image, str(labels_path))
    sitk.WriteImage(image_image, str(image_path))
    
    # Run with unit="bqml" (should preserve concentration, so sum will not be scaled by voxel volume ratio)
    output_bqml = tmp_path / "stats_bqml.json"
    compute_stats_file(
        labels_path=labels_path,
        image_path=image_path,
        image_type="pet",
        output_path=output_bqml,
        unit="bqml"
    )
    
    # Run with unit="bq" (should preserve total activity, scaling by 1/8)
    output_bq = tmp_path / "stats_bq.json"
    compute_stats_file(
        labels_path=labels_path,
        image_path=image_path,
        image_type="pet",
        output_path=output_bq,
        unit="bq"
    )
    
    with open(output_bqml, "r") as f:
        stats_bqml = json.load(f)
    with open(output_bq, "r") as f:
        stats_bq = json.load(f)
        
    # Sum in bq should be 1/8 of the sum in bqml due to preserve_total scaling
    # Voxel volume ratio is 1.0 / 8.0 = 0.125
    assert stats_bq["12"]["sum"] == pytest.approx(stats_bqml["12"]["sum"] * 0.125, rel=1e-3)

def test_stats_cli_unit_requirement(tmp_path):
    from click.testing import CliRunner
    from radot.cli import main
    
    labels_image = sitk.Image([10, 10, 10], sitk.sitkUInt16)
    image = sitk.Image([10, 10, 10], sitk.sitkFloat32)
    
    labels_path = tmp_path / "labels.nii.gz"
    image_path = tmp_path / "pet.nii.gz"
    sitk.WriteImage(labels_image, str(labels_path))
    sitk.WriteImage(image, str(image_path))
    
    runner = CliRunner()
    # If PET and no unit is provided (and no metadata exists):
    result = runner.invoke(main, ["stats", str(labels_path), "-i", str(image_path), "-t", "pet"])
    assert result.exit_code != 0
    assert "unit option is required" in result.output
    
    # If CT and no unit is provided, it should not raise error about unit option
    result_ct = runner.invoke(main, ["stats", str(labels_path), "-i", str(image_path), "-t", "ct"])
    assert "unit option is required" not in result_ct.output

def test_stats_metadata_fallback(tmp_path):
    from click.testing import CliRunner
    from radot.cli import main
    from radot.meta import manage_metadata

    # Create dummy images with labels
    labels_arr = np.zeros((10, 10, 10), dtype=np.uint16)
    labels_arr[3:7, 3:7, 3:7] = 5
    labels_image = sitk.GetImageFromArray(labels_arr)
    labels_image.SetSpacing((1.0, 1.0, 1.0))

    image_arr = np.zeros((10, 10, 10), dtype=np.float32)
    image_arr[3:7, 3:7, 3:7] = 10.0
    image = sitk.GetImageFromArray(image_arr)
    image.SetSpacing((1.0, 1.0, 1.0))

    labels_path = tmp_path / "labels.nii.gz"
    image_path = tmp_path / "pet_meta.nii.gz"
    sitk.WriteImage(labels_image, str(labels_path))
    sitk.WriteImage(image, str(image_path))

    # Set metadata
    manage_metadata(image_path, modality="pet", unit="suv")

    # Run stats CLI without specifying type or unit (should resolve from metadata)
    runner = CliRunner()
    result = runner.invoke(main, ["stats", str(labels_path), "-i", str(image_path)])
    assert result.exit_code == 0
    assert "Mean" in result.output

def test_stats_missing_metadata_errors(tmp_path):
    from click.testing import CliRunner
    from radot.cli import main

    labels_image = sitk.Image([10, 10, 10], sitk.sitkUInt16)
    image = sitk.Image([10, 10, 10], sitk.sitkFloat32)

    labels_path = tmp_path / "labels.nii.gz"
    image_path = tmp_path / "no_meta.nii.gz"
    sitk.WriteImage(labels_image, str(labels_path))
    sitk.WriteImage(image, str(image_path))

    runner = CliRunner()
    # Missing type/modality completely
    result = runner.invoke(main, ["stats", str(labels_path), "-i", str(image_path)])
    assert result.exit_code != 0
    assert "modality type is not specified" in result.output
