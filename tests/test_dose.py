import json
import pytest
import SimpleITK as sitk
from click.testing import CliRunner
from radot.cli import main
from radot.dose import compute_roi_dose

def test_compute_roi_dose_core():
    # Create dummy images with appropriate metadata tags
    image = sitk.Image([10, 10, 10], sitk.sitkFloat32)
    image.SetMetaData("time_from_injection", "24.0")
    image.SetMetaData("unit", "Bq")
    
    labels = sitk.Image([10, 10, 10], sitk.sitkUInt16)
    labels[2, 2, 2] = 1  # Mark one voxel with label 1
    
    # We call with rois=None so it doesn't filter out our "unknown_1" label
    res = compute_roi_dose(
        image=image,
        labels=labels,
        method="hanscheid2017",
        isotope="Lu177",
        rois=None,
        total_activity_mbq=250.0,
        method_params={"unknown_1": {"teff": 45.2}}
    )
    
    assert "unknown_1" in res
    # With a dummy empty image, concentration is 0, so mean dose is 0
    assert res["unknown_1"]["mean_dose_gy"] == 0.0
    assert res["unknown_1"]["roi_teff_h"] == 45.2

def test_physical_half_life_fallback():
    # Create dummy images
    image = sitk.Image([10, 10, 10], sitk.sitkFloat32)
    image.SetMetaData("time_from_injection", "24.0")
    image.SetMetaData("unit", "Bq")
    
    labels = sitk.Image([10, 10, 10], sitk.sitkUInt16)
    labels[2, 2, 2] = 1
    
    res = compute_roi_dose(
        image=image,
        labels=labels,
        method="hanscheid2017",
        isotope="Lu177",
        rois=None,
        total_activity_mbq=None,
        method_params=None  # No Teff provided, should fall back to physical half life
    )
    
    assert "unknown_1" in res
    # Physical half life of Lu177 is 159.528 hours
    assert res["unknown_1"]["roi_teff_h"] == pytest.approx(159.528)

def test_invalid_isotope_error():
    image = sitk.Image([10, 10, 10], sitk.sitkFloat32)
    image.SetMetaData("time_from_injection", "24.0")
    image.SetMetaData("unit", "Bq")
    
    labels = sitk.Image([10, 10, 10], sitk.sitkUInt16)
    
    with pytest.raises(ValueError, match="Hanscheid 2017 method .* specific to Lu177"):
        compute_roi_dose(
            image=image,
            labels=labels,
            method="hanscheid2017",
            isotope="Y90",
        )

def test_dose_cli(tmp_path):
    image_path = tmp_path / "spect.nii.gz"
    labels_path = tmp_path / "labels.nii.gz"
    ct_path = tmp_path / "ct.nii.gz"
    output_path = tmp_path / "dose_out.json"

    # Create dummy images with metadata tags
    img = sitk.Image([5, 5, 5], sitk.sitkFloat32)
    img.SetMetaData("time_from_injection", "24.0")
    img.SetMetaData("unit", "Bq")
    sitk.WriteImage(img, str(image_path))

    # Write a companion spect.json sidecar to define metadata
    with open(tmp_path / "spect.json", "w") as f:
        json.dump({"unit": "Bq", "time_from_injection": "24.0"}, f)
    
    lbls = sitk.Image([5, 5, 5], sitk.sitkUInt16)
    lbls[2, 2, 2] = 1
    sitk.WriteImage(lbls, str(labels_path))
    
    # Write a companion labels.json sidecar to define the Liver label
    labels_json_path = tmp_path / "labels.json"
    with open(labels_json_path, "w") as f:
        json.dump({"1": "Liver"}, f)
        
    ct_img = sitk.Image([5, 5, 5], sitk.sitkInt16)
    sitk.WriteImage(ct_img, str(ct_path))

    runner = CliRunner()
    result = runner.invoke(main, [
        "dose",
        str(image_path),
        "-l", str(labels_path),
        "--ct", str(ct_path),
        "-m", "hanscheid2017",
        "--isotope", "Lu177",
        "--total-activity", "250.0",
        "--rois", "Liver",
        "--teff", "Liver:45.2",
        "-o", str(output_path)
    ])

    assert result.exit_code == 0, f"CLI output: {result.output}"
    assert output_path.exists()

    with open(output_path, "r") as f:
        data = json.load(f)

    assert "Liver" in data
    assert data["Liver"]["roi_teff_h"] == 45.2
