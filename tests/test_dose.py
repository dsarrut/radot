import json
import pytest
import SimpleITK as sitk
from click.testing import CliRunner
from radot.cli import main
from radot.dose import compute_roi_dose

def test_compute_roi_dose_core():
    # Test compute_roi_dose function directly with dummy images
    image = sitk.Image([10, 10, 10], sitk.sitkFloat32)
    labels = sitk.Image([10, 10, 10], sitk.sitkUInt16)
    
    res = compute_roi_dose(
        image=image,
        labels=labels,
        method="madsen2018",
        isotope="Lu177",
        rois=["Liver", "Kidney_L"],
        total_activity_mbq=250.0,
        method_params={"Liver": {"teff": 45.2}}
    )
    
    assert "Liver" in res
    assert "Kidney_L" in res
    assert res["Liver"]["mean_dose_gy"] == 0.0

def test_dose_cli(tmp_path):
    image_path = tmp_path / "spect.nii.gz"
    labels_path = tmp_path / "labels.nii.gz"
    ct_path = tmp_path / "ct.nii.gz"
    output_path = tmp_path / "dose_out.json"

    # Create tiny dummy images
    sitk.WriteImage(sitk.Image([5, 5, 5], sitk.sitkFloat32), str(image_path))
    sitk.WriteImage(sitk.Image([5, 5, 5], sitk.sitkUInt16), str(labels_path))
    sitk.WriteImage(sitk.Image([5, 5, 5], sitk.sitkInt16), str(ct_path))

    runner = CliRunner()
    result = runner.invoke(main, [
        "dose",
        str(image_path),
        "-l", str(labels_path),
        "--ct", str(ct_path),
        "-m", "madsen2018",
        "--isotope", "Lu177",
        "--total-activity", "250.0",
        "--rois", "Liver,Kidney_L",
        "--teff", "Liver:45.2,Kidney_L:32.1",
        "-o", str(output_path)
    ])

    assert result.exit_code == 0
    assert output_path.exists()

    with open(output_path, "r") as f:
        data = json.load(f)

    assert "Liver" in data
    assert "Kidney_L" in data
