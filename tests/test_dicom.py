import json
import pytest
import SimpleITK as sitk
from pathlib import Path
from click.testing import CliRunner
from radot.cli import main
from radot.dicom import convert_dicom

def test_convert_dicom_single_file(tmp_path):
    dcm_file = tmp_path / "pet.dcm"
    
    # Create 2D image for single file
    image = sitk.Image([5, 5], sitk.sitkInt16)
    image.SetMetaData("0008|0060", "PT")
    image.SetMetaData("0010|1030", "75.0")
    # 250 MBq = 250,000,000 Bq
    image.SetMetaData("0018|1074", "250000000")
    image.SetMetaData("0054|1001", "BQML")
    image.SetMetaData("0018|1078", "20260622191000")
    image.SetMetaData("0008|002a", "20260622194529")
    
    writer = sitk.ImageFileWriter()
    writer.KeepOriginalImageUIDOn()
    writer.SetFileName(str(dcm_file))
    writer.Execute(image)
    
    output_nii = tmp_path / "pet.nii.gz"
    output_json = tmp_path / "pet.json"
    
    runner = CliRunner()
    result = runner.invoke(main, ["dicom", str(dcm_file), "-o", str(output_nii)])
    
    assert result.exit_code == 0
    assert output_nii.exists()
    assert output_json.exists()
    
    with open(output_json, "r") as f:
        meta = json.load(f)
        
    assert meta["modality"] == "pet"
    assert meta["unit"] == "bqml"
    assert meta["weight"] == 75.0
    assert meta["activity"] == 250.0
    assert meta["injection_date"] == "2026-06-22 19:10:00"
    assert meta["acquisition_date"] == "2026-06-22 19:45:29"


def test_convert_dicom_series_directory(tmp_path):
    series_dir = tmp_path / "ct_series"
    series_dir.mkdir()
    
    # Create a 3-slice series
    for i in range(3):
        slice_img = sitk.Image([5, 5], sitk.sitkInt16)
        slice_img.SetMetaData("0008|0060", "CT")
        slice_img.SetMetaData("0020|000d", "1.2.3.4.5.1")  # StudyInstanceUID
        slice_img.SetMetaData("0020|000e", "1.2.3.4.5.1.2")  # SeriesInstanceUID
        slice_img.SetMetaData("0020|0013", str(i + 1))  # InstanceNumber
        slice_img.SetMetaData("0020|0032", f"0\\0\\{i * 2.5}")  # ImagePositionPatient
        
        # Test fallback date/time combining (Series Date + Series Time)
        slice_img.SetMetaData("0008|0021", "20260622")
        slice_img.SetMetaData("0008|0031", "194529.000000")
        
        writer = sitk.ImageFileWriter()
        writer.KeepOriginalImageUIDOn()
        writer.SetFileName(str(series_dir / f"slice_{i}.dcm"))
        writer.Execute(slice_img)
        
    output_nii = tmp_path / "ct.nii.gz"
    output_json = tmp_path / "ct.json"
    
    runner = CliRunner()
    result = runner.invoke(main, ["dicom", str(series_dir), "-o", str(output_nii)])
    
    assert result.exit_code == 0
    assert output_nii.exists()
    assert output_json.exists()
    
    with open(output_json, "r") as f:
        meta = json.load(f)
        
    assert meta["modality"] == "ct"
    assert meta["unit"] == "hu"
    assert meta["acquisition_date"] == "2026-06-22 19:45:29"
    
    # Read back image to check size/dimensions
    img = sitk.ReadImage(str(output_nii))
    assert img.GetDimension() == 3
    assert img.GetSize() == (5, 5, 3)


def test_convert_dicom_no_output(tmp_path):
    dcm_file = tmp_path / "pet_print.dcm"
    
    image = sitk.Image([5, 5], sitk.sitkInt16)
    image.SetMetaData("0008|0060", "PT")
    image.SetMetaData("0010|1030", "75.0")
    image.SetMetaData("0018|1074", "250000000")
    image.SetMetaData("0054|1001", "BQML")
    image.SetMetaData("0018|1078", "20260622191000")
    image.SetMetaData("0008|002a", "20260622194529")
    
    writer = sitk.ImageFileWriter()
    writer.KeepOriginalImageUIDOn()
    writer.SetFileName(str(dcm_file))
    writer.Execute(image)
    
    runner = CliRunner()
    result = runner.invoke(main, ["dicom", str(dcm_file)])
    
    assert result.exit_code == 0
    assert "(0008,0060) Modality: PT" in result.output
    assert "(0054,1001) Units: BQML" in result.output
    assert "(0010,1030) Patient's Weight: 75.0" in result.output
    assert "(0018,1074) Radionuclide Total Dose: 250000000" in result.output
    assert "(0018,1078) Radiopharmaceutical Start DateTime: 20260622191000" in result.output
    assert "(0008,002a) Acquisition DateTime: 20260622194529" in result.output


def test_convert_dicom_all_tags(tmp_path):
    dcm_file = tmp_path / "pet_all.dcm"
    
    image = sitk.Image([5, 5], sitk.sitkInt16)
    image.SetMetaData("0008|0060", "PT")
    image.SetMetaData("0010|1030", "75.0")
    
    writer = sitk.ImageFileWriter()
    writer.KeepOriginalImageUIDOn()
    writer.SetFileName(str(dcm_file))
    writer.Execute(image)
    
    runner = CliRunner()
    result = runner.invoke(main, ["dicom", str(dcm_file), "--all-tags"])
    
    assert result.exit_code == 0
    # Modality and Weight should be output in (GGGG,EEEE) format
    assert "(0008,0060): PT" in result.output
    assert "(0010,1030): 75.0" in result.output


def test_convert_dicom_unsupported_output_format(tmp_path):
    dcm_file = tmp_path / "pet_fail.dcm"
    
    image = sitk.Image([5, 5], sitk.sitkInt16)
    writer = sitk.ImageFileWriter()
    writer.KeepOriginalImageUIDOn()
    writer.SetFileName(str(dcm_file))
    writer.Execute(image)
    
    runner = CliRunner()
    result = runner.invoke(main, ["dicom", str(dcm_file), "-o", str(tmp_path / "fail.dcm")])
    assert result.exit_code != 0
    assert "Writing back to DICOM format is not supported" in result.output


def test_dicom_list(tmp_path):
    # Create two subfolders with different series
    folder_a = tmp_path / "series_a"
    folder_b = tmp_path / "subdir" / "series_b"
    folder_a.mkdir(parents=True)
    folder_b.mkdir(parents=True)

    # Series A: Date 2026-06-22
    slice_a = sitk.Image([5, 5], sitk.sitkInt16)
    slice_a.SetMetaData("0008|0060", "CT")
    slice_a.SetMetaData("0020|000d", "1.2.3.4.5.1")
    slice_a.SetMetaData("0020|000e", "1.2.3.4.5.1.2")
    slice_a.SetMetaData("0020|0013", "1")
    slice_a.SetMetaData("0020|0032", "0\\0\\0")
    slice_a.SetMetaData("0008|0021", "20260622")
    slice_a.SetMetaData("0008|0031", "120000.000000")
    
    writer_a = sitk.ImageFileWriter()
    writer_a.KeepOriginalImageUIDOn()
    writer_a.SetFileName(str(folder_a / "slice_1.dcm"))
    writer_a.Execute(slice_a)

    # Series B: Date 2026-06-20 (older date, should appear first)
    # We make this a 3D image (e.g. size 5x5x3) in a single file to test multi-frame slice count
    slice_b = sitk.Image([5, 5, 3], sitk.sitkInt16)
    slice_b.SetMetaData("0008|0060", "PT")
    slice_b.SetMetaData("0020|000d", "1.2.3.4.5.2")
    slice_b.SetMetaData("0020|000e", "1.2.3.4.5.2.2")
    slice_b.SetMetaData("0020|0013", "1")
    slice_b.SetMetaData("0008|0021", "20260620")
    slice_b.SetMetaData("0008|0031", "093000.000000")
    # Mark as TOMO (sinogram)
    slice_b.SetMetaData("0008|0008", "ORIGINAL\\PRIMARY\\TOMO\\EMISSION")
    
    writer_b = sitk.ImageFileWriter()
    writer_b.KeepOriginalImageUIDOn()
    writer_b.SetFileName(str(folder_b / "slice_1.dcm"))
    writer_b.Execute(slice_b)

    runner = CliRunner()
    result = runner.invoke(main, ["dicom-list", str(tmp_path)])
    assert result.exit_code == 0
    
    # Assert sorted order: Series B (2026-06-20) should appear before Series A (2026-06-22)
    output = result.output
    assert "Date/Time" in output
    assert "Modality" in output
    assert "2026-06-20 09:30:00" in output
    assert "2026-06-22 12:00:00" in output
    
    # Check order in output string
    idx_b = output.index("2026-06-20 09:30:00")
    idx_a = output.index("2026-06-22 12:00:00")
    assert idx_b < idx_a

    # Verify slices count for 3D multi-frame:
    # Series B is 3D with 3 slices in 1 file -> slices count should display 3
    # Series A is 2D with 1 slice in 1 file -> slices count should display 1
    # Check the lines of output (split by newline and find corresponding rows)
    lines = [l.strip() for l in output.split("\n") if l.strip()]
    row_b = [l for l in lines if "2026-06-20 09:30:00" in l][0]
    row_a = [l for l in lines if "2026-06-22 12:00:00" in l][0]
    
    # Extract slices column (4th column)
    assert "3" in row_b
    assert "1" in row_a

    # Verify sinogram detection:
    assert "PET (Sino)" in row_b
    assert "CT" in row_a






