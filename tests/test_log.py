import logging
import pytest
from click.testing import CliRunner
import radot
from radot.cli import main
from radot.log import logger

def test_set_log_level():
    # Test setting standard levels
    radot.set_log_level("DEBUG")
    assert logger.level == logging.DEBUG
    
    radot.set_log_level("INFO")
    assert logger.level == logging.INFO
    
    # Test setting special QUIET level
    radot.set_log_level("QUIET")
    assert logger.level == logging.ERROR
    
    # Reset to default INFO
    radot.set_log_level("INFO")
    assert logger.level == logging.INFO


def test_cli_verbosity_flags(tmp_path):
    runner = CliRunner()
    
    # Test default level configuration
    radot.set_log_level("INFO")
    result_default = runner.invoke(main, ["hello", "--name", "Test"])
    assert result_default.exit_code == 0
    assert logger.level == logging.INFO
    
    # Test verbose sets logger to DEBUG
    result_verbose = runner.invoke(main, ["--verbose", "hello", "--name", "Test"])
    assert result_verbose.exit_code == 0
    assert logger.level == logging.DEBUG
    
    # Test quiet sets logger to ERROR
    result_quiet = runner.invoke(main, ["--quiet", "hello", "--name", "Test"])
    assert result_quiet.exit_code == 0
    assert logger.level == logging.ERROR

    # Create valid dummy files
    import SimpleITK as sitk
    import json
    image_path = tmp_path / "spect.nii.gz"
    labels_path = tmp_path / "labels.nii.gz"
    
    img = sitk.Image([2, 2, 2], sitk.sitkFloat32)
    img.SetMetaData("time_from_injection", "24.0")
    img.SetMetaData("unit", "Bq")
    sitk.WriteImage(img, str(image_path))
    
    with open(tmp_path / "spect.json", "w") as f:
        json.dump({"unit": "Bq", "time_from_injection": "24.0"}, f)
        
    lbls = sitk.Image([2, 2, 2], sitk.sitkUInt16)
    lbls[1, 1, 1] = 1
    sitk.WriteImage(lbls, str(labels_path))
    
    with open(tmp_path / "labels.json", "w") as f:
        json.dump({"1": "Liver"}, f)

    # Test subcommand-level quiet on dose command sets logger to ERROR
    result_sub_quiet = runner.invoke(main, [
        "dose",
        str(image_path),
        "-l", str(labels_path),
        "-m", "hanscheid2017",
        "--isotope", "Lu177",
        "-q"
    ])
    assert result_sub_quiet.exit_code == 0
    assert logger.level == logging.ERROR
    
    # Reset level
    radot.set_log_level("INFO")
