import pytest
from radot.db import (
    get_physical_half_life_h,
    normalize_phantom,
    normalize_isotope,
    resolve_organ_name,
    get_reference_mass,
    get_s_value
)
@pytest.fixture(autouse=True)
def mock_db_data(monkeypatch):
    def mock_load_masses(phantom):
        if phantom == "icrp_110_male":
            return {
                "Left kidney": 152.983,
                "Right kidney": 157.011,
                "Liver": 1799.776,
                "Spleen": 149.968
            }
        else:
            return {
                "Left kidney": 149.483,
                "Right kidney": 125.530,
                "Liver": 1399.995,
                "Spleen": 129.996
            }
            
    def mock_load_s_values(phantom, isotope):
        if isotope == "lu177":
            if phantom == "icrp_110_male":
                return {
                    "Left kidney": {
                        "Left kidney": 1.3499e-7,
                        "Right kidney": 9.0850e-8,
                        "Liver": 8.6622e-8,
                        "Spleen": 4.2378e-7
                    },
                    "Liver": {
                        "Left kidney": 9.1512e-8,
                        "Right kidney": 3.3894e-7,
                        "Liver": 1.7533e-5,
                        "Spleen": 7.6519e-8
                    }
                }
            else:
                return {
                    "Left kidney": {
                        "Left kidney": 1.6034e-4,
                        "Right kidney": 9.0396e-8,
                        "Liver": 8.6437e-8,
                        "Spleen": 4.2378e-7
                    },
                    "Liver": {
                        "Left kidney": 9.1512e-8,
                        "Right kidney": 3.0825e-7,
                        "Liver": 1.7533e-5,
                        "Spleen": 7.6519e-8
                    }
                }
        elif isotope == "y90":
            if phantom == "icrp_110_male":
                return {
                    "Left kidney": {
                        "Left kidney": 5.1234e-6,
                        "Right kidney": 2.1234e-7,
                        "Liver": 1.9876e-7,
                        "Spleen": 1.1122e-6
                    },
                    "Liver": {
                        "Left kidney": 2.1122e-7,
                        "Right kidney": 8.9876e-7,
                        "Liver": 6.1234e-5,
                        "Spleen": 2.1234e-7
                    }
                }
            else:
                return {
                    "Left kidney": {
                        "Left kidney": 5.9876e-6,
                        "Right kidney": 2.9876e-7,
                        "Liver": 2.1234e-7,
                        "Spleen": 1.2345e-6
                    },
                    "Liver": {
                        "Left kidney": 2.1234e-7,
                        "Right kidney": 9.1234e-7,
                        "Liver": 6.9876e-5,
                        "Spleen": 2.9876e-7
                    }
                }
        return {}
        
    monkeypatch.setattr("radot.db._load_masses", mock_load_masses)
    monkeypatch.setattr("radot.db._load_s_values", mock_load_s_values)

def test_get_physical_half_life():
    # Test standard known isotopes
    assert get_physical_half_life_h("Lu177") == pytest.approx(159.528)
    assert get_physical_half_life_h("lu-177") == pytest.approx(159.528)
    assert get_physical_half_life_h("Y90") == pytest.approx(64.1)
    assert get_physical_half_life_h("y-90") == pytest.approx(64.1)
    assert get_physical_half_life_h("I131") == pytest.approx(192.4968)
    assert get_physical_half_life_h("i-131") == pytest.approx(192.4968)
    # Ac225 should work now that icrp107-database is required
    assert get_physical_half_life_h("Ac225") == pytest.approx(240.0)

def test_unknown_isotope_error():
    # Test that a truly unknown isotope raises ValueError
    with pytest.raises(ValueError, match="is not available in the icrp107-database package"):
        get_physical_half_life_h("Xyz999")

def test_normalize_phantom():
    assert normalize_phantom("AM") == "icrp_110_male"
    assert normalize_phantom("male") == "icrp_110_male"
    assert normalize_phantom("icrp_110_male") == "icrp_110_male"
    assert normalize_phantom("AF") == "icrp_110_female"
    assert normalize_phantom("female") == "icrp_110_female"
    assert normalize_phantom("icrp_110_female") == "icrp_110_female"
    
    with pytest.raises(ValueError, match="Unknown phantom"):
        normalize_phantom("unknown_phantom")

def test_normalize_isotope():
    assert normalize_isotope("Lu-177") == "lu177"
    assert normalize_isotope("lu177") == "lu177"
    assert normalize_isotope("Y_90") == "y90"

def test_resolve_organ_name():
    # Test TG-263 nomenclature mappings
    assert resolve_organ_name("Kidney_L") == "Left kidney"
    assert resolve_organ_name("kidney_r") == "Right kidney"
    assert resolve_organ_name("bladder") == "Urinary bladder wall"
    assert resolve_organ_name("liver") == "Liver"
    
    # Test fallback formatting
    assert resolve_organ_name("some_random_organ") == "Some random organ"

def test_get_reference_mass():
    # Check male masses from fake database
    assert get_reference_mass("male", "Kidney_L") == pytest.approx(152.983)
    assert get_reference_mass("AM", "liver") == pytest.approx(1799.776)
    
    # Check female masses from fake database
    assert get_reference_mass("female", "Kidney_R") == pytest.approx(125.530)
    assert get_reference_mass("AF", "spleen") == pytest.approx(129.996)

def test_get_s_value():
    # Check male S-values from fake matrix
    assert get_s_value("male", "lu177", "left_kidney", "right_kidney") == pytest.approx(9.0850e-8)
    assert get_s_value("AM", "Y90", "liver", "spleen") == pytest.approx(2.1234e-7)
    
    # Check female S-values from fake matrix
    assert get_s_value("female", "lu177", "left_kidney", "left_kidney") == pytest.approx(1.6034e-4)
    assert get_s_value("AF", "y90", "liver", "liver") == pytest.approx(6.9876e-5)

def test_db_errors():
    with pytest.raises(ValueError, match="not found in phantom"):
        get_reference_mass("male", "non_existent_organ")
        
    with pytest.raises(ValueError, match="not found in S-value matrix"):
        get_s_value("male", "lu177", "non_existent_organ", "left_kidney")
        
    with pytest.raises(ValueError, match="not found in S-value matrix"):
        get_s_value("male", "lu177", "left_kidney", "non_existent_organ")

def test_cli_opendose(monkeypatch):
    from click.testing import CliRunner
    from radot.cli import main
    
    runner = CliRunner()
    
    # Test single isotope download (which runs either real or fallback offline downloader)
    result = runner.invoke(main, ["opendose", "-m", "AM", "-r", "Lu177", "--force"])
    assert result.exit_code == 0
    assert "icrp_110_male" in result.output
    
    # Test error cases
    result = runner.invoke(main, ["opendose", "-m", "AM"])
    assert result.exit_code != 0
    assert "Missing option '-r' / '--radionuclide'" in result.output

    # Test upload error case (missing destination)
    result = runner.invoke(main, ["opendose", "--upload"])
    assert result.exit_code != 0
    assert "Destination path (--upload-dest) is required" in result.output

    # Test upload success case (mocking subprocess.run)
    import subprocess
    def mock_run(cmd, check=True):
        assert cmd[0] == "scp"
        assert cmd[1] == "-r"
        # Since it lists all subfolders of opendose, verify the destination is the last argument
        assert cmd[-1] == "user@server:path"
        return
    monkeypatch.setattr(subprocess, "run", mock_run)
    result = runner.invoke(main, ["opendose", "--upload", "--upload-dest", "user@server://path"])
    assert result.exit_code == 0
    assert "Uploading local database contents using" in result.output

def test_ondemand_downloader(monkeypatch, tmp_path):
    import urllib.request
    from radot.db import _ensure_file_downloaded, DEFAULT_DATABASE_URL
    
    # Temporarily set packaging resource directories to a temporary test directory
    import importlib.resources
    def mock_files(package):
        class MockTraversable:
            def joinpath(self, *parts):
                path = tmp_path
                for p in parts:
                    path = path / p
                return path
        return MockTraversable()
    monkeypatch.setattr(importlib.resources, "files", mock_files)
    
    requested_urls = []
    def mock_urlopen(req, timeout=None):
        requested_urls.append(req.full_url)
        class MockResponse:
            def read(self):
                return b'{"mock_key": 12.3}'
            def __enter__(self):
                return self
            def __exit__(self, exc_type, exc_val, exc_tb):
                pass
        return MockResponse()
    monkeypatch.setattr(urllib.request, "urlopen", mock_urlopen)
    
    # Trigger on-demand download of masses.json
    local_file = _ensure_file_downloaded("icrp_110_male", "masses.json")
    
    assert local_file.name == "masses.json"
    assert local_file.exists()
    assert len(requested_urls) == 1
    assert requested_urls[0] == f"{DEFAULT_DATABASE_URL}icrp_110_male/masses.json"



