import pytest
import json
from pathlib import Path
from radot.db import (
    get_physical_half_life_h,
    normalize_phantom,
    normalize_isotope,
    resolve_organ_name,
    get_reference_mass,
    get_s_value
)


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
    assert get_s_value("male", "lu177", "left_kidney", "right_kidney") == pytest.approx(5.4199e-8)
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
    import radot.db

    # Mock download_opendose_data to prevent live network calls during tests
    def mock_download(phantom, isotope=None, force=False):
        import gzip
        from radot.db import normalize_phantom, opendose_cache_dir
        norm_phantom = normalize_phantom(phantom)
        print(f"Mock downloading for phantom: {norm_phantom}")

        dest_dir = opendose_cache_dir() / norm_phantom
        dest_dir.mkdir(parents=True, exist_ok=True)

        if isotope:
            norm_iso = isotope.lower()
            with gzip.open(dest_dir / f"{norm_iso}.json.gz", "wt", encoding="utf-8") as f:
                json.dump({"Left kidney": {"Left kidney": 1.5e-4}}, f)

        with open(dest_dir / "masses.json", "w", encoding="utf-8") as f:
            json.dump({"Left kidney": 150.0}, f)

    monkeypatch.setattr(radot.db, "download_opendose_data", mock_download)

    runner = CliRunner()

    # Test single isotope download (website, mocked)
    result = runner.invoke(main, ["opendose", "-m", "AM", "-r", "Lu177", "--force"])
    assert result.exit_code == 0
    assert "icrp_110_male" in result.output

    # Test error cases
    result = runner.invoke(main, ["opendose", "-m", "AM"])
    assert result.exit_code != 0
    assert "Missing option '-r' / '--radionuclide'" in result.output


def _mock_archive(tmp_path, monkeypatch, sha256=None):
    """A small OpenDose archive (as on Zenodo), served by a mocked urlopen."""
    import gzip
    import hashlib
    import io
    import urllib.request
    import zipfile
    import radot.db

    top = radot.db.OPENDOSE_DATA_TOP_FOLDER
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        z.writestr(f"{top}/icrp_110_female/masses.json", json.dumps({"Liver": 1400.0}))
        z.writestr(
            f"{top}/icrp_110_female/y90.json.gz",
            gzip.compress(json.dumps({"Liver": {"Liver": 7.0e-5}}).encode()),
        )
        z.writestr(f"{top}/README.txt", "test")
    content = buf.getvalue()
    monkeypatch.setattr(
        radot.db, "OPENDOSE_DATA_SHA256", sha256 or hashlib.sha256(content).hexdigest()
    )

    requested_urls = []

    def mock_urlopen(req, timeout=None):
        requested_urls.append(req.full_url)
        return io.BytesIO(content)

    monkeypatch.setattr(urllib.request, "urlopen", mock_urlopen)
    return requested_urls


def test_data_dir(monkeypatch, tmp_path):
    from radot.db import data_dir, opendose_cache_dir

    monkeypatch.setenv("RADOT_DATA_DIR", str(tmp_path / "radot_data"))
    assert data_dir() == tmp_path / "radot_data"
    assert opendose_cache_dir() == tmp_path / "radot_data" / "opendose"
    monkeypatch.delenv("RADOT_DATA_DIR")
    assert data_dir().name == "radot"


def test_ondemand_downloader(monkeypatch, tmp_path):
    """A file not shipped with radot comes from the Zenodo archive, extracted into the cache."""
    import radot.db
    from radot.db import _ensure_file_downloaded, OPENDOSE_DATA_URL

    monkeypatch.setenv("RADOT_DATA_DIR", str(tmp_path))
    requested_urls = _mock_archive(tmp_path, monkeypatch)
    monkeypatch.setattr(
        radot.db, "download_opendose_data", lambda *a, **k: pytest.fail("website used")
    )

    local_file = _ensure_file_downloaded("icrp_110_female", "y90.json.gz")
    assert local_file == tmp_path / "opendose" / "icrp_110_female" / "y90.json.gz"
    assert (tmp_path / "opendose" / "icrp_110_female" / "masses.json").exists()
    assert requested_urls == [OPENDOSE_DATA_URL]

    # Second call: from the cache, no download
    assert _ensure_file_downloaded("icrp_110_female", "y90.json.gz") == local_file
    assert len(requested_urls) == 1


def test_ondemand_downloader_checksum_and_fallback(monkeypatch, tmp_path):
    """A wrong archive checksum is refused; the OpenDose website is then used for the file."""
    import radot.db
    from radot.db import _ensure_file_downloaded, download_opendose_archive

    monkeypatch.setenv("RADOT_DATA_DIR", str(tmp_path))
    _mock_archive(tmp_path, monkeypatch, sha256="0" * 64)
    with pytest.raises(RuntimeError, match="Checksum mismatch"):
        download_opendose_archive()
    assert not (tmp_path / "opendose" / "icrp_110_female" / "y90.json.gz").exists()

    calls = []

    def mock_website(phantom, isotope=None, force=False):
        import gzip
        calls.append((phantom, isotope))
        d = tmp_path / "opendose" / phantom
        d.mkdir(parents=True, exist_ok=True)
        with gzip.open(d / f"{isotope}.json.gz", "wt") as f:
            json.dump({"Liver": {"Liver": 7.0e-5}}, f)

    monkeypatch.setattr(radot.db, "download_opendose_data", mock_website)
    local_file = _ensure_file_downloaded("icrp_110_female", "y90.json.gz")
    assert local_file.exists()
    assert calls == [("icrp_110_female", "y90")]


def test_shipped_data_found_first(monkeypatch, tmp_path):
    """masses.json is shipped with radot: no download."""
    import urllib.request
    from radot.db import _ensure_file_downloaded, opendose_package_dir

    monkeypatch.setenv("RADOT_DATA_DIR", str(tmp_path))
    monkeypatch.setattr(urllib.request, "urlopen", lambda *a, **k: pytest.fail("downloaded"))
    f = _ensure_file_downloaded("icrp_110_male", "masses.json")
    assert f == opendose_package_dir() / "icrp_110_male" / "masses.json"


def test_cli_opendose_status():
    from click.testing import CliRunner
    from radot.cli import main
    
    runner = CliRunner()
    result = runner.invoke(main, ["opendose", "status"])
    assert result.exit_code == 0
    assert "Local OpenDose Database Inventory" in result.output
    assert "Total space used:" in result.output

    # Test printing downloaded values for a given model and radionuclide
    result = runner.invoke(main, ["opendose", "list", "-m", "AM", "-r", "lu177"])
    assert result.exit_code == 0
    assert "S-values for Phantom: icrp_110_male, Radionuclide: lu177" in result.output
    assert "Left kidney" in result.output

    # Test printing reference masses for a model
    result = runner.invoke(main, ["opendose", "list", "-m", "AM"])
    assert result.exit_code == 0
    assert "Reference Masses for Phantom: icrp_110_male" in result.output
    assert "Left kidney" in result.output





