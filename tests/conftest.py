import pytest
import importlib.resources
import tempfile
from pathlib import Path
import shutil
import _pytest.monkeypatch

@pytest.fixture(scope="session", autouse=True)
def mock_opendose_dir_globally():
    """
    Globally mock importlib.resources.files to point to a temporary folder for all tests.
    This prevents unit tests from writing mock data into the live repository.
    """
    # Create a temporary directory for the session
    temp_dir = Path(tempfile.mkdtemp())
    
    # Define a mock package root and create data/opendose inside it
    mock_pkg_root = temp_dir / "radot"
    mock_opendose = mock_pkg_root / "data" / "opendose"
    mock_opendose.mkdir(parents=True, exist_ok=True)
    
    # Copy the actual masses.json (the reference files) from the live package
    live_db = Path(str(importlib.resources.files("radot").joinpath("data").joinpath("opendose")))
    
    if live_db.exists():
        for phantom_dir in live_db.iterdir():
            if phantom_dir.is_dir() and not phantom_dir.name.startswith("."):
                dest_phantom = mock_opendose / phantom_dir.name
                dest_phantom.mkdir(parents=True, exist_ok=True)
                masses_file = phantom_dir / "masses.json"
                if masses_file.exists():
                    shutil.copy(masses_file, dest_phantom / "masses.json")
                    
    # Now monkeypatch importlib.resources.files
    orig_files = importlib.resources.files
    
    def mock_files(package):
        if package == "radot":
            return mock_pkg_root
        return orig_files(package)
        
    mp = _pytest.monkeypatch.MonkeyPatch()
    mp.setattr(importlib.resources, "files", mock_files)
    
    yield
    
    mp.undo()
    # Cleanup session temp directory
    shutil.rmtree(temp_dir)
