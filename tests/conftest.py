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
    
    # Copy the actual masses.json (the reference files) from the live package
    # so that tests that read them still function correctly with real reference data.
    live_db = Path(str(importlib.resources.files("radot").joinpath("data").joinpath("opendose")))
    
    if live_db.exists():
        for phantom_dir in live_db.iterdir():
            if phantom_dir.is_dir() and not phantom_dir.name.startswith("."):
                dest_phantom = temp_dir / phantom_dir.name
                dest_phantom.mkdir(parents=True, exist_ok=True)
                masses_file = phantom_dir / "masses.json"
                if masses_file.exists():
                    shutil.copy(masses_file, dest_phantom / "masses.json")
                    
    # Now monkeypatch importlib.resources.files
    orig_files = importlib.resources.files
    
    def mock_files(package):
        if package == "radot":
            class MockTraversable:
                def joinpath(self, *parts):
                    if len(parts) >= 2 and parts[0] == "data" and parts[1] == "opendose":
                        # Redirect to our temp directory
                        path = temp_dir
                        for p in parts[2:]:
                            path = path / p
                        return path
                    # Otherwise return original path
                    return orig_files(package).joinpath(*parts)
            return MockTraversable()
        return orig_files(package)
        
    mp = _pytest.monkeypatch.MonkeyPatch()
    mp.setattr(importlib.resources, "files", mock_files)
    
    yield
    
    mp.undo()
    # Cleanup session temp directory
    shutil.rmtree(temp_dir)
