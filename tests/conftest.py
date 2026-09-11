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

    # Copy the actual masses.json and parent JSON mappings from the live package
    live_data = Path(str(importlib.resources.files("radot").joinpath("data")))

    if live_data.exists():
        # Copy standard JSON maps (e.g. ts_to_tg263.json)
        for p in live_data.iterdir():
            if p.is_file() and p.name.endswith(".json"):
                shutil.copy(p, mock_pkg_root / "data" / p.name)

        # Copy masses.json reference databases
        live_db = live_data / "opendose"
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


@pytest.fixture(autouse=True)
def mock_db_data(monkeypatch):
    """Globally mock S-values and masses lookup to keep all tests 100% offline."""

    def mock_load_masses(phantom):
        if phantom == "icrp_110_male":
            return {
                "Left kidney": 152.983,
                "Right kidney": 157.011,
                "Liver": 1799.776,
                "Spleen": 149.968,
            }
        else:
            return {
                "Left kidney": 149.483,
                "Right kidney": 125.530,
                "Liver": 1399.995,
                "Spleen": 129.996,
            }

    def mock_load_s_values(phantom, isotope):
        if isotope == "lu177":
            if phantom == "icrp_110_male":
                return {
                    "Left kidney": {
                        "Left kidney": 1.5662e-4,
                        "Right kidney": 5.4199e-8,
                        "Liver": 4.6791e-8,
                        "Spleen": 3.4509e-7,
                    },
                    "Liver": {
                        "Left kidney": 4.7165e-8,
                        "Right kidney": 2.2932e-7,
                        "Liver": 1.3713e-5,
                        "Spleen": 4.1688e-8,
                    },
                    "Spleen": {
                        "Left kidney": 3.4509e-7,
                        "Right kidney": 2.1234e-7,
                        "Liver": 4.1688e-8,
                        "Spleen": 1.597e-4,
                    },
                }
            else:
                return {
                    "Left kidney": {
                        "Left kidney": 1.6034e-4,
                        "Right kidney": 9.0396e-8,
                        "Liver": 8.6437e-8,
                        "Spleen": 4.2378e-7,
                    },
                    "Liver": {
                        "Left kidney": 8.6622e-8,
                        "Right kidney": 3.0864e-7,
                        "Liver": 1.7533e-5,
                        "Spleen": 7.6519e-8,
                    },
                }
        elif isotope == "y90":
            if phantom == "icrp_110_male":
                return {
                    "Left kidney": {
                        "Left kidney": 5.1234e-6,
                        "Right kidney": 2.1234e-7,
                        "Liver": 1.9876e-7,
                        "Spleen": 1.1122e-6,
                    },
                    "Liver": {
                        "Left kidney": 2.1122e-7,
                        "Right kidney": 8.9876e-7,
                        "Liver": 7.93e-5,
                        "Spleen": 2.1234e-7,
                    },
                }
            else:
                return {
                    "Left kidney": {
                        "Left kidney": 5.9876e-6,
                        "Right kidney": 2.9876e-7,
                        "Liver": 2.1234e-7,
                        "Spleen": 1.2345e-6,
                    },
                    "Liver": {
                        "Left kidney": 2.1234e-7,
                        "Right kidney": 9.1234e-7,
                        "Liver": 6.9876e-5,
                        "Spleen": 2.9876e-7,
                    },
                }
        elif isotope == "i131":
            if phantom == "icrp_110_male":
                return {
                    "Left kidney": {
                        "Left kidney": 2.279e-4,
                        "Right kidney": 5.586e-7,
                        "Liver": 4.688e-7,
                        "Spleen": 3.317e-6,
                    },
                    "Liver": {
                        "Left kidney": 4.70e-7,
                        "Right kidney": 2.15e-6,
                        "Liver": 2.274e-5,
                        "Spleen": 4.182e-7,
                    },
                }
            else:
                return {
                    "Left kidney": {
                        "Left kidney": 2.333e-4,
                        "Right kidney": 9.283e-7,
                        "Liver": 8.613e-7,
                        "Spleen": 4.065e-6,
                    },
                    "Liver": {
                        "Left kidney": 8.65e-7,
                        "Right kidney": 2.745e-6,
                        "Liver": 2.841e-5,
                        "Spleen": 7.657e-7,
                    },
                }
        elif isotope == "ac225":
            if phantom == "icrp_110_male":
                return {
                    "Left kidney": {
                        "Left kidney": 6.201e-3,
                        "Right kidney": 2.215e-8,
                        "Liver": 1.926e-8,
                        "Spleen": 1.457e-7,
                    },
                    "Liver": {
                        "Left kidney": 1.949e-8,
                        "Right kidney": 9.514e-8,
                        "Liver": 5.272e-4,
                        "Spleen": 1.719e-8,
                    },
                }
            else:
                return {
                    "Left kidney": {
                        "Left kidney": 6.346e-3,
                        "Right kidney": 3.688e-8,
                        "Liver": 3.552e-8,
                        "Spleen": 1.802e-7,
                    },
                    "Liver": {
                        "Left kidney": 3.558e-8,
                        "Right kidney": 1.215e-7,
                        "Liver": 6.777e-4,
                        "Spleen": 3.139e-8,
                    },
                }
        return {}

    monkeypatch.setattr("radot.db._load_masses", mock_load_masses)
    monkeypatch.setattr("radot.db._load_s_values", mock_load_s_values)
