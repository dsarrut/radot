"""Dose results file: {"radot_results": 1, "radot_version": ..., "rois": {...}} (dose.files)."""
import json

import pytest

from radot import __version__
from radot.dose import RESULTS_FORMAT, read_results, write_results

TABLE = {"Liver": {"mean_dose_gy": 2.0, "roi_mass_g": 1500.0, "roi_volume_ml": 1450.0}}


def test_write_read_round_trip(tmp_path):
    path = tmp_path / "sub" / "dose.json"
    write_results(TABLE, path)
    doc = json.loads(path.read_text())
    assert doc == {"radot_results": RESULTS_FORMAT, "radot_version": __version__, "rois": TABLE}
    assert read_results(path) == TABLE


def test_read_old_top_level_table_errors(tmp_path):
    path = tmp_path / "old.json"
    path.write_text(json.dumps(TABLE))
    with pytest.raises(ValueError, match="not a radot results file"):
        read_results(path)


def test_read_unknown_format_version_errors(tmp_path):
    path = tmp_path / "future.json"
    path.write_text(json.dumps({"radot_results": 99, "rois": TABLE}))
    with pytest.raises(ValueError, match="results format 99"):
        read_results(path)


def test_extra_top_level_keys_allowed(tmp_path):
    path = tmp_path / "dose.json"
    path.write_text(json.dumps({"radot_results": RESULTS_FORMAT, "rois": TABLE, "provenance": {"x": 1}}))
    assert read_results(path) == TABLE
