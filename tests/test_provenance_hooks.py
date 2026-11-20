"""Hooks for the provenance (docs_wip 110 §1.1): each value or decision says where it comes from.

- metadata: sidecar, file header or in-memory image (meta.get_metadata_with_source)
- ROI name -> OpenDose organ: the nomenclature stage (nomenclature.resolve_roi_name_with_stage)
- organ without OpenDose data: an exception class, not a message (db.OrganNotFoundError)
- Teff: given or the method default, one function (dose.models.resolve_teff_h), logged once per run
- OpenDose files: origin and SHA-256 of the files read (db.opendose_files_used)
- option vs config: one merge, with the source of each value (dose.config.merge_dose_config)
"""
import json
import logging

import numpy as np
import pytest
import SimpleITK as sitk

from radot import db
from radot.db import OrganNotFoundError, get_reference_mass, get_s_value, opendose_files_used
from radot.dose import compute_roi_dose, compute_roi_dose_mtp
from radot.dose.config import merge_dose_config
from radot.dose.models import default_teff_h, log_teff_sources, resolve_teff_h
from radot.meta import get_metadata_case_insensitive, get_metadata_with_source
from radot.nomenclature import MAPPING_STAGES, resolve_roi_name_to_opendose, resolve_roi_name_with_stage

from test_dose_input_errors import _labels_arr, _mtp_study, _write
from test_dose_tia import make_image, make_labels
from radot.dose import read_results


# ----------------------------------------------------------------------------------------------
# metadata source
# ----------------------------------------------------------------------------------------------


def test_metadata_source_sidecar_header_image(tmp_path):
    image = sitk.Image([2, 2, 2], sitk.sitkFloat32)
    image.SetMetaData("radionuclide", "Y90")
    path = tmp_path / "spect.nrrd"
    sitk.WriteImage(image, str(path))
    # file header (no sidecar)
    assert get_metadata_with_source(None, "Radionuclide", path) == ("Y90", "header")
    # the sidecar wins
    (tmp_path / "spect.json").write_text(json.dumps({"radionuclide": "Lu177"}))
    assert get_metadata_with_source(image, "radionuclide", path) == ("Lu177", "sidecar")
    # in-memory image only
    assert get_metadata_with_source(image, "radionuclide") == ("Y90", "image")
    assert get_metadata_with_source(image, "weight", path) == (None, None)
    assert get_metadata_case_insensitive(image, "radionuclide", path) == "Lu177"


# ----------------------------------------------------------------------------------------------
# nomenclature stage
# ----------------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name, organ, stage",
    [
        ("Left kidney", "Left kidney C+M+P", "alias"),
        ("liver", "Liver", "opendose"),
        ("kidney_left", "Left kidney C+M+P", "ts_tg263"),
        ("Kidney_L", "Left kidney C+M+P", "tg263"),
        ("Liver_met1", "Liver_met1", "raw"),
    ],
)
def test_roi_name_stage(name, organ, stage):
    assert resolve_roi_name_with_stage(name) == (organ, stage)
    assert resolve_roi_name_to_opendose(name) == organ
    assert stage in MAPPING_STAGES


def test_roi_name_stage_ts_bypass(monkeypatch):
    # with the shipped tables every stage-2 entry is already resolved at stage 1 (same organ): the rule
    # is tested with an added entry
    from radot import nomenclature

    monkeypatch.setitem(nomenclature.TS_TO_OD_FALLBACK, "some_ts_structure", "Liver")
    assert resolve_roi_name_with_stage("Some_TS_Structure") == ("Liver", "ts_bypass")


# ----------------------------------------------------------------------------------------------
# OrganNotFoundError
# ----------------------------------------------------------------------------------------------


def test_organ_not_found_error_class():
    for call in (lambda: get_s_value("male", "Lu177", "Liver_met1", "Liver_met1"),
                 lambda: get_reference_mass("male", "Liver_met1")):
        with pytest.raises(OrganNotFoundError):
            call()
    assert issubclass(OrganNotFoundError, ValueError)


def test_mtp_organ_without_opendose_data_skipped_in_loop_error_when_requested(tmp_path, caplog):
    spects, labels = _mtp_study(tmp_path, labels_sidecar={"1": "Liver", "2": "Liver_met1"})
    with caplog.at_level(logging.WARNING, logger="radot"):
        res = compute_roi_dose_mtp(spects, labels, method="jackson2013", isotope="Lu177", labels_path=labels)
    assert "Liver" in res and "Liver_met1" not in res
    assert "Skipping dosimetry for ROI 'Liver_met1'" in caplog.text
    with pytest.raises(OrganNotFoundError):
        compute_roi_dose_mtp(spects, labels, method="jackson2013", isotope="Lu177", labels_path=labels,
                             rois=["Liver_met1"])


# ----------------------------------------------------------------------------------------------
# Teff
# ----------------------------------------------------------------------------------------------


def test_resolve_teff_given_default_and_none():
    params = {"Liver": {"teff": 60.0}, "phantom": "male"}
    assert resolve_teff_h("hanscheid2018", "Liver", params, 24.0, "Lu177") == (60.0, "given")
    assert resolve_teff_h("hanscheid2018", "Spleen", params, 24.0, "Lu177") == (48.0, "default")
    t, src = resolve_teff_h("hanscheid2017", "Spleen", params, 24.0, "Lu177")
    assert src == "default" and t == pytest.approx(db.get_physical_half_life_h("Lu177"))
    with pytest.raises(ValueError, match="teff\\) must be provided for ROI 'Spleen'"):
        resolve_teff_h("madsen2018", "Spleen", params, 24.0, "Lu177")
    assert default_teff_h("hanscheid2018", 10.0, "Lu177") == 20.0
    with pytest.raises(ValueError, match="no default"):
        default_teff_h("madsen2018", 10.0, "Lu177")


def test_teff_sources_logged_once_per_run(caplog):
    with caplog.at_level(logging.INFO, logger="radot"):
        log_teff_sources("hanscheid2018", ["Liver", "Spleen", "Kidney_L"], {"Liver": {"teff": 60.0}}, 24.0, "Lu177")
    line = [r for r in caplog.records if "Teff default" in r.message]
    assert len(line) == 1 and line[0].levelno == logging.INFO
    assert "2 ROI(s)" in line[0].message and "48 h: Spleen, Kidney_L; given for Liver" in line[0].message
    caplog.clear()
    with caplog.at_level(logging.INFO, logger="radot"):
        log_teff_sources("hanscheid2017", ["Liver"], None, 24.0, "Lu177")
    assert caplog.records[0].levelno == logging.WARNING and "physical half-life" in caplog.records[0].message


def test_dose_logs_default_teff(caplog):
    image = make_image(1e9, "bq", time_from_injection=24.0)
    with caplog.at_level(logging.INFO, logger="radot"):
        compute_roi_dose(image=image, labels=make_labels(), method="hanscheid2018", isotope="Lu177",
                         labels_map={1: "Left kidney"})
    assert "Teff default for 1 ROI(s), 2 x time from injection" in caplog.text


# ----------------------------------------------------------------------------------------------
# OpenDose files read
# ----------------------------------------------------------------------------------------------


def test_opendose_files_used_origin_and_sha256(tmp_path, monkeypatch):
    monkeypatch.setattr(db, "_opendose_files_used", {})
    db._ensure_file_downloaded("icrp_110_male", "masses.json")
    cache = db.opendose_cache_dir() / "icrp_110_male"
    cache.mkdir(parents=True, exist_ok=True)
    (cache / "y90.json.gz").write_bytes(b"not a real file")
    db._ensure_file_downloaded("icrp_110_male", "y90.json.gz")
    used = {u["file"]: u for u in opendose_files_used()}
    assert used["masses.json"]["origin"] == "shipped"
    assert used["y90.json.gz"]["origin"] == "cache"
    import hashlib

    assert used["y90.json.gz"]["sha256"] == hashlib.sha256(b"not a real file").hexdigest()


# ----------------------------------------------------------------------------------------------
# option vs config
# ----------------------------------------------------------------------------------------------


def test_merge_dose_config_sources():
    cfg = {"method": "madsen2018", "isotope": "Lu177", "rois": ["Liver"],
           "method_params": {"Liver": {"teff": 45.0}, "Spleen": {"teff": 50.0}, "phantom": "male"}}
    given = {"method": "hanscheid2018", "isotope": None, "method_params": {"Liver": {"teff": 60.0}}}
    merged, sources = merge_dose_config(given, cfg, given_source="option")
    assert merged["method"] == "hanscheid2018" and sources["method"] == "option"
    assert merged["isotope"] == "Lu177" and sources["isotope"] == "config"
    assert merged["labels_path"] is None and "labels_path" not in sources
    assert merged["method_params"] == {"Liver": {"teff": 60.0}, "Spleen": {"teff": 50.0}, "phantom": "male"}
    assert sources["method_params.Liver"] == "option" and sources["method_params.Spleen"] == "config"
    assert merge_dose_config({"method": "ldm"}, None)[1] == {"method": "argument"}
    with pytest.raises(ValueError, match="unknown keys"):
        merge_dose_config({"methd": "ldm"}, cfg)


def test_cli_option_wins_over_config(tmp_path):
    from click.testing import CliRunner
    from radot.cli import main

    labels = _write(tmp_path / "labels.nii.gz", _labels_arr(), {"1": "Liver", "2": "Spleen"})
    sc = {"modality": "spect", "unit": "bqml", "time_from_injection": 24.0, "radionuclide": "Lu177",
          "sex": "male", "decay_correction": "start"}
    spect = _write(tmp_path / "spect.nii.gz", np.where(_labels_arr() > 0, 1e4, 0).astype(np.float32), sc)
    cfg = tmp_path / "dose.json"
    cfg.write_text(json.dumps({"image": str(spect), "labels": str(labels), "method": "hanscheid2018",
                               "teff": {"Liver": 30.0}}))
    out = tmp_path / "out.json"
    result = CliRunner().invoke(main, ["dose", str(cfg), "--teff", "Liver:60", "-o", str(out)])
    assert result.exit_code == 0, result.output
    rois = read_results(out)
    assert rois["Liver"]["roi_teff_h"] == pytest.approx(60.0)
    assert rois["Spleen"]["roi_teff_h"] == pytest.approx(48.0)
