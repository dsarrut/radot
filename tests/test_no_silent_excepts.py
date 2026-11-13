"""Values that were silently dropped or kept as text by an `except: pass` (consolidation, 2026-11-13)."""
import json

import numpy as np
import pytest
import SimpleITK as sitk


def _img(tmp_path, name, sidecar=None, arr=None, dtype=np.float32):
    arr = np.ones((4, 4, 4), dtype=dtype) if arr is None else arr
    path = tmp_path / f"{name}.nii.gz"
    sitk.WriteImage(sitk.GetImageFromArray(arr), str(path))
    if sidecar is not None:
        (tmp_path / f"{name}.json").write_text(json.dumps(sidecar))
    return path


# Numbers in the sidecar: never kept as text, never dropped

def test_metadata_number():
    from radot.meta import metadata_number

    assert metadata_number("70.5", "weight") == 70.5 and metadata_number(None, "weight") is None
    with pytest.raises(ValueError, match="weight '70kg' in the metadata is not a number"):
        metadata_number("70kg", "weight")


def test_meta_write_with_non_numeric_weight_errors(tmp_path):
    from radot.meta import write_metadata

    p = _img(tmp_path, "s", {"unit": "bqml", "weight": "70kg"})
    with pytest.raises(ValueError, match="weight '70kg'"):
        write_metadata(p, height=175.0)


def test_meta_read_keeps_text_with_warning(tmp_path, caplog):
    """radot meta must still show a broken value, so that it can be fixed."""
    from radot.meta import read_metadata

    p = _img(tmp_path, "s", {"unit": "bqml", "weight": "70kg"})
    assert read_metadata(p)["metadata"]["weight"] == "70kg"
    assert "weight '70kg'" in caplog.text


def test_convert_units_non_numeric_weight_errors(tmp_path):
    from radot.convert import convert_image_units_file

    p = _img(tmp_path, "s", {"modality": "spect", "unit": "bqml", "weight": "70kg"})
    with pytest.raises(ValueError, match="weight '70kg'"):
        convert_image_units_file(p, "bq", output_path=tmp_path / "o.nii.gz")


def test_check_non_numeric_weight_fails(tmp_path):
    from radot.check import Level, run_checks

    p = _img(tmp_path, "s", {"modality": "spect", "unit": "suv", "weight": "70kg", "activity": 7400,
                             "radionuclide": "Lu177", "injection_date": "2026-01-01 10:00:00",
                             "acquisition_date": "2026-01-02 10:00:00"})
    fail = [f for f in run_checks(p) if f.code == "meta.required" and f.level == Level.FAIL]
    assert any("weight '70kg' is not a number" in f.message for f in fail)


def test_dicom_non_numeric_weight_tag_warned(tmp_path, caplog):
    from radot.dicom import convert_dicom

    dcm = tmp_path / "pet.dcm"
    img = sitk.Image([5, 5], sitk.sitkInt16)
    for k, v in {"0008|0060": "PT", "0054|1001": "BQML", "0010|1030": "unknown",
                 "0008|002a": "20260622194529"}.items():
        img.SetMetaData(k, v)
    w = sitk.ImageFileWriter()
    w.KeepOriginalImageUIDOn()
    w.SetFileName(str(dcm))
    w.Execute(img)
    convert_dicom(dcm, tmp_path / "pet.nii.gz")
    assert "weight" not in json.loads((tmp_path / "pet.json").read_text())
    assert "Patient's Weight (0010,1030) 'unknown' is not a number" in caplog.text


# Errors that propagate

def test_resolve_organ_name_bad_phantom_errors():
    from radot.db import resolve_organ_name

    with pytest.raises(ValueError, match="Unknown phantom"):
        resolve_organ_name("Liver", "child")


def test_stats_unreadable_image_sidecar_errors(tmp_path):
    from radot.stats import compute_stats

    p = _img(tmp_path, "ct", None)
    (tmp_path / "ct.json").write_text("{ broken")
    labels = sitk.GetImageFromArray(np.ones((4, 4, 4), dtype=np.uint8))
    with pytest.raises(ValueError, match="cannot be read"):
        compute_stats(labels, sitk.ReadImage(str(p)), image_type="ct", image_path=p)


# Labels sidecar: one rule (integer keys are label ids, other keys are metadata) in dose, stats, check

def _labels(tmp_path, sidecar):
    arr = np.zeros((4, 4, 4), dtype=np.uint8)
    arr[:2] = 1
    return _img(tmp_path, "labels", sidecar, arr=arr, dtype=np.uint8)


def test_labels_sidecar_with_metadata_keys_same_everywhere(tmp_path):
    from radot.check import CheckContext
    from radot.dose.compute import load_labels_map
    from radot.stats import compute_stats_file

    lp = _labels(tmp_path, {"modality": "labels", "1": "Liver"})
    assert load_labels_map(lp) == {1: "Liver"}
    assert CheckContext(activity=[], labels=lp).labels_map() == {1: "Liver"}
    ct = _img(tmp_path, "ct", {"modality": "ct", "unit": "hu"})
    out = tmp_path / "stats.json"
    compute_stats_file(labels_path=lp, image_path=ct, output_path=out)
    assert json.loads(out.read_text())["1"]["roi"] == "Liver"


def test_labels_sidecar_without_label_id_errors_in_stats(tmp_path):
    """Was: an inverted map {"Liver": 1} silently gave no names in radot stats."""
    from radot.stats import compute_stats_file

    lp = _labels(tmp_path, {"Liver": 1})
    ct = _img(tmp_path, "ct", {"modality": "ct", "unit": "hu"})
    with pytest.raises(ValueError, match="has no ROI"):
        compute_stats_file(labels_path=lp, image_path=ct, output_path=tmp_path / "s.json")


def test_check_labels_sidecar_without_label_id_fails(tmp_path):
    from radot.check import Level, run_checks

    lp = _labels(tmp_path, {"Liver": 1})
    fail = [f for f in run_checks(labels=lp) if f.code == "roi.names" and f.level == Level.FAIL]
    assert fail and "has no ROI" in fail[0].message
