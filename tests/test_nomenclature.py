import json
import pytest
from click.testing import CliRunner
from radot.cli import main
from radot.nomenclature import (
    resolve_roi_name_to_opendose,
    get_all_mappings,
    find_closest_matches,
    verify_sidecar_mappings,
)

def test_resolve_roi_name_to_opendose():
    # Stage 1: TS to TG-263 to OpenDose
    assert resolve_roi_name_to_opendose("kidney_left") == "Left kidney"
    assert resolve_roi_name_to_opendose("Kidney_L") == "Left kidney"
    assert resolve_roi_name_to_opendose("Left kidney") == "Left kidney"

    # Stage 2: Bypass TG-263 (Direct TS to OpenDose fallback)
    assert resolve_roi_name_to_opendose("adrenal_gland_left") == "Adrenals"
    assert resolve_roi_name_to_opendose("rib_left_5") == "Skeleton"
    assert resolve_roi_name_to_opendose("autochthon_left") == "Muscle"

    # Stage 3: Raw retain fallback
    assert resolve_roi_name_to_opendose("Tumor_Lesion_1") == "Tumor_Lesion_1"


def test_get_all_mappings():
    mappings = get_all_mappings()
    assert len(mappings) > 0
    # verify standard structure tuple (ts_name, tg263, od_name)
    assert mappings[0][0] is not None
    assert mappings[0][1] is not None


def test_find_closest_matches():
    matches = find_closest_matches("kidney_l")
    assert len(matches) > 0
    match_names = [m[0] for m in matches]
    assert "kidney_left" in match_names or "Kidney_L" in match_names


def test_verify_sidecar_mappings():
    labels_map = {
        1: "kidney_left",
        2: "adrenal_gland_right",
        3: "custom_tumor_lesion"
    }
    report = verify_sidecar_mappings(labels_map)
    assert 1 in report["mapped"]
    assert report["mapped"][1] == ("kidney_left", "Left kidney")
    assert 2 in report["mapped"]
    assert report["mapped"][2] == ("adrenal_gland_right", "Adrenals")
    assert 3 in report["unmapped"]
    assert report["unmapped"][3] == "custom_tumor_lesion"


def test_cli_roi_list():
    runner = CliRunner()
    result = runner.invoke(main, ["roi", "list"])
    assert result.exit_code == 0
    assert "kidney_left" in result.output
    assert "Kidney_L" in result.output

    # Source TS option
    result_ts = runner.invoke(main, ["roi", "list", "--source", "ts"])
    assert result_ts.exit_code == 0
    assert "kidney_left" in result_ts.output


def test_cli_roi_match():
    runner = CliRunner()
    result_exact = runner.invoke(main, ["roi", "match", "kidney_left"])
    assert result_exact.exit_code == 0
    assert "is registered" in result_exact.output
    assert "Left kidney" in result_exact.output

    result_fuzzy = runner.invoke(main, ["roi", "match", "kiney_l"])
    assert result_fuzzy.exit_code == 0
    assert "is not explicitly registered" in result_fuzzy.output
    assert "Kidney_L" in result_fuzzy.output or "kidney_left" in result_fuzzy.output


def test_cli_roi_verify(tmp_path):
    runner = CliRunner()
    sidecar_path = tmp_path / "sidecar.json"
    
    # 1. Successful verification
    with open(sidecar_path, "w") as f:
        json.dump({"1": "kidney_left", "2": "liver"}, f)
        
    result_ok = runner.invoke(main, ["roi", "verify", str(sidecar_path)])
    assert result_ok.exit_code == 0
    assert "Verification successful" in result_ok.output

    # 2. Failed verification
    with open(sidecar_path, "w") as f:
        json.dump({"1": "kidney_left", "3": "my_custom_metastasis"}, f)
        
    result_fail = runner.invoke(main, ["roi", "verify", str(sidecar_path)])
    assert result_fail.exit_code != 0
    assert "Verification failed" in result_fail.output
    assert "my_custom_metastasis" in result_fail.output
