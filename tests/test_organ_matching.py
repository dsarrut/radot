"""Organ names: exact OpenDose keys only, no fuzzy matching (on the shipped OpenDose data)."""

import gzip
import json
from pathlib import Path

import pytest

import radot
import radot.db as db
from radot.nomenclature import OD_ALIASES, TG263_TO_OD, TS_TO_OD_FALLBACK

SHIPPED = Path(radot.__file__).parent / "data" / "opendose"
MALE_ONLY = {"Prostate"}


def _masses(phantom):
    return json.loads((SHIPPED / phantom / "masses.json").read_text())


def _lu177(phantom):
    with gzip.open(SHIPPED / phantom / "lu177.json.gz", "rt") as f:
        return json.load(f)


@pytest.fixture
def real_opendose(monkeypatch):
    """The shipped masses and Lu177 S-values instead of the conftest mock."""
    monkeypatch.setattr(db, "_load_masses", _masses)
    monkeypatch.setattr(db, "_load_s_values", lambda phantom, isotope: _lu177(phantom))


@pytest.mark.parametrize("phantom", ["icrp_110_male", "icrp_110_female"])
def test_nomenclature_targets_are_opendose_keys(phantom):
    # every OpenDose name produced by the nomenclature is an exact key of the masses and S-values
    masses, s = _masses(phantom), _lu177(phantom)
    targets = set(TG263_TO_OD.values()) | set(TS_TO_OD_FALLBACK.values()) | set(OD_ALIASES.values())
    if phantom == "icrp_110_female":
        targets -= MALE_ONLY
    missing = sorted(t for t in targets if t not in masses or t not in s or t not in s[t])
    assert not missing, f"not OpenDose keys of {phantom}: {missing}"


@pytest.mark.parametrize(
    "name, organ",
    [
        ("Kidney_L", "Left kidney C+M+P"),
        ("kidney_right", "Right kidney C+M+P"),
        ("Left kidney", "Left kidney C+M+P"),
        ("Eye_R", "Eye bulb, right"),
        ("adrenal_gland_left", "Adrenal, left"),
        ("liver", "Liver"),
        ("LIVER", "Liver"),
    ],
)
def test_resolved_names(real_opendose, name, organ):
    assert db.resolve_organ_name(name, "male") == organ
    assert db.get_reference_mass("male", name) == _masses("icrp_110_male")[organ]
    assert db.get_s_value("male", "lu177", name, name) == _lu177("icrp_110_male")[organ][organ]


@pytest.mark.parametrize("name", ["Liver_met1", "Liver lesion", "Spleen_2", "Kidney", "Lung", "T", "Ri", "Small intestine"])
def test_no_fuzzy_match(real_opendose, name):
    # a lesion named after its organ, or a partial name, must not get an organ S-value or mass
    with pytest.raises(ValueError, match="not found"):
        db.get_reference_mass("male", name)
    with pytest.raises(ValueError, match="not found"):
        db.get_s_value("male", "lu177", name, name)


def test_find_matching_key_case_and_separators():
    d = {"Left kidney C+M+P": 1, "Eye bulb, left": 2}
    assert db._find_matching_key(d, "left_kidney  c+m+p") == "Left kidney C+M+P"
    assert db._find_matching_key(d, "Left kidney") is None
    assert db._find_matching_key(d, "Eye bulb") is None


@pytest.mark.parametrize(
    "name, expected",
    [
        ("kidney left", "Did you mean: Kidney_L (→ Left kidney C+M+P)?"),
        ("Splen", "Did you mean: Spleen?"),
        ("Small intestine", "Did you mean: Bowel_Sm (→ Small intestine wall), Small intestine contents?"),
    ],
)
def test_not_found_suggests_closest_names(real_opendose, name, expected):
    # suggestions in the error message only: the name is never resolved to them
    with pytest.raises(ValueError, match="Did you mean") as e:
        db.get_s_value("male", "lu177", name, name)
    assert expected in str(e.value) and "radot roi match" in str(e.value)
    with pytest.raises(ValueError, match="Did you mean"):
        db.get_reference_mass("male", name)


def test_not_found_no_organ_suggested_for_a_lesion(real_opendose):
    with pytest.raises(ValueError) as e:
        db.get_s_value("male", "lu177", "Lesion_1", "Lesion_1")
    assert "Did you mean" not in str(e.value) and "has no OpenDose organ" in str(e.value)
