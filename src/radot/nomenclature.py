import json
import importlib.resources
import difflib
from pathlib import Path
from typing import Dict, Tuple, List, Optional, Any


def _load_json_data(filename: str) -> dict:
    p = importlib.resources.files("radot").joinpath("data").joinpath(filename)
    with open(str(p), "r", encoding="utf-8") as f:
        return json.load(f)


# Load the raw mapping dictionaries
TS_TO_TG263: Dict[str, List] = _load_json_data("ts_to_tg263.json")
TS_TO_OD_FALLBACK: Dict[str, str] = _load_json_data("ts_to_od_fallback.json")
TS_GROUPS: Dict[str, str] = _load_json_data("ts_groups.json")
OPENDOSE_PASSIVE_TARGETS: List[str] = _load_json_data("opendose_passive_targets.json")

# Load TG-263 to OpenDose mapping and invert it for reverse lookup
TG263_TO_OD: Dict[str, str] = _load_json_data("tg263_to_od.json")
OD_TO_TG263: Dict[str, str] = {v: k for k, v in TG263_TO_OD.items()}


# Reserved label-ID ranges shared across radot and downstream consumers (e.g. pet_segment).
# These must stay non-overlapping so merged label volumes never collide.
PRIMARY_ORGAN_ID_RANGE: Tuple[int, int] = (1, 999)
"""IDs assigned to organs/vessels/muscles/bones registered in TS_TO_TG263."""

DISCRETE_TISSUE_IDS: set = {1000, 1001, 1002, 1003}
"""Fixed IDs produced by the 'tissue_4_types' task: subcutaneous_fat, torso_fat,
skeletal_muscle, intermuscular_fat (in that order)."""

FALLBACK_UNMAPPED_ID_RANGE: Tuple[int, int] = (9000, 9999)
"""Reserved for auto-assigned IDs given to TotalSegmentator classes not yet
present in TS_TO_TG263. Never reuse this range for a permanent mapping."""


def assert_no_id_range_overlap() -> None:
    """
    Validate that the reserved ID ranges/sets above do not overlap each other,
    and that no ID currently registered in TS_TO_TG263 falls outside its
    designated range. Raises AssertionError on any collision.
    """
    primary_lo, primary_hi = PRIMARY_ORGAN_ID_RANGE
    fallback_lo, fallback_hi = FALLBACK_UNMAPPED_ID_RANGE
    primary_range = set(range(primary_lo, primary_hi + 1))
    fallback_range = set(range(fallback_lo, fallback_hi + 1))

    assert primary_range.isdisjoint(DISCRETE_TISSUE_IDS), (
        "PRIMARY_ORGAN_ID_RANGE overlaps DISCRETE_TISSUE_IDS"
    )
    assert primary_range.isdisjoint(fallback_range), (
        "PRIMARY_ORGAN_ID_RANGE overlaps FALLBACK_UNMAPPED_ID_RANGE"
    )
    assert DISCRETE_TISSUE_IDS.isdisjoint(fallback_range), (
        "DISCRETE_TISSUE_IDS overlaps FALLBACK_UNMAPPED_ID_RANGE"
    )

    registered_ids = {val[0] for val in TS_TO_TG263.values()}
    offending = registered_ids - primary_range
    assert not offending, (
        f"TS_TO_TG263 contains ID(s) outside PRIMARY_ORGAN_ID_RANGE: {sorted(offending)}"
    )


assert_no_id_range_overlap()


def resolve_roi_name_to_opendose(name: str) -> str:
    """
    Resolve user-supplied ROI names (TotalSegmentator, TG-263, or custom names)
    to official OpenDose reference organ names using the 3-stage cascade.
    """
    clean_name = name.strip()
    lower_name = clean_name.lower()

    # 0. Check if already a standard OpenDose name (case-insensitive)
    for od_name in OD_TO_TG263:
        if od_name.lower() == lower_name:
            return od_name

    # 1. Stage 1: The Standard TG-263 Route
    # A. Check if name is a TotalSegmentator name
    for ts_name, val in TS_TO_TG263.items():
        if ts_name.lower() == lower_name:
            tg263_label = val[1]
            for tg_key, od_val in TG263_TO_OD.items():
                if tg_key.lower() == tg263_label.lower():
                    return od_val

    # B. Check if name is a TG-263 label directly
    for tg_key, od_val in TG263_TO_OD.items():
        if tg_key.lower() == lower_name:
            return od_val

    # 2. Stage 2: The Direct TS-to-OD Bypass Bridge
    for ts_name, od_val in TS_TO_OD_FALLBACK.items():
        if ts_name.lower() == lower_name:
            return od_val

    # 3. Stage 3: Absolute Raw Retain Fallback & Safeguard
    return clean_name


def get_all_mappings() -> List[Tuple[str, str, str]]:
    """Return all explicit registered mappings: (ts_name, tg263_label, od_target)."""
    mappings = []
    for ts_name, val in TS_TO_TG263.items():
        tg263_label = val[1]
        od_target = TG263_TO_OD.get(tg263_label, "")
        if not od_target:
            od_target = TS_TO_OD_FALLBACK.get(ts_name, "")
        mappings.append((ts_name, tg263_label, od_target))
    return sorted(mappings, key=lambda x: x[0])


def find_closest_matches(query: str, n: int = 10) -> List[Tuple[str, str]]:
    """Find closest nomenclature matches for a query string using string distance."""
    query_clean = query.strip().lower()
    candidates = []

    for name in TS_TO_TG263.keys():
        candidates.append((name, "TotalSegmentator"))
    for name in TG263_TO_OD.keys():
        candidates.append((name, "TG-263"))
    for name in OD_TO_TG263.keys():
        candidates.append((name, "OpenDose"))
    for name in TS_TO_OD_FALLBACK.keys():
        if not any(c[0] == name for c in candidates):
            candidates.append((name, "TotalSegmentator"))

    scored = []
    for cand_name, source in candidates:
        ratio = difflib.SequenceMatcher(None, query_clean, cand_name.lower()).ratio()
        if ratio >= 0.4:
            scored.append((ratio, cand_name, source))

    scored.sort(key=lambda x: x[0], reverse=True)

    results = []
    seen = set()
    for ratio, name, source in scored:
        if name not in seen:
            seen.add(name)
            results.append((name, source))
            if len(results) >= n:
                break
                
    # Sort by source priority: TG-263 (1), TotalSegmentator (2), OpenDose (3)
    source_priority = {"TG-263": 1, "TotalSegmentator": 2, "OpenDose": 3}
    results.sort(key=lambda x: source_priority.get(x[1], 99))
    return results


def verify_sidecar_mappings(labels_map: Dict[int, str]) -> Dict[str, dict]:
    """Verify labels sidecar mapping, identifying mapped and unmapped dosimetry structures."""
    standard_od_names = set(OD_TO_TG263.keys()) | set(TS_TO_OD_FALLBACK.values())
    report = {"mapped": {}, "unmapped": {}}
    for lid, name in labels_map.items():
        resolved = resolve_roi_name_to_opendose(name)
        if resolved in standard_od_names:
            report["mapped"][lid] = (name, resolved)
        else:
            report["unmapped"][lid] = name
    return report
