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
OD_TO_TG263: Dict[str, str] = _load_json_data("od_to_tg263.json")
TS_TO_OD_FALLBACK: Dict[str, str] = _load_json_data("ts_to_od_fallback.json")

# Invert OD_TO_TG263 to map TG-263 label -> OpenDose reference organ name
TG263_TO_OD: Dict[str, str] = {v: k for k, v in OD_TO_TG263.items()}


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


def find_closest_matches(query: str, n: int = 3) -> List[Tuple[str, str]]:
    """Find closest nomenclature matches for a query string using string distance."""
    candidates = []
    for name in TS_TO_TG263.keys():
        candidates.append((name, "TotalSegmentator"))
    for name in TG263_TO_OD.keys():
        candidates.append((name, "TG-263"))
    for name in OD_TO_TG263.keys():
        candidates.append((name, "OpenDose"))
    for name in TS_TO_OD_FALLBACK.keys():
        if (name, "TotalSegmentator") not in candidates:
            candidates.append((name, "TotalSegmentator"))

    names_list = [c[0] for c in candidates]
    matches = difflib.get_close_matches(query.strip(), names_list, n=n, cutoff=0.4)

    results = []
    for m in matches:
        for cand_name, source in candidates:
            if cand_name == m:
                results.append((cand_name, source))
                break
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
