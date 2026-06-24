import json
import re
import importlib.resources
import icrp107_database
from typing import Dict, Any, Optional
from pathlib import Path

def get_physical_half_life_h(isotope: str) -> float:
    """
    Look up the physical half-life (Tp) of an isotope in hours.
    Loads data from the required 'icrp107-database' package.
    """
    iso_clean = re.sub(r'[-_\s]', '', isotope).lower()
    
    # Locate the icrp107 directory inside the installed package
    icrp_dir = importlib.resources.files("icrp107_database").joinpath("icrp107")
    
    matched_path = None
    for p in icrp_dir.iterdir():
        if p.name.endswith(".json"):
            name_without_ext = p.name[:-5]
            if re.sub(r'[-_\s]', '', name_without_ext).lower() == iso_clean:
                matched_path = p
                break
                
    if not matched_path:
        raise ValueError(f"RADOT Database: Isotope '{isotope}' is not available in the icrp107-database package.")
        
    try:
        with matched_path.open("r", encoding="utf-8") as f:
            data = json.load(f)
            # If utility.py/JSON double-serialized it as string:
            if isinstance(data, str):
                data = json.loads(data)
    except Exception as e:
        raise ValueError(f"RADOT Database: Failed to read database file for isotope '{isotope}': {e}")
        
    half_life = data.get("half_life")
    time_unit = data.get("time_unit")
    
    if half_life is None or time_unit is None:
        raise ValueError(f"RADOT Database: Isotope '{isotope}' entry in icrp107-database does not contain half_life or time_unit.")
        
    # Convert time units to hours
    unit_lower = time_unit.lower().strip()
    if unit_lower in ("y", "yr", "years", "year"):
        factor = 365.25 * 24.0
    elif unit_lower in ("d", "day", "days"):
        factor = 24.0
    elif unit_lower in ("h", "hr", "hrs", "hour", "hours"):
        factor = 1.0
    elif unit_lower in ("m", "min", "mins", "minute", "minutes"):
        factor = 1.0 / 60.0
    elif unit_lower in ("s", "sec", "secs", "second", "seconds"):
        factor = 1.0 / 3600.0
    elif unit_lower in ("ms", "msec", "milliseconds"):
        factor = 1.0 / (3600.0 * 1000.0)
    elif unit_lower in ("us", "microsecond", "microseconds"):
        factor = 1.0 / (3600.0 * 1e6)
    else:
        raise ValueError(f"RADOT Database: Unknown time unit '{time_unit}' for isotope '{isotope}'.")
        
    return float(half_life) * factor


# Caches for lazy loading OpenDose data
_masses_cache: Dict[str, Dict[str, float]] = {}
_s_values_cache: Dict[tuple, Dict[str, Dict[str, float]]] = {}

# Nomenclature dictionary mapping AAPM TG-263 to OpenDose reference names
NOMENCLATURE_MAP = {
    "kidney_l": "Left kidney",
    "kidney_r": "Right kidney",
    "kidneyl": "Left kidney",
    "kidneyr": "Right kidney",
    "left_kidney": "Left kidney",
    "right_kidney": "Right kidney",
    "leftkidney": "Left kidney",
    "rightkidney": "Right kidney",
    "leftkidneycmp": "Left kidney",
    "rightkidneycmp": "Right kidney",
    "liver": "Liver",
    "spleen": "Spleen",
    "stomach": "Stomach wall",
    "stomach_wall": "Stomach wall",
    "stomachwall": "Stomach wall",
    "bladder": "Urinary bladder wall",
    "urinary_bladder": "Urinary bladder wall",
    "urinary_bladder_wall": "Urinary bladder wall",
    "urinarybladderwall": "Urinary bladder wall",
    "lung_l": "Left lung lobe",
    "lung_r": "Right lung lobe",
    "lungl": "Left lung lobe",
    "lungr": "Right lung lobe",
    "left_lung": "Left lung lobe",
    "right_lung": "Right lung lobe",
}

def normalize_phantom(phantom: str) -> str:
    """Normalize phantom name to standard OpenDose folder names."""
    p = phantom.lower().strip()
    if p in ("am", "male", "icrp_110_male", "icrp 110 am"):
        return "icrp_110_male"
    elif p in ("af", "female", "icrp_110_female", "icrp 110 af"):
        return "icrp_110_female"
    else:
        raise ValueError(f"RADOT Database: Unknown phantom '{phantom}'. Available: 'male' (AM), 'female' (AF)")

def normalize_isotope(isotope: str) -> str:
    """Normalize isotope name to standard filename format (e.g. lu177, y90)."""
    return re.sub(r'[-_\s]', '', isotope).lower()

def resolve_organ_name(user_name: str, phantom: Optional[str] = None) -> str:
    """
    Resolve user-supplied organ names (including TG-263 nomenclature)
    to the standard OpenDose reference organ name.
    If phantom is provided, checks the organs list from that phantom's masses.json.
    """
    clean_name = user_name.strip().lower()
    
    # Try the explicit TG-263 mapping first
    if clean_name in NOMENCLATURE_MAP:
        return NOMENCLATURE_MAP[clean_name]
        
    # Standardize spaces and underscores/hyphens for fallback check
    standardized = re.sub(r'[-_\s]+', ' ', clean_name)
    
    # If phantom is provided, we can look up the exact name in masses.json case-insensitively
    if phantom is not None:
        try:
            norm_phantom = normalize_phantom(phantom)
            masses = _load_masses(norm_phantom)
            for official_name in masses.keys():
                if official_name.lower() == clean_name or official_name.lower() == standardized:
                    return official_name
        except Exception:
            pass
            
    # Fallback to simple title/capitalization or guess
    return standardized.capitalize()

def _load_masses(phantom: str) -> Dict[str, float]:
    """Lazy load masses.json for the specified phantom."""
    global _masses_cache
    if phantom not in _masses_cache:
        try:
            p_dir = importlib.resources.files("radot").joinpath("data").joinpath("opendose").joinpath(phantom)
            masses_file = p_dir.joinpath("masses.json")
            with masses_file.open("r", encoding="utf-8") as f:
                _masses_cache[phantom] = json.load(f)
        except Exception as e:
            raise FileNotFoundError(f"RADOT Database: Failed to load masses.json for phantom '{phantom}': {e}")
    return _masses_cache[phantom]

def _load_s_values(phantom: str, isotope: str) -> Dict[str, Dict[str, float]]:
    """Lazy load S-value JSON for the specified phantom and isotope."""
    global _s_values_cache
    key = (phantom, isotope)
    if key not in _s_values_cache:
        try:
            p_dir = importlib.resources.files("radot").joinpath("data").joinpath("opendose").joinpath(phantom)
            s_file = p_dir.joinpath(f"{isotope}.json")
            with s_file.open("r", encoding="utf-8") as f:
                _s_values_cache[key] = json.load(f)
        except Exception as e:
            raise FileNotFoundError(f"RADOT Database: Failed S-value matrix file for '{isotope}' on '{phantom}': {e}")
    return _s_values_cache[key]

def get_reference_mass(phantom: str, organ: str) -> float:
    """Get the reference mass for an organ in grams under the specified phantom."""
    norm_phantom = normalize_phantom(phantom)
    resolved_organ = resolve_organ_name(organ, phantom=norm_phantom)
    masses = _load_masses(norm_phantom)
    
    # Try resolved organ first, fallback case-insensitive lookup
    if resolved_organ in masses:
        return masses[resolved_organ]
        
    for k, v in masses.items():
        if k.lower() == resolved_organ.lower():
            return v
            
    raise ValueError(f"RADOT Database: Organ '{organ}' (resolved as '{resolved_organ}') not found in phantom '{phantom}' masses.")

def get_s_value(phantom: str, isotope: str, source_organ: str, target_organ: str) -> float:
    """Get the S-value from a source organ to a target organ for a given phantom and isotope."""
    norm_phantom = normalize_phantom(phantom)
    norm_isotope = normalize_isotope(isotope)
    resolved_src = resolve_organ_name(source_organ, phantom=norm_phantom)
    resolved_tgt = resolve_organ_name(target_organ, phantom=norm_phantom)
    
    s_matrix = _load_s_values(norm_phantom, norm_isotope)
    
    # Case-insensitive source lookup
    src_key = None
    for k in s_matrix.keys():
        if k.lower() == resolved_src.lower():
            src_key = k
            break
            
    if src_key is None:
        raise ValueError(f"RADOT Database: Source organ '{source_organ}' (resolved as '{resolved_src}') not found in S-value matrix.")
        
    tgt_matrix = s_matrix[src_key]
    
    # Case-insensitive target lookup
    for k, v in tgt_matrix.items():
        if k.lower() == resolved_tgt.lower():
            return v
            
    raise ValueError(f"RADOT Database: Target organ '{target_organ}' (resolved as '{resolved_tgt}') not found in S-value matrix for source '{resolved_src}'.")

def download_opendose_data(phantom: str, isotope: Optional[str] = None, group: Optional[str] = None, force: bool = False) -> None:
    """
    Placeholder downloader function for OpenDose data.
    """
    norm_phantom = normalize_phantom(phantom)
    
    isotopes_to_download = []
    if group == "core":
        isotopes_to_download = ["lu-177", "y-90", "i-131", "ac-225", "ra-223", "ga-68", "f-18"]
    elif isotope:
        isotopes_to_download = [isotope]
        
    print(f"[Placeholder] Downloading for phantom: {norm_phantom}")
    print(f"[Placeholder] Isotopes to download: {isotopes_to_download}")
    print(f"[Placeholder] Force overwrite: {force}")
