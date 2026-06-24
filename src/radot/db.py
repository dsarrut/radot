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

