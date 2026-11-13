import json
import re
import gzip
import hashlib
import importlib.resources
import logging
import os
import shutil
import sys
import tempfile
import urllib.request
import zipfile

logger = logging.getLogger("radot")
import icrp107_database
from typing import Dict, Any, Optional
from pathlib import Path


def nuclide_key(name: str) -> str:
    """Comparison key of a radionuclide name, for every usual spelling: 'Tb149', 'Tb-149', 'tb 149',
    '149Tb', '149-Tb', '^149^Tb' -> 'tb149'; 'Tc99m', 'Tc-99m', '99mTc' -> 'tc99m'. A mass number first is
    read as mass + symbol ('99Mo' -> 'mo99'), then as metastable ('99mTc' -> 'tc99m')."""
    clean = re.sub(r"[-_\s^]", "", str(name)).lower()
    m = re.fullmatch(r"(\d+)([a-z]{1,2})", clean)
    if m:
        return f"{m.group(2)}{m.group(1)}"
    m = re.fullmatch(r"(\d+)m([a-z]{1,2})", clean)
    if m:
        return f"{m.group(2)}{m.group(1)}m"
    return clean


def _icrp107_file(isotope: str):
    """ICRP 107 data file of a radionuclide (any spelling, see nuclide_key), or None."""
    key = nuclide_key(isotope)
    icrp_dir = importlib.resources.files("icrp107_database").joinpath("icrp107")
    for p in icrp_dir.iterdir():
        if p.name.endswith(".json") and nuclide_key(p.name[:-5]) == key:
            return p
    return None


def get_physical_half_life_h(isotope: str) -> float:
    """
    Look up the physical half-life (Tp) of an isotope in hours.
    Loads data from the required 'icrp107-database' package.
    """
    matched_path = _icrp107_file(isotope)

    if not matched_path:
        raise ValueError(
            f"RADOT Database: Isotope '{isotope}' is not available in the icrp107-database package."
        )

    try:
        with matched_path.open("r", encoding="utf-8") as f:
            data = json.load(f)
            # If utility.py/JSON double-serialized it as string:
            if isinstance(data, str):
                data = json.loads(data)
    except Exception as e:
        raise ValueError(
            f"RADOT Database: Failed to read database file for isotope '{isotope}': {e}"
        )

    half_life = data.get("half_life")
    time_unit = data.get("time_unit")

    if half_life is None or time_unit is None:
        raise ValueError(
            f"RADOT Database: Isotope '{isotope}' entry in icrp107-database does not contain half_life or time_unit."
        )

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
        raise ValueError(
            f"RADOT Database: Unknown time unit '{time_unit}' for isotope '{isotope}'."
        )

    return float(half_life) * factor


# Caches for lazy loading OpenDose data
_masses_cache: Dict[str, Dict[str, float]] = {}
_s_values_cache: Dict[tuple, Dict[str, Dict[str, float]]] = {}


def normalize_phantom(phantom: str) -> str:
    """Normalize phantom name to standard OpenDose folder names."""
    p = phantom.lower().strip()
    if p in ("am", "male", "icrp_110_male", "icrp 110 am"):
        return "icrp_110_male"
    elif p in ("af", "female", "icrp_110_female", "icrp 110 af"):
        return "icrp_110_female"
    else:
        raise ValueError(
            f"RADOT Database: Unknown phantom '{phantom}'. Available: 'male' (AM), 'female' (AF)"
        )


def normalize_isotope(isotope: str) -> str:
    """Normalize isotope name to standard filename format (e.g. lu177, y90; '177Lu' too)."""
    clean = nuclide_key(isotope)
    match = re.match(r"^([a-z]+[0-9]+m?)", clean)
    if match:
        return match.group(1)
    return clean


def resolve_organ_name(user_name: str, phantom: Optional[str] = None) -> str:
    """
    Resolve user-supplied organ names (including TG-263 nomenclature)
    to the standard OpenDose reference organ name.
    If phantom is provided, checks the organs list from that phantom's masses.json.
    """
    from radot.nomenclature import resolve_roi_name_to_opendose

    resolved = resolve_roi_name_to_opendose(user_name)

    # Standardize spaces and underscores/hyphens for fallback check
    standardized = re.sub(r"[-_\s]+", " ", resolved)

    # If phantom is provided, we can look up the exact name in masses.json case-insensitively
    if phantom is not None:
        masses = _load_masses(normalize_phantom(phantom))
        for official_name in masses.keys():
            if (
                official_name.lower() == resolved.lower()
                or official_name.lower() == standardized.lower()
            ):
                return official_name

    # Fallback or resolved name
    if resolved == user_name:
        return standardized.capitalize()
    return resolved


# OpenDose S-values and masses (CC BY-SA 4.0), Zenodo doi:10.5281/zenodo.23841128
OPENDOSE_DATA_URL = "https://zenodo.org/records/23841128/files/radot_opendose_v1.0.zip?download=1"
OPENDOSE_DATA_SHA256 = "a5fc2357bc65e948b204fb588973124bb8a27f90072ddb812ba6423df270b0bf"
OPENDOSE_DATA_TOP_FOLDER = "radot_opendose_v1.0"


def data_dir() -> Path:
    """
    Folder of the downloaded radot data (writable): $RADOT_DATA_DIR, else the user cache
    (~/.cache/radot on Linux, ~/Library/Caches/radot on macOS, %LOCALAPPDATA%\\radot on Windows).
    """
    env = os.environ.get("RADOT_DATA_DIR")
    if env:
        return Path(env).expanduser()
    if sys.platform == "win32":
        base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Caches"
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
    return base / "radot"


def opendose_cache_dir() -> Path:
    """Downloaded OpenDose data: <data_dir>/opendose/<phantom>/."""
    return data_dir() / "opendose"


def opendose_package_dir() -> Path:
    """OpenDose data shipped with radot (read-only: masses and Lu177)."""
    return Path(str(importlib.resources.files("radot").joinpath("data").joinpath("opendose")))


def _find_opendose_file(phantom: str, filename: str) -> Optional[Path]:
    for folder in (opendose_package_dir(), opendose_cache_dir()):
        f = folder / phantom / filename
        if f.exists():
            return f
    return None


def download_opendose_archive(force: bool = False) -> Path:
    """
    Download the OpenDose data archive (Zenodo, all phantoms and radionuclides) into the cache,
    after checking its SHA-256. Returns the cache folder.
    """
    dest = opendose_cache_dir()
    with tempfile.TemporaryDirectory() as tmp:
        zip_path = Path(tmp) / "opendose.zip"
        logger.info(f"Downloading the OpenDose data from {OPENDOSE_DATA_URL} ...")
        req = urllib.request.Request(OPENDOSE_DATA_URL, headers={"User-Agent": "radot"})
        with urllib.request.urlopen(req, timeout=60) as response, open(zip_path, "wb") as f:
            shutil.copyfileobj(response, f)
        sha256 = hashlib.sha256(zip_path.read_bytes()).hexdigest()
        if sha256 != OPENDOSE_DATA_SHA256:
            raise RuntimeError(f"Checksum mismatch for {OPENDOSE_DATA_URL}: {sha256}")
        with zipfile.ZipFile(zip_path) as z:
            for info in z.infolist():
                if info.is_dir():
                    continue
                rel = Path(info.filename).relative_to(OPENDOSE_DATA_TOP_FOLDER)
                if ".." in rel.parts or rel.is_absolute():
                    raise RuntimeError(f"Unexpected path in archive: {info.filename}")
                target = dest / rel
                if target.exists() and not force:
                    continue
                target.parent.mkdir(parents=True, exist_ok=True)
                with z.open(info) as src, open(target, "wb") as dst:
                    shutil.copyfileobj(src, dst)
    logger.info(f"OpenDose data in {dest}")
    return dest


def _ensure_file_downloaded(phantom: str, filename: str) -> Path:
    """
    Path of an OpenDose file: shipped with radot, else in the cache, else downloaded (Zenodo
    archive, then the OpenDose website for this file only).
    """
    local_file = _find_opendose_file(phantom, filename)
    if local_file is not None:
        return local_file

    try:
        download_opendose_archive()
    except Exception as e:
        archive_error = e
        logger.warning(
            f"Download of the OpenDose archive failed: {e}. Trying the OpenDose website..."
        )
    else:
        archive_error = None
        local_file = _find_opendose_file(phantom, filename)
        if local_file is not None:
            return local_file

    try:
        if filename == "masses.json":
            download_opendose_data(phantom, isotope="lu177", force=True)
        elif filename.endswith(".json.gz"):
            download_opendose_data(phantom, isotope=filename[: -len(".json.gz")], force=True)
        else:
            raise ValueError(f"Unknown database file name format: {filename}")
        local_file = _find_opendose_file(phantom, filename)
        if local_file is None:
            raise FileNotFoundError(f"OpenDose download did not produce {filename}")
    except Exception as od_err:
        raise FileNotFoundError(
            f"RADOT Database: '{filename}' for phantom '{phantom}' is not available locally "
            f"and could not be downloaded.\n"
            f"  - OpenDose archive (Zenodo): {archive_error or 'does not contain it'}\n"
            f"  - OpenDose website: {od_err}\n"
            f"  Cache folder: {opendose_cache_dir()} (set RADOT_DATA_DIR to change it)."
        )
    return local_file


def _load_masses(phantom: str) -> Dict[str, float]:
    """Lazy load masses.json for the specified phantom."""
    global _masses_cache
    if phantom not in _masses_cache:
        try:
            masses_file = _ensure_file_downloaded(phantom, "masses.json")
            with open(masses_file, "r", encoding="utf-8") as f:
                _masses_cache[phantom] = json.load(f)
        except Exception as e:
            raise FileNotFoundError(
                f"RADOT Database: Failed to load masses.json for phantom '{phantom}': {e}"
            )
    return _masses_cache[phantom]


def _load_s_values(phantom: str, isotope: str) -> Dict[str, Dict[str, float]]:
    """Lazy load S-value JSON for the specified phantom and isotope."""
    global _s_values_cache
    key = (phantom, isotope)
    if key not in _s_values_cache:
        try:
            s_file_gz = _ensure_file_downloaded(phantom, f"{isotope}.json.gz")
            with gzip.open(s_file_gz, "rt", encoding="utf-8") as f:
                _s_values_cache[key] = json.load(f)
        except Exception as e:
            raise FileNotFoundError(
                f"RADOT Database: Failed S-value matrix file for '{isotope}' on '{phantom}': {e}"
            )
    return _s_values_cache[key]


def _normalize_organ_key(name: str) -> str:
    """Case and separators only (spaces, '_', '-'): 'left_kidney C+M+P' == 'Left kidney C+M+P'."""
    return re.sub(r"[-_\s]+", " ", name.strip()).lower()


def _find_matching_key(d: dict, key: str) -> Optional[str]:
    """The key of d equal to key, ignoring case and separators; None otherwise.

    No fuzzy matching (prefix, suffix stripping): a name that is not an OpenDose organ must not
    silently get the S-value or mass of another one (e.g. 'Liver_met1' is not 'Liver'). Names
    that differ from the OpenDose keys are mapped explicitly in the nomenclature tables.
    """
    if key in d:
        return key
    norm = _normalize_organ_key(key)
    for k in d.keys():
        if _normalize_organ_key(k) == norm:
            return k
    return None


def _closest_hint(name: str, od_keys, n: int = 3, cutoff: float = 0.6) -> str:
    """' Did you mean: ...?' for a 'not found' error. Suggestions only, never applied.

    The query is compared with every accepted name (TotalSegmentator, TG-263, OpenDose); each match
    is shown once, as its TG-263 name with the OpenDose organ it leads to, or as the OpenDose name
    when no TG-263 name maps to it; each with its naming system written out. Only names leading to an organ of od_keys are suggested.
    """
    import difflib
    from radot.nomenclature import TS_TO_TG263, TG263_TO_OD, TS_TO_OD_FALLBACK

    od_keys = set(od_keys)
    od_to_tg = {}  # one TG-263 name per OpenDose organ (the first in the table)
    for tg, od in TG263_TO_OD.items():
        od_to_tg.setdefault(od, tg)

    def show(od):
        tg = od_to_tg.get(od)
        if tg is None:
            return f"{od} (OpenDose organ)"
        if tg == od:
            return f"{tg} (TG-263 and OpenDose)"
        return f"{tg} (TG-263, OpenDose organ '{od}')"

    # accepted name -> OpenDose organ it leads to
    leads = {tg: od for tg, od in TG263_TO_OD.items()}
    leads.update({ts: TG263_TO_OD.get(tg) for ts, (_, tg) in TS_TO_TG263.items()})
    leads.update(TS_TO_OD_FALLBACK)
    leads.update({od: od for od in od_keys})

    query = re.sub(r"[-_\s]+", " ", name.strip().lower())
    scored = sorted(
        (
            (difflib.SequenceMatcher(None, query, re.sub(r"[-_\s]+", " ", c.lower())).ratio(), od)
            for c, od in leads.items()
            if od in od_keys
        ),
        reverse=True,
    )
    suggestions = []
    best = scored[0][0] if scored else 0.0
    for ratio, od in scored:
        if ratio < max(cutoff, best - 0.15) or len(suggestions) == n:
            break
        if show(od) not in suggestions:
            suggestions.append(show(od))
    lesion = " A lesion or custom structure has no OpenDose organ (do not rename it to one)."
    if not suggestions:
        return lesion + " See radot roi match."
    return f" Did you mean: {', '.join(suggestions)}?{lesion} See radot roi match."


def get_reference_mass(phantom: str, organ: str) -> float:
    """Get the reference mass for an organ in grams under the specified phantom."""
    norm_phantom = normalize_phantom(phantom)
    resolved_organ = resolve_organ_name(organ, phantom=norm_phantom)
    masses = _load_masses(norm_phantom)

    match_key = _find_matching_key(masses, resolved_organ)
    if match_key is not None:
        return masses[match_key]

    raise ValueError(
        f"RADOT Database: Organ '{organ}' (resolved as '{resolved_organ}') not found in phantom '{phantom}' masses."
        + _closest_hint(organ, masses.keys())
    )


def get_s_value(
    phantom: str, isotope: str, source_organ: str, target_organ: str
) -> float:
    """Get the S-value from a source organ to a target organ for a given phantom and isotope."""
    norm_phantom = normalize_phantom(phantom)
    norm_isotope = normalize_isotope(isotope)
    resolved_src = resolve_organ_name(source_organ, phantom=norm_phantom)
    resolved_tgt = resolve_organ_name(target_organ, phantom=norm_phantom)

    s_matrix = _load_s_values(norm_phantom, norm_isotope)

    src_key = _find_matching_key(s_matrix, resolved_src)
    if src_key is None:
        raise ValueError(
            f"RADOT Database: Source organ '{source_organ}' (resolved as '{resolved_src}') not found in S-value matrix."
            + _closest_hint(source_organ, s_matrix.keys())
        )

    tgt_matrix = s_matrix[src_key]

    tgt_key = _find_matching_key(tgt_matrix, resolved_tgt)
    if tgt_key is not None:
        return tgt_matrix[tgt_key]

    raise ValueError(
        f"RADOT Database: Target organ '{target_organ}' (resolved as '{resolved_tgt}') not found in S-value matrix for source '{resolved_src}'."
        + _closest_hint(target_organ, tgt_matrix.keys())
    )


from html.parser import HTMLParser
import urllib.request
import urllib.parse


class DropdownParser(HTMLParser):
    def __init__(self, select_id_keyword: str):
        super().__init__()
        self.select_id_keyword = select_id_keyword
        self.in_select = False
        self.current_option_value = None
        self.options = {}  # value -> text

    def handle_starttag(self, tag, attrs):
        attrs_dict = dict(attrs)
        select_id = attrs_dict.get("id", "").lower()
        select_name = attrs_dict.get("name", "").lower()
        if tag == "select" and (
            self.select_id_keyword in select_id or self.select_id_keyword in select_name
        ):
            self.in_select = True
        elif tag == "option" and self.in_select:
            self.current_option_value = attrs_dict.get("value")

    def handle_endtag(self, tag):
        if tag == "select" and self.in_select:
            self.in_select = False

    def handle_data(self, data):
        if self.in_select and self.current_option_value is not None:
            text = data.strip()
            if text and "Select a " not in text and "Select " not in text:
                self.options[self.current_option_value] = text
            self.current_option_value = None


class TableParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.in_table = False
        self.in_tbody = False
        self.in_row = False
        self.in_cell = False
        self.current_cell_data = []
        self.current_row_data = []
        self.rows = []  # list of lists of strings

    def handle_starttag(self, tag, attrs):
        attrs_dict = dict(attrs)
        table_id = attrs_dict.get("id", "")
        if tag == "table" and (
            table_id == "svalues-table"
            or attrs_dict.get("id") == "div_table"
            or "table" in attrs_dict.get("class", "")
        ):
            self.in_table = True
        elif tag == "tbody" and self.in_table:
            self.in_tbody = True
        elif tag == "tr" and self.in_tbody:
            self.in_row = True
            self.current_row_data = []
        elif tag in ("td", "th") and self.in_row:
            self.in_cell = True
            self.current_cell_data = []

    def handle_endtag(self, tag):
        if tag == "table" and self.in_table:
            self.in_table = False
        elif tag == "tbody" and self.in_tbody:
            self.in_tbody = False
        elif tag == "tr" and self.in_row:
            self.in_row = False
            if self.current_row_data:
                self.rows.append(self.current_row_data)
        elif tag in ("td", "th") and self.in_cell:
            self.in_cell = False
            self.current_row_data.append("".join(self.current_cell_data).strip())

    def handle_data(self, data):
        if self.in_cell:
            self.current_cell_data.append(data)


def download_opendose_data(
    phantom: str, isotope: Optional[str] = None, force: bool = False
) -> None:
    """
    Downloads S-values and target organ reference masses from the OpenDose website (fallback
    when the Zenodo archive is not reachable, or for a radionuclide it does not contain).
    Saves the S-value matrix as <cache>/opendose/<phantom>/<isotope>.json.gz and masses.json.
    """
    norm_phantom = normalize_phantom(phantom)

    if norm_phantom == "icrp_110_male":
        phantom_url_name = "ICRP_110_AM"
    else:
        phantom_url_name = "ICRP_110_AF"

    isotopes_to_download = []
    if isotope:
        isotopes_to_download = [isotope]

    dest_dir = opendose_cache_dir() / norm_phantom
    dest_dir.mkdir(parents=True, exist_ok=True)

    url_options = ["https://opendose.polsl.pl/svalues/", "https://opendose.org/svalues"]
    html = None
    base_url = None
    for u in url_options:
        try:
            req = urllib.request.Request(u, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=10) as response:
                html = response.read().decode("utf-8")
                base_url = u
                break
        except Exception:
            continue

    if html is None:
        raise ConnectionError(
            f"RADOT Database: Failed to connect to OpenDose registry url options: {url_options}"
        )

    src_parser = DropdownParser("source")
    src_parser.feed(html)
    source_organs = src_parser.options

    iso_parser = DropdownParser("isotope")
    iso_parser.feed(html)
    isotopes_map = iso_parser.options

    for iso_name in isotopes_to_download:
        norm_iso = normalize_isotope(iso_name)
        target_file = dest_dir / f"{norm_iso}.json.gz"

        isotope_opt_text = None
        for opt_text in isotopes_map.values():
            if normalize_isotope(opt_text) == norm_iso:
                isotope_opt_text = opt_text
                break

        if not isotope_opt_text:
            logger.warning(
                f"Isotope '{iso_name}' not available on OpenDose website. Skipping."
            )
            continue

        s_matrix = {}
        masses = {}

        # Load existing data if file exists and we are not forcing a clean overwrite
        if target_file.exists() and not force:
            try:
                with gzip.open(target_file, "rt", encoding="utf-8") as f:
                    s_matrix = json.load(f)
            except Exception as e:
                logger.warning(f"Partial download {target_file} unreadable ({e}): downloaded again.")

        masses_file = dest_dir / "masses.json"
        if masses_file.exists() and not force:
            try:
                with open(masses_file, "r", encoding="utf-8") as f:
                    masses = json.load(f)
            except Exception as e:
                logger.warning(f"Partial download {masses_file} unreadable ({e}): downloaded again.")

        # If the file is already complete (has all source organs) and not forcing, we can skip
        if target_file.exists() and not force and len(s_matrix) >= len(source_organs):
            logger.info(
                f"S-value file for '{iso_name}' is already complete. Use --force to overwrite."
            )
            continue

        for src_name in source_organs.values():
            if not force and src_name in s_matrix and s_matrix[src_name]:
                # Skip fetching if organ data already exists
                continue

            logger.info(f"Fetching S-values for source organ: {src_name}...")

            phantom_part = urllib.parse.quote(phantom_url_name)
            source_part = urllib.parse.quote(src_name).replace("/", "%252F")
            isotope_part = urllib.parse.quote(isotope_opt_text)

            if not base_url.endswith("/"):
                query_url = f"{base_url}/{phantom_part}/{source_part}/{isotope_part}/"
            else:
                query_url = f"{base_url}{phantom_part}/{source_part}/{isotope_part}/"

            try:
                req_get = urllib.request.Request(
                    query_url, headers={"User-Agent": "Mozilla/5.0"}
                )
                with urllib.request.urlopen(req_get, timeout=15) as res:
                    res_html = res.read().decode("utf-8")
            except Exception as e:
                logger.error(
                    f"Failed to fetch S-values for {src_name} ({e}). Skipping organ."
                )
                continue

            tbl_parser = TableParser()
            tbl_parser.feed(res_html)

            if not tbl_parser.rows:
                logger.warning(
                    f"No S-value data parsed for source '{src_name}' (skipping)."
                )
                continue

            s_matrix[src_name] = {}
            for row in tbl_parser.rows:
                if len(row) >= 4:
                    target_name = row[0]
                    try:
                        s_val = float(row[1])
                        mass_val = float(row[3])
                        s_matrix[src_name][target_name] = s_val
                        masses[target_name] = mass_val
                    except ValueError:
                        continue

            # Save incrementally after each organ is fetched to prevent data loss on interruption
            if masses:
                with open(masses_file, "w", encoding="utf-8") as f:
                    json.dump(masses, f, indent=4)
            if s_matrix:
                with gzip.open(target_file, "wt", encoding="utf-8") as f:
                    json.dump(s_matrix, f, indent=4)

        logger.info(
            f"Successfully saved OpenDose S-values for '{iso_name}' to: {target_file}"
        )


def get_non_penetrating_energy_constant_mev(isotope: str) -> float:
    """
    Calculate the non-penetrating decay energy constant (Delta_local) in MeV
    for a chosen isotope by summing the energy * yield for non-penetrating emissions:
    alpha, beta-, beta+, auger, IE, and alpha recoil.
    """
    import importlib.resources
    import json
    import re

    matched_path = _icrp107_file(isotope)

    if not matched_path:
        raise ValueError(
            f"Isotope '{isotope}' is not available in the ICRP-107 database."
        )

    with matched_path.open("r", encoding="utf-8") as f:
        data = json.load(f)
        if isinstance(data, str):
            data = json.loads(data)

    emissions = data.get("emissions", {})
    non_penetrating_keys = ["alpha", "beta-", "beta+", "auger", "IE", "alpha recoil"]

    total_energy_mev = 0.0
    for key in non_penetrating_keys:
        em_list = emissions.get(key, [])
        for item in em_list:
            if len(item) >= 2:
                energy, yield_val = item[0], item[1]
                total_energy_mev += energy * yield_val

    return total_energy_mev


_decay_chains_cache = None


def _load_decay_chains() -> dict:
    """Load curated radioactive decay chains from package data."""
    global _decay_chains_cache
    if _decay_chains_cache is None:
        import importlib.resources
        import json

        try:
            p = importlib.resources.files("radot").joinpath("data", "decay_chains.json")
            with open(str(p), "r", encoding="utf-8") as f:
                _decay_chains_cache = json.load(f)
        except Exception as e:
            # No silent fallback to "no chain": the cascade would become single-nuclide
            raise RuntimeError(
                f"radot data file decay_chains.json cannot be read ({e}); reinstall radot."
            ) from e
    return _decay_chains_cache


def get_decay_chain(isotope: str) -> Optional[dict]:
    """
    Retrieve decay chain metadata and daughters for an isotope if available.
    Matches case-insensitively with flexible formatting (e.g. 'Ac225', 'Ac-225', '225Ac').
    """
    import re

    clean = nuclide_key(isotope)
    data = _load_decay_chains()
    chains = data.get("chains", {})
    for parent_key, chain_info in chains.items():
        if nuclide_key(parent_key) == clean:
            return chain_info
    return None


# Secular equilibrium: a daughter counts in the cascade when its half-life is below the parent half-life
# divided by this ratio (T1/2 daughter << T1/2 parent), and so are all the nuclides between them
EQUILIBRIUM_HALF_LIFE_RATIO = 10.0
# (dZ, dA) of one decay step: alpha, beta-, EC/beta+, isomeric transition, neutron, beta- n
_DECAY_STEPS = {(-2, -4), (1, 0), (-1, 0), (0, 0), (0, -1), (1, -1)}


def cascade_members(chain: dict) -> tuple[list[dict], list[dict]]:
    """
    Split the unstable nuclides of a decay chain (parent first) into those in secular equilibrium with the
    parent (counted by the cascade) and those excluded. A daughter is kept when its half-life is below
    T1/2 parent / EQUILIBRIUM_HALF_LIFE_RATIO and all its possible parents are kept: a long-lived daughter
    does not reach equilibrium in the patient, nor do its own daughters, even short-lived ones (Th227 ->
    Ra223 (11.4 d) -> Rn219 (4 s)). The decay tree is rebuilt from the (Z, A) change of each decay mode
    (the chain file lists the nuclides in decay order, with cumulative yields).
    """
    daughters = chain["daughters"]
    t_parent = chain.get("half_life_s") or daughters[0].get("half_life_s")
    kept_names, kept, excluded = set(), [], []
    for i, d in enumerate(daughters):
        if d.get("is_stable", False):
            continue
        if i == 0:
            kept_names.add(d["nuclide"])
            kept.append(d)
            continue
        parents = [
            e["nuclide"]
            for e in daughters[:i]
            if e["nuclide"] != d["nuclide"]
            and not e.get("is_stable", False)
            and (d["z"] - e["z"], d["a"] - e["a"]) in _DECAY_STEPS
        ]
        t = d.get("half_life_s")
        short = t is not None and t_parent is not None and t < t_parent / EQUILIBRIUM_HALF_LIFE_RATIO
        if short and parents and all(p in kept_names for p in parents):
            kept_names.add(d["nuclide"])
            kept.append(d)
        else:
            excluded.append(d)
    return kept, excluded


def get_cascade_energy_constant_mev(
    isotope: str, cascade: bool = True
) -> tuple[float, dict[str, float]]:
    """
    Calculate non-penetrating decay energy constant (Delta_local) in MeV.

    If cascade=True and the isotope has a curated decay chain (e.g. Ac225, Ra223, Pb212),
    computes the secular equilibrium weighted sum over the progeny in equilibrium with the parent:
        Delta_local,cascade = sum(Y_i * Delta_local(i))
    where Y_i is the cumulative branching yield to daughter i. Daughters that are not short-lived
    compared with the parent, and their own daughters, are excluded (cascade_members), with a WARNING.

    Returns:
        (total_delta_mev, contributions_dict)
    """
    chain = get_decay_chain(isotope) if cascade else None
    if not chain or "daughters" not in chain:
        single_delta = get_non_penetrating_energy_constant_mev(isotope)
        return single_delta, {isotope: single_delta}

    kept, excluded = cascade_members(chain)
    if excluded:
        logger.warning(
            f"Cascade of {isotope}: not in secular equilibrium with {isotope} (half-life not below "
            f"1/{EQUILIBRIUM_HALF_LIFE_RATIO:g} of {chain.get('half_life_str', 'the parent')}), excluded: "
            + ", ".join(f"{d['nuclide']} ({d.get('half_life_str', '?')}, yield {float(d.get('yield', 1)):.3g})" for d in excluded)
        )
    total_delta = 0.0
    contributions = {}
    for daughter in kept:
        nuclide = daughter["nuclide"]
        branch_yield = float(daughter.get("yield", 1.0))
        try:
            delta_i = get_non_penetrating_energy_constant_mev(nuclide)
        except ValueError as e:
            raise ValueError(
                f"Cascade of {isotope}: no ICRP 107 data for the daughter {nuclide} "
                f"(cumulative yield {branch_yield:g}), the cascade energy would miss it. "
                f"Without the cascade (--no-cascade), only {isotope} is counted."
            ) from e
        contrib = branch_yield * delta_i
        total_delta += contrib
        contributions[nuclide] = contrib

    return total_delta, contributions
