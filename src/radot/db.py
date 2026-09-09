import json
import re
import gzip
import importlib.resources
import logging

logger = logging.getLogger("radot")
import icrp107_database
from typing import Dict, Any, Optional
from pathlib import Path


def get_physical_half_life_h(isotope: str) -> float:
    """
    Look up the physical half-life (Tp) of an isotope in hours.
    Loads data from the required 'icrp107-database' package.
    """
    iso_clean = re.sub(r"[-_\s]", "", isotope).lower()

    # Locate the icrp107 directory inside the installed package
    icrp_dir = importlib.resources.files("icrp107_database").joinpath("icrp107")

    matched_path = None
    for p in icrp_dir.iterdir():
        if p.name.endswith(".json"):
            name_without_ext = p.name[:-5]
            if re.sub(r"[-_\s]", "", name_without_ext).lower() == iso_clean:
                matched_path = p
                break

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
    """Normalize isotope name to standard filename format (e.g. lu177, y90)."""
    clean = re.sub(r"[-_\s]", "", isotope).lower()
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
        try:
            norm_phantom = normalize_phantom(phantom)
            masses = _load_masses(norm_phantom)
            for official_name in masses.keys():
                if (
                    official_name.lower() == resolved.lower()
                    or official_name.lower() == standardized.lower()
                ):
                    return official_name
        except Exception:
            pass

    # Fallback or resolved name
    if resolved == user_name:
        return standardized.capitalize()
    return resolved


DEFAULT_DATABASE_URL = "https://www.creatis.insa-lyon.fr/~dsarrut/radot/opendose/"


def _ensure_file_downloaded(phantom: str, filename: str) -> Path:
    """Ensure the requested file is present locally, downloading it on-demand if missing."""
    p_dir = (
        importlib.resources.files("radot")
        .joinpath("data")
        .joinpath("opendose")
        .joinpath(phantom)
    )
    dest_dir = Path(str(p_dir))
    dest_dir.mkdir(parents=True, exist_ok=True)

    local_file = dest_dir / filename
    if not local_file.exists():
        url = f"{DEFAULT_DATABASE_URL}{phantom}/{filename}"
        logger.info(
            f"On-demand download of {filename} for {phantom} from web registry..."
        )
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
            with urllib.request.urlopen(req, timeout=20) as response:
                content = response.read()
            with open(local_file, "wb") as f:
                f.write(content)
        except Exception as e:
            logger.warning(
                f"Download of {filename} from central registry failed: {e}. Trying to download from OpenDose directly..."
            )
            try:
                if filename == "masses.json":
                    download_opendose_data(phantom, isotope="lu177", force=True)
                elif filename.endswith(".json.gz"):
                    iso_name = filename[:-8]
                    download_opendose_data(phantom, isotope=iso_name, force=True)
                else:
                    raise ValueError(f"Unknown database file name format: {filename}")

                if not local_file.exists():
                    raise FileNotFoundError(
                        f"OpenDose download did not produce the expected file {filename}"
                    )
            except Exception as od_err:
                raise FileNotFoundError(
                    f"RADOT Database: Failed to download '{filename}' for phantom '{phantom}'.\n"
                    f"  - First download (Creatis) error: {e}\n"
                    f"  - Second download (OpenDose) error: {od_err}"
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
            import importlib.resources
            from pathlib import Path

            p_dir = (
                importlib.resources.files("radot")
                .joinpath("data")
                .joinpath("opendose")
                .joinpath(phantom)
            )
            local_json = Path(str(p_dir)) / f"{isotope}.json"

            if local_json.exists():
                with open(local_json, "r", encoding="utf-8") as f:
                    _s_values_cache[key] = json.load(f)
            else:
                s_file_gz = _ensure_file_downloaded(phantom, f"{isotope}.json.gz")
                try:
                    with gzip.open(s_file_gz, "rt", encoding="utf-8") as f:
                        _s_values_cache[key] = json.load(f)
                except Exception:
                    # In case they decompressed the .json.gz file in-place keeping the .gz extension
                    with open(s_file_gz, "r", encoding="utf-8") as f:
                        _s_values_cache[key] = json.load(f)
        except Exception as e:
            raise FileNotFoundError(
                f"RADOT Database: Failed S-value matrix file for '{isotope}' on '{phantom}': {e}"
            )
    return _s_values_cache[key]


def _find_matching_key(d: dict, key: str) -> Optional[str]:
    """Helper to find a key in a dictionary with flexible matches."""
    if key in d:
        return key
    for k in d.keys():
        if k.lower() == key.lower():
            return k

    def normalize_key(name: str) -> str:
        name = name.lower()
        name = re.sub(r"\s+c\+m\+p$", "", name)
        name = re.sub(r"\s+lobe$", "", name)
        name = re.sub(r"\s+wall$", "", name)
        name = re.sub(r"\s+contents$", "", name)
        return re.sub(r"[-_\s]+", "", name)

    norm_key = normalize_key(key)
    for k in d.keys():
        if normalize_key(k) == norm_key:
            return k

    for k in d.keys():
        if k.lower().startswith(key.lower()) or key.lower().startswith(k.lower()):
            return k

    return None


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
        )

    tgt_matrix = s_matrix[src_key]

    tgt_key = _find_matching_key(tgt_matrix, resolved_tgt)
    if tgt_key is not None:
        return tgt_matrix[tgt_key]

    raise ValueError(
        f"RADOT Database: Target organ '{target_organ}' (resolved as '{resolved_tgt}') not found in S-value matrix for source '{resolved_src}'."
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
    Downloads S-values and target organ reference masses from OpenDose registry.
    Saves S-value matrix as <isotope>.json and masses as masses.json.
    """
    norm_phantom = normalize_phantom(phantom)

    if norm_phantom == "icrp_110_male":
        phantom_url_name = "ICRP_110_AM"
    else:
        phantom_url_name = "ICRP_110_AF"

    isotopes_to_download = []
    if isotope:
        isotopes_to_download = [isotope]

    dest_dir = Path(
        str(
            importlib.resources.files("radot")
            .joinpath("data")
            .joinpath("opendose")
            .joinpath(norm_phantom)
        )
    )
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
            except Exception:
                pass

        masses_file = dest_dir / "masses.json"
        if masses_file.exists() and not force:
            try:
                with open(masses_file, "r", encoding="utf-8") as f:
                    masses = json.load(f)
            except Exception:
                pass

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

    iso_clean = re.sub(r"[-_\s]", "", isotope).lower()
    icrp_dir = importlib.resources.files("icrp107_database").joinpath("icrp107")

    matched_path = None
    for p in icrp_dir.iterdir():
        if p.name.endswith(".json"):
            name_without_ext = p.name[:-5]
            if re.sub(r"[-_\s]", "", name_without_ext).lower() == iso_clean:
                matched_path = p
                break

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
        except Exception:
            _decay_chains_cache = {"chains": {}}
    return _decay_chains_cache


def get_decay_chain(isotope: str) -> Optional[dict]:
    """
    Retrieve decay chain metadata and daughters for an isotope if available.
    Matches case-insensitively with flexible formatting (e.g. 'Ac225', 'Ac-225', '225Ac').
    """
    import re

    clean = re.sub(r"[-_\s]", "", isotope).lower()
    data = _load_decay_chains()
    chains = data.get("chains", {})
    for parent_key, chain_info in chains.items():
        if re.sub(r"[-_\s]", "", parent_key).lower() == clean:
            return chain_info
    return None


def get_cascade_energy_constant_mev(
    isotope: str, cascade: bool = True
) -> tuple[float, dict[str, float]]:
    """
    Calculate non-penetrating decay energy constant (Delta_local) in MeV.

    If cascade=True and the isotope has a curated decay chain (e.g. Ac225, Ra223, Pb212),
    computes the secular equilibrium weighted sum over all progeny in the cascade:
        Delta_local,cascade = sum(Y_i * Delta_local(i))
    where Y_i is the cumulative branching yield to daughter i.

    Returns:
        (total_delta_mev, contributions_dict)
    """
    chain = get_decay_chain(isotope) if cascade else None
    if not chain or "daughters" not in chain:
        single_delta = get_non_penetrating_energy_constant_mev(isotope)
        return single_delta, {isotope: single_delta}

    total_delta = 0.0
    contributions = {}
    for daughter in chain["daughters"]:
        if daughter.get("is_stable", False):
            continue
        nuclide = daughter["nuclide"]
        branch_yield = float(daughter.get("yield", 1.0))
        try:
            delta_i = get_non_penetrating_energy_constant_mev(nuclide)
            contrib = branch_yield * delta_i
            total_delta += contrib
            contributions[nuclide] = contrib
        except Exception:
            pass

    return total_delta, contributions
