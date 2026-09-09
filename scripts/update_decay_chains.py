#!/usr/bin/env python3
"""
Maintainer script to retrieve authoritative radioactive decay chain data
from the IAEA Nuclear Data Services (LiveChart of Nuclides REST API)
and generate `src/radot/data/decay_chains.json`.

Authoritative Source:
  IAEA Nuclear Data Section (NDS) - LiveChart of Nuclides
  https://nds.iaea.org/relnsd/v0/data?fields=ground_states&nuclides=...
"""

import csv
import datetime
import io
import json
from pathlib import Path
import re
import sys
import urllib.request

ELEMENTS = [
    "",
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
    "Rf",
    "Db",
    "Sg",
    "Bh",
    "Hs",
    "Mt",
    "Ds",
]
ELEM_MAP = {sym.lower(): i for i, sym in enumerate(ELEMENTS) if sym}

# Curated list of medical radionuclides with cascades or daughter ingrowth
DEFAULT_RADIONUCLIDES = [
    # Alpha therapy cascades
    "Ac225",
    "Ra223",
    "Pb212",
    "Bi212",
    "At211",
    "Th227",
    "Tb149",
    # Generator / beta cascades
    "Sr90",
    "Ge68",
    "Sr82",
    "Mo99",
    "W188",
]


def parse_nuclide(name: str) -> tuple[int, int]:
    """Parse nuclide name (e.g. 'Ac225', 'Ac-225', '225Ac') into (Z, A)."""
    clean = name.strip()
    m = re.match(r"^([A-Za-z]+)[-_]?(\d+)$", clean)
    if m:
        sym, a = m.group(1), int(m.group(2))
        return ELEM_MAP[sym.lower()], a
    m = re.match(r"^(\d+)[-_]?([A-Za-z]+)$", clean)
    if m:
        a, sym = int(m.group(1)), m.group(2)
        return ELEM_MAP[sym.lower()], a
    raise ValueError(f"Cannot parse radionuclide name: {name}")


def fetch_nuclide_info(z: int, a: int) -> dict:
    """Fetch nuclide properties from IAEA LiveChart API."""
    sym = ELEMENTS[z].lower()
    nuclide_query = f"{a}{sym}"
    url = f"https://nds.iaea.org/relnsd/v0/data?fields=ground_states&nuclides={nuclide_query}"
    req = urllib.request.Request(url, headers={"User-Agent": "radot-decay-chains/1.0"})
    with urllib.request.urlopen(req, timeout=15) as resp:
        text = resp.read().decode("utf-8")
        reader = csv.DictReader(io.StringIO(text))
        rows = list(reader)
        return rows[0] if rows else {}


def get_daughter_z_a(z: int, a: int, decay_mode: str) -> tuple[int, int] | None:
    """Determine daughter (Z, A) given parent (Z, A) and decay mode."""
    m = decay_mode.upper().strip()
    if m in ["A", "ALPHA"]:
        return z - 2, a - 4
    elif m in ["B-", "BMINUS", "B", "2B-"]:
        return z + 1, a
    elif m in ["B+", "EC", "EC+B+", "B+EC", "2EC"]:
        return z - 1, a
    elif m in ["IT"]:
        return z, a
    elif m in ["P"]:
        return z - 1, a - 1
    elif m in ["N"]:
        return z, a - 1
    return None


def trace_decay_chain(
    start_z: int, start_a: int, min_branch_pct: float = 0.05
) -> list[dict]:
    """
    Recursively trace decay chain starting from (start_z, start_a) down to stable nuclides.
    Returns list of progeny in order of decay with cumulative yields in secular equilibrium.
    """
    queue = [(start_z, start_a, 1.0)]  # (z, a, cumulative_yield)
    nuclide_yields: dict[str, dict] = {}
    order = []

    while queue:
        z, a, cum_yield = queue.pop(0)
        sym = ELEMENTS[z]
        nuc_id = f"{sym}{a}"

        # Fetch IAEA properties
        info = fetch_nuclide_info(z, a)
        if not info:
            continue

        hl_str = info.get("half_life", "").strip()
        unit_hl = info.get("unit_hl", "").strip()
        hl_sec_str = info.get("half_life_sec", "").strip()
        hl_sec = float(hl_sec_str) if hl_sec_str else None

        is_stable = (hl_str.upper() == "STABLE") or (
            hl_sec is not None and hl_sec > 1e16
        )

        if nuc_id not in nuclide_yields:
            nuclide_yields[nuc_id] = {
                "nuclide": nuc_id,
                "element": sym,
                "a": a,
                "z": z,
                "half_life_str": f"{hl_str} {unit_hl}".strip(),
                "half_life_s": hl_sec,
                "is_stable": is_stable,
                "cumulative_yield": cum_yield,
            }
            order.append(nuc_id)
        else:
            nuclide_yields[nuc_id]["cumulative_yield"] += cum_yield

        # Stop tracing if stable
        if is_stable:
            continue

        # Trace branches (decay_1, decay_2, decay_3)
        for i in [1, 2, 3]:
            d_mode = info.get(f"decay_{i}", "").strip()
            pct_str = info.get(f"decay_{i}_%", "").strip()
            if d_mode and pct_str:
                try:
                    pct = float(pct_str)
                except ValueError:
                    continue
                if pct < min_branch_pct:
                    continue
                branch = pct / 100.0
                d_coord = get_daughter_z_a(z, a, d_mode)
                if d_coord and d_coord != (z, a):
                    dz, da = d_coord
                    queue.append((dz, da, cum_yield * branch))

    return [nuclide_yields[k] for k in order]


def main():
    repo_root = Path(__file__).resolve().parent.parent
    target_json = repo_root / "src" / "radot" / "data" / "decay_chains.json"
    target_json.parent.mkdir(parents=True, exist_ok=True)

    radionuclides = sys.argv[1:] if len(sys.argv) > 1 else DEFAULT_RADIONUCLIDES

    print("=" * 70)
    print("RADOT Decay Chains Generator (IAEA LiveChart REST API)")
    print("=" * 70)

    chains_data = {}
    for rad in radionuclides:
        try:
            z, a = parse_nuclide(rad)
            sym = ELEMENTS[z]
            std_name = f"{sym}{a}"
            print(f"\nProcessing {std_name} (Z={z}, A={a})...")
            chain = trace_decay_chain(z, a)

            # Filter daughters (all radioactive progeny + terminal stable)
            daughters = []
            for item in chain:
                daughters.append(
                    {
                        "nuclide": item["nuclide"],
                        "z": item["z"],
                        "a": item["a"],
                        "yield": round(item["cumulative_yield"], 6),
                        "half_life_str": item["half_life_str"],
                        "half_life_s": item["half_life_s"],
                        "is_stable": item["is_stable"],
                    }
                )
                print(
                    f"  -> {item['nuclide']:7s}: yield={item['cumulative_yield']:.4f}, "
                    f"T1/2={item['half_life_str']}"
                )

            chains_data[std_name] = {
                "parent": std_name,
                "z": z,
                "a": a,
                "half_life_str": chain[0]["half_life_str"] if chain else "",
                "half_life_s": chain[0]["half_life_s"] if chain else None,
                "daughters": daughters,
            }
        except Exception as e:
            print(f"[-] Error processing {rad}: {e}", file=sys.stderr)

    output = {
        "_metadata": {
            "source": "IAEA Nuclear Data Services (LiveChart of Nuclides REST API)",
            "url": "https://nds.iaea.org/relnsd/v0/data",
            "generated_at_utc": datetime.datetime.now(
                datetime.timezone.utc
            ).isoformat(),
            "description": "Authoritative radioactive decay chains and secular equilibrium yields for medical radionuclides",
        },
        "chains": chains_data,
    }

    with open(target_json, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2)

    print("\n" + "=" * 70)
    print(f"Successfully saved {len(chains_data)} decay chains to: {target_json}")
    print("=" * 70)


if __name__ == "__main__":
    main()
