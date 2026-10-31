#!/usr/bin/env python3
"""
Fetch OpenDose S-values and organ masses from the OpenDose website, in the radot format
(<out>/<phantom>/<radionuclide>.json.gz, <out>/<phantom>/masses.json), to build the radot
OpenDose data archive (Zenodo).

The website has no export: one page per (phantom, source organ, radionuclide), about 6 s each on
the server side, so the pages are fetched in parallel (--workers, keep it small: the server is
shared). Resumable: the organs already fetched are skipped, so the script can be stopped and
restarted. Standard library only (curl is used when Python cannot verify the server certificate).

    python scripts/fetch_opendose.py                          # default list, 4 workers, data/opendose_fetch/
    python scripts/fetch_opendose.py -r Cu67 Tb161 -w 2
    python scripts/fetch_opendose.py --check                  # completeness of the fetched files only

Licence of the data: CC BY-SA 4.0 (OpenDose collaboration), research and education only.
"""

import argparse
import gzip
import json
import shutil
import ssl
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from html.parser import HTMLParser
from pathlib import Path

BASE_URL = "https://opendose.polsl.pl/svalues/"
PHANTOMS = {"icrp_110_male": "ICRP_110_AM", "icrp_110_female": "ICRP_110_AF"}

# Content of radot_opendose v1.0 (the first 11 were extracted in June 2026)
RADIONUCLIDES = [
    "Ac-225", "At-211", "Ho-166", "I-125", "I-131", "Lu-177", "Pb-212", "Ra-223", "Sm-153",
    "Tb-161", "Y-90",
    # beta therapy
    "Cu-67", "Er-169", "Re-186", "Re-188", "Sc-47", "Sr-89", "P-32", "Ag-111",
    # alpha therapy
    "Bi-213", "Th-227", "Tb-149", "Bi-212",
    # imaging and theranostic pairs
    "F-18", "Ga-68", "Cu-64", "Zr-89", "In-111", "Tc-99m", "I-123", "I-124", "Y-86", "Sc-44",
    "Tb-152", "Tb-155", "Pb-203",
]


def radot_name(name: str) -> str:
    """Lu-177 -> lu177, Tc-99m -> tc99m (radot file names)."""
    return name.replace("-", "").replace(" ", "").lower()


_use_curl = False
_lock = threading.Lock()


def get(url: str, timeout: int = 120, retries: int = 4) -> str:
    """GET a page; falls back to curl if Python cannot verify the certificate."""
    global _use_curl
    for attempt in range(retries):
        try:
            if not _use_curl:
                req = urllib.request.Request(url, headers={"User-Agent": "radot-fetch-opendose"})
                with urllib.request.urlopen(req, timeout=timeout) as r:
                    return r.read().decode("utf-8", errors="replace")
            out = subprocess.run(
                ["curl", "-sSfL", "--max-time", str(timeout), url],
                capture_output=True, check=True,
            )
            return out.stdout.decode("utf-8", errors="replace")
        except urllib.error.URLError as e:
            if isinstance(e.reason, ssl.SSLError) and shutil.which("curl"):
                with _lock:
                    if not _use_curl:
                        print("Python cannot verify the certificate: using curl", flush=True)
                    _use_curl = True
                continue
            err = e
        except Exception as e:  # timeouts, curl errors
            err = e
        time.sleep(5 * (attempt + 1))
    raise RuntimeError(f"{url}: {err}")


class SelectParser(HTMLParser):
    """Options of the <select> whose id or name contains a keyword."""

    def __init__(self, keyword):
        super().__init__()
        self.keyword, self.inside, self.value, self.options = keyword, False, None, {}

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == "select" and (self.keyword in a.get("id", "").lower() or self.keyword in a.get("name", "").lower()):
            self.inside = True
        elif tag == "option" and self.inside:
            self.value = a.get("value")

    def handle_endtag(self, tag):
        if tag == "select":
            self.inside = False

    def handle_data(self, data):
        text = data.strip()
        if self.inside and self.value is not None and text and not text.startswith("Select"):
            self.options[self.value] = text
            self.value = None


class TableParser(HTMLParser):
    """Rows of the S-value table: target organ, S-value, ..., mass."""

    def __init__(self):
        super().__init__()
        self.depth, self.in_body, self.row, self.cell, self.rows = 0, False, None, None, []

    def handle_starttag(self, tag, attrs):
        if tag == "tbody":
            self.in_body = True
        elif tag == "tr" and self.in_body:
            self.row = []
        elif tag in ("td", "th") and self.row is not None:
            self.cell = []

    def handle_endtag(self, tag):
        if tag == "tbody":
            self.in_body = False
        elif tag in ("td", "th") and self.cell is not None:
            self.row.append("".join(self.cell).strip())
            self.cell = None
        elif tag == "tr" and self.row is not None:
            if self.row:
                self.rows.append(self.row)
            self.row = None

    def handle_data(self, data):
        if self.cell is not None:
            self.cell.append(data)


def fetch_source(phantom_url, source, isotope_label):
    """S-values {target: S} and masses {target: g} of one source organ (empty if none)."""
    url = (f"{BASE_URL}{urllib.parse.quote(phantom_url)}/"
           f"{urllib.parse.quote(source).replace('/', '%252F')}/{urllib.parse.quote(isotope_label)}/")
    t = TableParser()
    t.feed(get(url))
    s_values, masses = {}, {}
    for row in t.rows:
        if len(row) >= 4:
            try:
                s_values[row[0]] = float(row[1])
                masses[row[0]] = float(row[3])
            except ValueError:
                continue
    return s_values, masses


def load_json(path: Path) -> dict:
    if not path.exists():
        return {}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as f:
        return json.load(f)


def save_json(path: Path, data: dict):
    tmp = path.with_name(path.name + ".tmp")
    if path.suffix == ".gz":
        with open(tmp, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as g:
            g.write(json.dumps(data, indent=4).encode())
    else:
        tmp.write_text(json.dumps(data, indent=4))
    tmp.replace(path)


def fetch(radionuclides, out: Path, workers: int):
    html = get(BASE_URL)
    p = SelectParser("isotope")
    p.feed(html)
    # lu177 -> "Lu-177 (t½=6.647d)" (exact name: Bi-212 and Bi-212n are different)
    labels = {radot_name(lab.split(" (")[0]): lab for lab in p.options.values()}
    p = SelectParser("source")
    p.feed(html)
    sources = sorted(set(p.options.values()))
    radionuclides = [radot_name(r) for r in radionuclides]
    missing = [r for r in radionuclides if r not in labels]
    if missing:
        sys.exit(f"Not on OpenDose: {missing}")
    print(f"{len(sources)} source organs on the website, {len(radionuclides)} radionuclides, {workers} workers", flush=True)

    for phantom, phantom_url in PHANTOMS.items():
        folder = out / phantom
        folder.mkdir(parents=True, exist_ok=True)
        masses_file = folder / "masses.json"
        masses = load_json(masses_file)
        empty_file = folder / "empty_sources.json"  # sources without data on this phantom
        empty = load_json(empty_file)
        for name in radionuclides:
            target = folder / f"{name}.json.gz"
            matrix = load_json(target)
            done = set(matrix) | set(empty.get(name, []))
            todo = [s for s in sources if s not in done]
            if not todo:
                print(f"{phantom} {name}: complete ({len(matrix)} sources)", flush=True)
                continue
            t0 = time.time()
            print(f"{phantom} {name}: {len(todo)} source organs to fetch", flush=True)
            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = {pool.submit(fetch_source, phantom_url, s, labels[name]): s for s in todo}
                for i, fut in enumerate(as_completed(futures), 1):
                    src = futures[fut]
                    try:
                        s_values, m = fut.result()
                    except Exception as e:
                        print(f"  {src}: FAILED ({e}), rerun to retry", flush=True)
                        continue
                    if s_values:
                        matrix[src] = s_values
                        masses.update(m)
                    else:
                        empty.setdefault(name, []).append(src)
                    if i % 20 == 0 or i == len(todo):
                        save_json(target, matrix)
                        save_json(masses_file, masses)
                        save_json(empty_file, empty)
                        rate = (time.time() - t0) / i
                        print(f"  {i}/{len(todo)} ({rate:.1f} s/organ, {rate * (len(todo) - i) / 60:.0f} min left)", flush=True)
            save_json(target, matrix)
            save_json(masses_file, masses)
            save_json(empty_file, empty)


def check(out: Path) -> bool:
    """Every S-value file has the same source and target organs as masses.json, no gap."""
    ok = True
    for phantom in PHANTOMS:
        folder = out / phantom
        masses = load_json(folder / "masses.json")
        files = sorted(folder.glob("*.json.gz"))
        print(f"== {phantom}: {len(masses)} organs in masses.json, {len(files)} radionuclides")
        for f in files:
            m = load_json(f)
            targets = set().union(*[set(v) for v in m.values()]) if m else set()
            bad = len(m) != len(masses) or targets != set(masses) or any(len(v) != len(masses) for v in m.values())
            ok &= not bad
            print(f"  {'BAD' if bad else 'ok '} {f.name:16s} {len(m)} sources, {len(targets)} targets")
    return ok


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-r", "--radionuclides", nargs="+", default=RADIONUCLIDES, help="e.g. Lu-177 Y-90 (default: the radot list)")
    parser.add_argument("-o", "--output", default=None, help="output folder (default: <repo>/data/opendose_fetch)")
    parser.add_argument("-w", "--workers", type=int, default=4, help="parallel requests (default 4)")
    parser.add_argument("--check", action="store_true", help="only check the fetched files")
    args = parser.parse_args()
    out = Path(args.output) if args.output else Path(__file__).resolve().parents[1] / "data" / "opendose_fetch"
    if not args.check:
        fetch(args.radionuclides, out, args.workers)
    sys.exit(0 if check(out) else 1)


if __name__ == "__main__":
    main()
