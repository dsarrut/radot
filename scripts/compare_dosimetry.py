#!/usr/bin/env python3
"""
Compare dosimetry results across 5 methods for p007:
1. Jackson 2013 (MTP reference)
2. Madsen 2018 (STP with Jackson Teff)
3. Hänscheid 2017 (STP with Jackson Teff)
4. Hänscheid 2018 (STP default Teff = 2*t)
5. LDM (STP TP2 with physical decay)

Outputs a Markdown comparison table with absolute doses and relative % differences wrt Jackson 2013.
"""

import argparse
import json
from pathlib import Path
import sys


METHODS = [
    {
        "id": "jackson2013",
        "name": "Jackson 2013 (MTP)",
        "config": "dose_jackson2013.json",
        "results": "dose_jackson2013_results.json",
        "is_ref": True,
    },
    {
        "id": "madsen2018",
        "name": "Madsen 2018 (Teff)",
        "config": "dose_madsen2018.json",
        "results": "dose_madsen2018_results.json",
        "is_ref": False,
    },
    {
        "id": "hanscheid2017",
        "name": "Hänscheid 2017 (Teff)",
        "config": "dose_hanscheid2017.json",
        "results": "dose_hanscheid2017_results.json",
        "is_ref": False,
    },
    {
        "id": "hanscheid2018",
        "name": "Hänscheid 2018 (2*t)",
        "config": "dose_hanscheid2018.json",
        "results": "dose_hanscheid2018_results.json",
        "is_ref": False,
    },
    {
        "id": "ldm",
        "name": "LDM (TP2)",
        "config": "dose_ldm.json",
        "results": "dose_ldm_results.json",
        "is_ref": False,
    },
]


def load_or_compute_results(directory: Path, method_info: dict) -> dict:
    results_path = directory / method_info["results"]
    if results_path.exists():
        with open(results_path, "r") as f:
            return json.load(f)

    config_path = directory / method_info["config"]
    if config_path.exists():
        try:
            from radot.dose import compute_roi_dose_file
            print(f"[INFO] Computing {method_info['name']} from {config_path.name}...", file=sys.stderr)
            res = compute_roi_dose_file(config_path)
            with open(results_path, "w") as f:
                json.dump(res, f, indent=4)
            return res
        except Exception as e:
            print(f"[WARNING] Failed to compute {config_path.name}: {e}", file=sys.stderr)
            return {}

    # Check fallback for jackson: dose_p007.json
    if method_info["id"] == "jackson2013":
        fallback = directory / "dose_p007.json"
        if fallback.exists():
            with open(fallback, "r") as f:
                return json.load(f)

    print(f"[WARNING] Results or config not found for {method_info['name']} in {directory}", file=sys.stderr)
    return {}


def generate_markdown_table(directory: Path) -> str:
    data = {}
    all_rois = set()

    for m in METHODS:
        res = load_or_compute_results(directory, m)
        data[m["id"]] = res
        for roi in res.keys():
            all_rois.add(roi)

    sorted_rois = sorted(list(all_rois))
    ref_id = "jackson2013"
    ref_data = data.get(ref_id, {})

    headers = [
        "ROI",
        "Jackson 2013 (MTP)",
        "Madsen 2018 (Teff)",
        "Hänscheid 2017 (Teff)",
        "Hänscheid 2018 (2*t)",
        "LDM (TP2)",
    ]

    lines = []
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("|" + "|".join(["---"] * len(headers)) + "|")

    for roi in sorted_rois:
        row = [f"**{roi}**"]
        ref_dose = ref_data.get(roi, {}).get("mean_dose_gy")

        # Ref column
        if ref_dose is not None:
            row.append(f"**{ref_dose:.3f} Gy**")
        else:
            row.append("N/A")

        # Other methods
        for m in METHODS:
            if m["is_ref"]:
                continue
            m_res = data.get(m["id"], {}).get(roi, {})
            dose = m_res.get("mean_dose_gy")
            if dose is None:
                row.append("N/A")
            elif ref_dose is None or ref_dose == 0:
                row.append(f"{dose:.3f} Gy")
            else:
                diff_rel = ((dose - ref_dose) / ref_dose) * 100.0
                sign = "+" if diff_rel > 0 else ""
                row.append(f"{dose:.3f} Gy ({sign}{diff_rel:.1f}%)")

        lines.append("| " + " | ".join(row) + " |")

    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Compare 5 dosimetry methods wrt Jackson 2013.")
    parser.add_argument(
        "directory",
        nargs="?",
        default=".",
        help="Directory containing the patient dose configs and results (default: current directory).",
    )
    parser.add_argument(
        "-o", "--output",
        help="Optional output markdown file.",
    )
    args = parser.parse_args()

    target_dir = Path(args.directory).resolve()
    table_md = generate_markdown_table(target_dir)

    print("\n### Absorbed Dose Comparison on " + target_dir.name + "\n")
    print(table_md)
    print()

    if args.output:
        out_path = Path(args.output).resolve()
        with open(out_path, "w") as f:
            f.write(table_md + "\n")
        print(f"Table saved to: {out_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
