#!/usr/bin/env python3
"""
Generate docs/8_validation_snmmi.md: radot doses compared with the SNMMI 177Lu Dosimetry
Challenge participants.

Reads the stored radot results (results/t5_radot.json, written by t5.py) and the participants'
statistics (challenge_stats.json): the challenge data are not needed to regenerate the page.

Usage:
    python validation/snmmi/report.py          # write docs/8_validation_snmmi.md
    python validation/snmmi/report.py --check  # exit 1 if the page is out of date
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cases import PATIENTS, RESULTS_DIR  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "docs" / "8_validation_snmmi.md"
T5_RESULTS = RESULTS_DIR / "t5_radot.json"
REPO_URL = "https://github.com/dsarrut/radot/tree/main/validation/snmmi"
DATA_URL = "https://deepblue.lib.umich.edu/data/collections/hm50ts030"


def fmt_dose(r, method):
    if r[method] is None:
        return "—"
    return f"{r[method]:.2f}"


def fmt_diff(r, method):
    if r.get(method) is None:
        return "—"
    flag = "" if r[f"{method}_in_iqr"] else " ✗"
    return f"{r[f'{method}_diff_pct']:+.1f} %{flag}"


def fmt_unc(r):
    return "—" if r.get("gate_unc_pct") is None else f"{r['gate_unc_pct']:.1f} %"


def methods_in(t5):
    """svalue and ldm, and gate when the Monte Carlo results are present."""
    has_gate = any(r.get("gate") is not None for results in t5.values() for r in results.values())
    return ("svalue", "ldm", "gate") if has_gate else ("svalue", "ldm")


def checked(name, voi):
    """VOIs compared without ambiguity (the liver of A depends on the liver definition)."""
    return not (voi == "Liver (no lesions)" or (voi == "Liver" and PATIENTS[name].challenge_id == "A"))


def mass_section(t5, order):
    """Dose with the CT mass, and with the VOI volume x 1.00 or 1.06 g/mL (the dose scales as 1 / mass)."""
    densities = (1.00, 1.06)

    def rows():
        for name, results in t5.items():
            for voi in sorted(results, key=lambda v: (order.index(v) if v in order else len(order), v)):
                r = results[voi]
                method = "svalue" if r["svalue"] is not None else "ldm"
                median = r["challenge"]["median"]
                diffs = [100.0 * (r[method] / median - 1.0)]
                for rho in densities:
                    dose = r[method] * r["mass_g"] / (rho * r["volume_ml"])
                    diffs.append(100.0 * (dose / median - 1.0))
                yield name, voi, r, method, diffs

    # VOIs whose definition is not ambiguous (the liver VOIs contain lesions)
    not_liver = [d for name, voi, _, _, d in rows() if not voi.startswith("Liver")]
    lo = [min(d[i] for d in not_liver) for i in range(3)]
    hi = [max(d[i] for d in not_liver) for i in range(3)]
    densities_ct = [r["mass_g"] / r["volume_ml"] for name, voi, r, _, _ in rows() if not voi.startswith("Liver")]

    lines = [
        "",
        "## Effect of the mass",
        "",
        "radot computes the VOI mass from the CT (density = 1 + HU/1000, CT resampled to the TIA grid). "
        "Many dosimetry tools use the VOI volume times a fixed density instead, often 1.00 or 1.06 g/mL. "
        "The participants' masses are not published. All the methods divide by the mass, so the dose "
        "scales as 1 / mass: below, the radot dose (`svalue`, or `ldm` for the lesions) recomputed with "
        "the other mass conventions, as the difference with the participants' median.",
        "",
        "| Patient | VOI | Volume (mL) | CT density (g/mL) | Diff. CT mass | Diff. volume × 1.00 | Diff. volume × 1.06 |",
        "|---|---|--:|--:|--:|--:|--:|",
    ]
    for name, voi, r, method, d in rows():
        lines.append(
            f"| {PATIENTS[name].challenge_id} | {voi} | {r['volume_ml']:.0f} | {r['mass_g'] / r['volume_ml']:.3f} "
            f"| {d[0]:+.1f} % | {d[1]:+.1f} % | {d[2]:+.1f} % |"
        )
    lines += [
        "",
        f"- The CT density of the VOIs other than the liver is {min(densities_ct):.3f}–{max(densities_ct):.3f} g/mL "
        "(kidneys about 1.015: low HU and fat at the edges). Going from 1.00 to 1.06 g/mL changes the "
        "dose by 6 %, as much as the participants' spread (QCD 3–7 %).",
        f"- For the kidneys, spleen and lesions, the difference with the median is {lo[0]:+.1f} to {hi[0]:+.1f} % "
        f"with the CT mass, {lo[1]:+.1f} to {hi[1]:+.1f} % with volume × 1.00, and {lo[2]:+.1f} to "
        f"{hi[2]:+.1f} % with volume × 1.06: with a density of 1 g/mL, radot is close to the median.",
        "- The liver does not follow (B: closer to the median with 1.06 g/mL, or with 1.00 g/mL without the "
        "lesions): the median likely mixes mass conventions and liver definitions.",
        "- The CT computation itself is not a source of difference: resampling the CT with or without "
        "smoothing, or computing the mass on the original CT grid, changes the organ masses by less than 0.5 %.",
        "- radot keeps the CT mass (patient-specific). The difference with the median is a convention, "
        "not an error of the dose calculation.",
    ]
    return lines


def render(t5):
    methods = methods_in(t5)
    has_gate = "gate" in methods
    diffs = [
        r[f"{m}_diff_pct"]
        for name, results in t5.items()
        for voi, r in results.items()
        if checked(name, voi)
        for m in methods
        if r.get(m) is not None
    ]
    outside = [
        f"{PATIENTS[name].challenge_id} {voi} ({m})"
        for name, results in t5.items()
        for voi, r in results.items()
        if checked(name, voi)
        for m in methods
        if r.get(m) is not None and not r[f"{m}_in_iqr"]
    ]

    def vs(method, ref):
        """Range of the differences (%) between two methods over the checked VOIs."""
        d = [
            100.0 * (r[method] / r[ref] - 1.0)
            for name, results in t5.items()
            for voi, r in results.items()
            if checked(name, voi) and r.get(method) is not None and r.get(ref) is not None
        ]
        return min(d), max(d)

    liver_a = t5["patient_4"]["Liver"]
    liver_a_nl = t5["patient_4"]["Liver (no lesions)"]
    c = liver_a["challenge"]

    lines = [
        "<!-- Generated by validation/snmmi/report.py. Do not edit by hand. -->",
        "",
        "# Validation against the SNMMI Dosimetry Challenge",
        "",
        "radot doses are computed on the clinical images of the "
        "[SNMMI 177Lu Dosimetry Challenge](https://doi.org/10.2967/jnumed.122.265094) "
        "and compared with the results of the challenge participants.",
        "",
        "**There is no ground truth**: the reference is the distribution of the participants' results "
        "(median, interquartile range IQR, quartile coefficient of dispersion QCD = (Q3 − Q1) / (Q3 + Q1)), "
        "from Brosch-Lenz et al., J Nucl Med 2023;64:1109, Supplemental Table 1. "
        "This is complementary to the [validation against Monte Carlo](7_validation_opengate.md): "
        "OpenGATE checks the physics with a known activity, the challenge checks the workflow on real "
        "images against the community.",
        "",
        "## The challenge",
        "",
        "Two patients treated with ¹⁷⁷Lu-DOTATATE, imaged with SPECT/CT at 4 time points "
        "(Siemens Intevo, xSPECT Quant, images in Bq/mL).",
        "",
        "| Patient | Sex | Injected (MBq) | Imaging days | Organs | Lesions |",
        "|---|---|--:|---|---|--:|",
    ]
    for p in PATIENTS.values():
        lines.append(
            f"| {p.challenge_id} | {p.sex} | {p.injected_mbq:.0f} | {', '.join(str(d) for d in p.days)} "
            f"| {', '.join(p.organs)} | {len(p.lesions)} |"
        )
    lines += [
        "",
        "Tasks:",
        "",
        "| Task | Input | Left to the participant | radot |",
        "|---|---|---|---|",
        "| T1 | serial SPECT/CT | segmentation, registration, fit, dose | not yet |",
        "| T2, T3 | planar (+ 1 SPECT/CT) | everything | out of scope |",
        "| T4 | serial SPECT/CT + VOIs at each time point | time integration, dose | not yet |",
        "| T5 | time-integrated activity (TIA) map + VOIs | dose only | below |",
        "",
        "## T5: dose from the TIA map",
        "",
        "- **Input**: the TIA map provided by the challenge (Bq·s/mL, registered to the first SPECT, "
        "unit `bqsml` in the sidecar) and the organ and lesion VOIs of the first time point.",
        "- **Mass**: from the CT of the first time point resampled to the TIA grid "
        "(density from the Hounsfield units).",
        "- **`svalue`**: TIA in the VOI × OpenDose self S-value, scaled to the VOI mass "
        "(ICRP-110 male phantom for A, female for B). Organs only: there are no S-values for lesions.",
        "- **`ldm`**: local energy deposition of the voxel TIA, mean in the VOI.",
    ]
    if has_gate:
        lines += [
            "- **`gate`**: [OpenGATE](https://github.com/OpenGATE/opengate) Monte Carlo (`opengate_dose_rate`, "
            "Lu-177, electrons and photons) with the TIA map as source and the CT on the TIA grid; the "
            "dose-rate map is used with `radot dose --dose-rate-map` (dose = map × TIA / simulated activity). "
            "Includes the photon self- and cross-doses. GATE unc.: statistical uncertainty (one standard "
            "deviation) of the VOI mean.",
        ]
    lines += [
        "- **Total kidney**: mass-weighted mean of the left and right kidney doses.",
        "- **Difference**: (radot − participants' median) / median; ✗ when radot is outside the "
        "participants' interquartile range [Q1–Q3].",
        "",
        "Absorbed doses in Gy for the whole treatment cycle (the TIA map is integrated over time).",
        "",
        "| Patient | VOI | Mass (g) | TIA (MBq·h) | "
        + " | ".join(f"`{m}` (Gy) | Diff." for m in methods)
        + (" | GATE unc." if has_gate else "")
        + " | Participants: median [Q1–Q3] (Gy) | QCD |",
        "|---|---|--:|--:|" + "--:|--:|" * len(methods) + ("--:|" if has_gate else "") + "--:|--:|",
    ]
    order = ["Liver", "Liver (no lesions)", "Spleen", "Kidney_L", "Kidney_R", "Kidneys (total)"]
    for name, results in t5.items():
        pid = PATIENTS[name].challenge_id
        for voi in sorted(results, key=lambda v: (order.index(v) if v in order else len(order), v)):
            r = results[voi]
            ch = r["challenge"]
            lines.append(
                f"| {pid} | {voi} | {r['mass_g']:.0f} | {r['tia_mbq_h']:.0f} | "
                + " | ".join(f"{fmt_dose(r, m)} | {fmt_diff(r, m)}" for m in methods)
                + (f" | {fmt_unc(r)}" if has_gate else "")
                + f" | {ch['median']:.2f} [{ch['q1']:.2f}–{ch['q3']:.2f}] | {ch['qcd_pct']:.0f} % |"
            )
    lines += [
        "",
        "Findings:",
        "",
        f"- All organs and lesions (except the liver of A, see below) are within "
        f"{max(abs(d) for d in diffs):.1f} % of the participants' median"
        + (f"; outside the IQR: {', '.join(outside)}." if outside else ", inside the IQR."),
        f"- **Liver of A**: the liver VOI contains lesion 1. The participants' spread (QCD {c['qcd_pct']:.0f} %, "
        f"the largest of T5) matches two liver definitions: the VOI as given gives {liver_a['svalue']:.2f} Gy "
        f"(median {c['median']:.2f} Gy), the liver without the lesions {liver_a_nl['svalue']:.2f} Gy "
        f"(Q1 {c['q1']:.2f} Gy). For B the lesions inside the liver are small.",
        "- The doses of A are a few percent below the median, those of B close to it. "
        + (
            "The Monte Carlo (`gate`, with photons and cross-doses) is as low as `svalue`: the cause is "
            "not the dose model but most likely the mass convention (see [Effect of the mass](#effect-of-the-mass))."
            if has_gate
            else "Possible causes: the mass convention (see [Effect of the mass](#effect-of-the-mass)), the "
            "cross-doses from neighbouring organs (not included in `svalue`), the phantom."
        ),
        "- `ldm` is 1–2 % below `svalue`: no photon self-dose.",
    ]
    if has_gate:
        lo_l, hi_l = vs("gate", "ldm")
        lo_s, hi_s = vs("gate", "svalue")
        lines += [
            f"- `gate` vs `ldm`: {lo_l:+.1f} to {hi_l:+.1f} %; `gate` vs `svalue` (organs): "
            f"{lo_s:+.1f} to {hi_s:+.1f} %. The Monte Carlo adds the photons and the cross-doses, "
            "and transports the electrons beyond the voxel and the VOI boundary.",
        ]
    lines += mass_section(t5, order)
    lines += [
        "",
        "## What this validates",
        "",
        "- The **dose step** on real images: TIA image as input, VOI masses from the CT, "
        "OpenDose S-values with mass scaling, local deposition"
        + (", and an external Monte Carlo dose-rate map (`--dose-rate-map`)." if has_gate else "."),
        "- Agreement with the community: radot falls in the narrow spread of the T5 participants "
        "(QCD 3–7 % for these targets).",
        "",
        "## What this does not validate",
        "",
        "- **Time integration** (T4) and **segmentation** (T1): not yet done.",
        "- **Accuracy**: the participants' median is not a ground truth; a common bias of the "
        "participants and radot would not be seen.",
        "- **Lesion doses with S-value methods**: no S-values for lesions.",
        "",
        "## Data and how to run",
        "",
        f"The images and VOIs are from the University of Michigan [Deep Blue collection]({DATA_URL}) "
        "(CC BY-NC 4.0, research only). They are not redistributed: download the zips by hand, see "
        f"[`validation/snmmi/`]({REPO_URL}).",
        "",
        "> We acknowledge Yuni Dewaraja at the University of Michigan for providing access to Lu-177 "
        "patient imaging data and segmentations obtained with support from R01EB022075 awarded by NIBIB "
        "and R01CA240706 awarded by NCI, NIH. Data were shared via the University of Michigan Deep Blue "
        "Data sharing repository.",
        "",
        "```bash",
        "python validation/snmmi/prepare_data.py    # extract the zips and convert to NIfTI",
        "python validation/snmmi/run_gate.py        # optional: OpenGATE dose-rate maps (needs opengate)",
        "python validation/snmmi/t5.py              # T5 doses, writes results/t5_radot.json",
        "pytest validation/snmmi                    # tests (T5 skipped without the data)",
        "python validation/snmmi/report.py          # regenerate this page",
        "```",
        "",
    ]
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument(
        "--check",
        action="store_true",
        help="Do not write; exit with status 1 if the page is out of date.",
    )
    args = parser.parse_args()

    page = render(json.loads(T5_RESULTS.read_text(encoding="utf-8")))

    if args.check:
        current = OUTPUT.read_text(encoding="utf-8") if OUTPUT.exists() else ""
        if current != page:
            print(f"{OUTPUT.relative_to(ROOT)} is out of date: run python validation/snmmi/report.py")
            return 1
        print(f"{OUTPUT.relative_to(ROOT)} is up to date")
        return 0

    OUTPUT.write_text(page, encoding="utf-8")
    print(f"Wrote {OUTPUT.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
