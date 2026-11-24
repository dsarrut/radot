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
T4_RESULTS = RESULTS_DIR / "t4_radot.json"
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


def t4_section(t4, order):
    """T4: radot integrate (time integration per VOI) then the dose, against the participants' T4."""
    from t4 import METHODS

    lines = [
        "",
        "## T4: time integration and dose from the serial SPECT/CT",
        "",
        "- **Input**: the 4 SPECT (Bq/mL) and the organ and lesion VOIs of each time point (no registration "
        "needed: each VOI is drawn on its own scan). Patient A's images are decay-corrected to the injection "
        "(`decay_correction` `admin`): radot brings them back to the activity at each acquisition time.",
        "- **Time integration**: `radot integrate`, one TIA per VOI, with each method below (4 time points: "
        "with the noise level estimated, `f2` has a standard error and `f4a` none; no model selection, the "
        "AICc needs more points).",
        "- **Dose**: organs, `radot dose -m svalue` on the TIA table (mass from the CT of time point 1, "
        "as T5); lesions (no S-value), the radot T5 `ldm` dose of the VOI scaled by the TIA ratio T4 / T5 "
        "map (marked *).",
        "- **TIA vs map**: radot TIA vs the TIA of the challenge map in the same VOI (the map is the "
        "organisers' time integration, the T5 input).",
        "- **Fit after**: `radot dose -m ldm --dose-rate` at each time point with the CT of that time point, "
        "then `radot integrate` on the 4 dose-rate images (VOI mean dose rate weighted by the masses of that "
        "time point, `--ct`), with the same method. Every VOI, organs and lesions (LDM, no S-value).",
        "",
        "| Method | `radot integrate` |",
        "|---|---|",
    ]
    for label, opts in METHODS.items():
        lines.append(f"| {label} | " + " ".join(f"`--{k.replace('_', '-')} {v}`" for k, v in opts.items()) + " |")
    labels = list(METHODS)

    def dose_cell(r, is_organ):
        if "failed" in r:
            return "failed"
        k = "svalue" if is_organ else "ldm_scaled"
        flag = "" if r[f"{k}_in_iqr"] else " ✗"
        return f"{r[k]:.2f}{'' if is_organ else '*'} ({r[f'{k}_diff_pct']:+.0f} %{flag})"

    lines += [
        "",
        "Absorbed doses (Gy) and difference with the participants' T4 median (✗: outside [Q1–Q3]):",
        "",
        "| Patient | VOI | " + " | ".join(labels) + " | Participants: median [Q1–Q3] (Gy) | QCD |",
        "|---|---|" + "--:|" * len(labels) + "--:|--:|",
    ]
    inside = {label: [0, 0] for label in labels}
    for name, results in t4.items():
        pid = PATIENTS[name].challenge_id
        for voi in sorted(results, key=lambda v: (order.index(v) if v in order else len(order), v)):
            per = results[voi]
            ch = per["challenge"]
            is_organ = not voi.startswith("Lesion")
            for label in labels:
                r = per[label]
                if "failed" not in r and checked(name, voi):
                    k = "svalue" if is_organ else "ldm_scaled"
                    inside[label][0] += r[f"{k}_in_iqr"]
                    inside[label][1] += 1
            lines.append(
                f"| {pid} | {voi} | " + " | ".join(dose_cell(per[label], is_organ) for label in labels)
                + f" | {ch['median']:.2f} [{ch['q1']:.2f}–{ch['q3']:.2f}] | {ch['qcd_pct']:.0f} % |"
            )
    lines += [
        "",
        "TIA vs the challenge map (%):",
        "",
        "| Patient | VOI | " + " | ".join(labels) + " |",
        "|---|---|" + "--:|" * len(labels),
    ]
    for name, results in t4.items():
        pid = PATIENTS[name].challenge_id
        for voi in sorted(results, key=lambda v: (order.index(v) if v in order else len(order), v)):
            per = results[voi]
            lines.append(f"| {pid} | {voi} | " + " | ".join(
                "failed" if "failed" in per[label] else f"{per[label]['tia_vs_t5_map_pct']:+.0f}" for label in labels
            ) + " |")
    def after_cell(r):
        if "failed" in r or r.get("ldm_fit_after") is None:
            return "failed"
        flag = "" if r["ldm_fit_after_in_iqr"] else " ✗"
        return f"{r['ldm_fit_after']:.2f} ({r['ldm_fit_after_diff_pct']:+.0f} %{flag})"

    def mass_spread(per):
        m = next((per[label].get("masses_g") for label in labels if per[label].get("masses_g")), None)
        return "—" if not m else f"{100 * (max(m) - min(m)) / (sum(m) / len(m)):.0f} %"

    lines += [
        "",
        "Fit after (LDM dose rate at each time point, then integrated), absorbed doses (Gy) and difference "
        "with the participants' T4 median; mass spread: (max − min) / mean of the VOI mass over the 4 time "
        "points (masks and CT of each scan):",
        "",
        "| Patient | VOI | " + " | ".join(labels) + " | Participants: median [Q1–Q3] (Gy) | Mass spread |",
        "|---|---|" + "--:|" * len(labels) + "--:|--:|",
    ]
    inside_after = {label: [0, 0] for label in labels}
    for name, results in t4.items():
        pid = PATIENTS[name].challenge_id
        for voi in sorted(results, key=lambda v: (order.index(v) if v in order else len(order), v)):
            per = results[voi]
            ch = per["challenge"]
            for label in labels:
                r = per[label]
                if r.get("ldm_fit_after") is not None and checked(name, voi):
                    inside_after[label][0] += r["ldm_fit_after_in_iqr"]
                    inside_after[label][1] += 1
            lines.append(
                f"| {pid} | {voi} | " + " | ".join(after_cell(per[label]) for label in labels)
                + f" | {ch['median']:.2f} [{ch['q1']:.2f}–{ch['q3']:.2f}] | {mass_spread(per)} |"
            )

    def after_vs_first(label):
        """Range of fit after / fit first (LDM, ldm_scaled) - 1 (%) over the checked VOIs."""
        d = [100.0 * (per[label]["ldm_fit_after"] / per[label]["ldm_scaled"] - 1.0)
             for n, res in t4.items() for v, per in res.items()
             if checked(n, v) and per[label].get("ldm_fit_after") is not None and "failed" not in per[label]]
        return f"{min(d):+.0f} to {max(d):+.0f} %"

    failed = [f"{PATIENTS[n].challenge_id} {v} ({label}): {per[label]['failed'].rstrip('.')}"
              for n, res in t4.items() for v, per in res.items() for label in labels if "failed" in per[label]]
    spreads = [100 * (max(m) - min(m)) / (sum(m) / len(m)) for res in t4.values() for v, per in res.items()
               if not v.startswith("Lesion") and (m := per["f4a"].get("masses_g"))]
    spread_lo, spread_hi = min(spreads), max(spreads)
    la = t4["patient_4"]["Liver"]["trapezoid extrap/fit2"]
    dt_a = la["times_h"][-1] - la["times_h"][-2]
    liver_dm = 100 * (la["masses_g"][-1] / la["masses_g"][-2] - 1)
    liver_diff = 100 * (la["ldm_fit_after"] / la["ldm_scaled"] - 1)
    lines += [
        "",
        "Findings:",
        "",
        "- Inside the participants' interquartile range (VOIs without the liver of A): "
        + ", ".join(f"`{label}` {a} / {b}" for label, (a, b) in inside.items()) + "; fit after: "
        + ", ".join(f"`{label}` {a} / {b}" for label, (a, b) in inside_after.items()) + ".",
        "- Fit after vs fit first (both LDM; fit first: `ldm*`), over the same VOIs: "
        + "; ".join(f"`{label}` {after_vs_first(label)}" for label in labels) + ". The engine is linear in "
        "activity: the difference comes from the geometry of each time point. The organ masses change by "
        f"{spread_lo:.0f} to {spread_hi:.0f} % between scans (each VOI is drawn again on each scan), and the "
        "fit after carries this change into the dose-rate curve. A tail fitted on two close points "
        f"(`extrap/fit2`, patient A: {dt_a:.0f} h between the last two scans) amplifies it: for the liver of A, "
        f"{liver_dm:+.0f} % of mass between scans 3 and 4 gives {liver_diff:+.0f} % of dose (fit after vs fit "
        "first). The geometry of each time point is more exact only when the VOIs are consistent between "
        "scans.",
        "- `f4a` (uptake then washout) is the closest to the participants for the lesions and the kidneys of "
        "B; `f2` (no uptake) underestimates the VOIs whose uptake is still visible at the first scan.",
        "- `trapezoid ramp/phys` overestimates every VOI: the physical-decay tail ignores the biological "
        "clearance after the last scan (as on the [synthetic curves](9_validation_time_integration.md)); "
        "a fitted tail (`extrap/fit2`) removes most of it.",
        "- The challenge TIA map and the participants differ: for some lesions the map TIA is far from radot "
        "while the radot dose is at the participants' median. The map is one time integration among others, "
        "not a reference.",
    ]
    if failed:
        lines.append("- Refused (radot does not extrapolate silently): " + "; ".join(failed) + ".")
    return lines


def render(t5, t4=None):
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
        "| T4 | serial SPECT/CT + VOIs at each time point | time integration, dose | "
        + ("below |" if t4 else "not yet |"),
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
    if t4:
        lines += t4_section(t4, order)
    lines += [
        "",
        "## What this validates",
        "",
        "- The **dose step** on real images: TIA image as input, VOI masses from the CT, "
        "OpenDose S-values with mass scaling, local deposition"
        + (", and an external Monte Carlo dose-rate map (`--dose-rate-map`)." if has_gate else "."),
        "- Agreement with the community: radot falls in the narrow spread of the T5 participants "
        "(QCD 3–7 % for these targets).",
    ]
    if t4:
        lines.append(
            "- The **time integration** on real serial SPECT/CT (T4): `radot integrate` per VOI, fit first "
            "(TIA, then `radot dose` on the table) and fit after (LDM dose rate at each time point with its "
            "CT, then integrated)."
        )
    lines += [
        "",
        "## What this does not validate",
        "",
        ("- **Segmentation** (T1): not yet done." if t4 else
         "- **Time integration** (T4) and **segmentation** (T1): not yet done."),
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
        "python validation/snmmi/t4.py              # T4 integrate + doses, writes results/t4_radot.json",
        "pytest validation/snmmi                    # tests (T4, T5 skipped without the data)",
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

    t4 = json.loads(T4_RESULTS.read_text(encoding="utf-8")) if T4_RESULTS.exists() else None
    page = render(json.loads(T5_RESULTS.read_text(encoding="utf-8")), t4)

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
