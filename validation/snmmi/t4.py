#!/usr/bin/env python3
"""
SNMMI challenge task 4: absorbed doses from the 4 SPECT/CT and the VOIs of each time point (time
integration + dose).

For each patient, VOI and time-integration method (METHODS), the radot chain fit first at ROI level:
- `radot integrate`: the VOI activity at each time point (mask of that scan, activity at the
  acquisition time: patient A is decay-corrected to the injection), then the TIA;
- organs: `radot dose -m svalue` on the TIA table (OpenDose self S-value, mass from the CT of time
  point 1, as T5);
- every VOI (lesions have no S-value): `ldm_scaled` = the T5 radot LDM dose of the VOI × TIA(T4) /
  TIA(T5 map), i.e. the radot LDM dose per unit TIA of the same VOI and CT (results/t5_radot.json).

Compared with the participants' T4 doses (challenge_stats.json), and the TIA with the TIA of the
challenge map (T5, made by the organisers). Writes results/t4_radot.json and prints a table. Run
prepare_data.py and t5.py first.
"""

import json
import sys
from pathlib import Path

import numpy as np
import SimpleITK as sitk

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cases import (  # noqa: E402
    CHALLENGE_VOI,
    ORGAN_NAMES,
    PATIENTS,
    RESULTS_DIR,
    SCANS,
    challenge_stats,
    has_t4_nifti,
    mask_path,
    nifti_dir,
    spect_path,
)

from radot.dose import compute_roi_dose_tia_table  # noqa: E402
from radot.integrate import integrate_rois  # noqa: E402
from radot.log import configure_logging  # noqa: E402
from radot.resampling import resample_image  # noqa: E402

# label -> radot integrate options (4 time points: f2 with the noise estimated; f4a without SE)
METHODS = {
    "f2": {"model": "f2", "error_model": "proportional"},
    "f4a": {"model": "f4a", "error_model": "proportional"},
    "trapezoid ramp/phys": {"model": "trapezoid", "head": "ramp", "tail": "phys"},
    "trapezoid extrap/fit2": {"model": "trapezoid", "head": "extrap", "tail": "fit", "tail_points": 2},
}
T5_RESULTS = RESULTS_DIR / "t5_radot.json"


def read_mask(path: Path) -> np.ndarray:
    return sitk.GetArrayFromImage(sitk.ReadImage(str(path))) > 0


def label_image(mask: np.ndarray, like: sitk.Image) -> sitk.Image:
    img = sitk.GetImageFromArray(mask.astype(np.uint8))
    img.CopyInformation(like)
    return img


def voi_t4(spects, spect_paths, masks, ct, roi, is_organ, sex, t5_voi):
    """TIA and doses of one VOI for every method. masks: one boolean array per scan."""
    labels = [label_image(m, s) for m, s in zip(masks, spects)]
    out = {}
    for label, opts in METHODS.items():
        try:
            doc = integrate_rois(spects, labels, image_paths=spect_paths, labels_map={1: roi}, rois=[roi], **opts)
        except ValueError as e:
            out[label] = {"failed": str(e)}
            continue
        r = doc["rois"][roi]
        res = {
            "tia_mbq_h": r["tia_mbq_h"],
            "tia_se_mbq_h": r["tia_se_mbq_h"],
            "activities_mbq": r["activities_mbq"],
            "times_h": r["times_h"],
            "tia_vs_t5_map_pct": 100.0 * (r["tia_mbq_h"] / t5_voi["tia_mbq_h"] - 1.0),
            "ldm_scaled": t5_voi["ldm"] * r["tia_mbq_h"] / t5_voi["tia_mbq_h"],
            "svalue": None,
        }
        if is_organ:
            d = compute_roi_dose_tia_table(doc, labels[0], labels_map={1: roi}, ct_image=ct,
                                           method_params={"phantom": sex}, rois=[roi])[roi]
            res["svalue"] = d["mean_dose_gy"]
            res["svalue_fit_se"] = d["mean_dose_fit_se_gy"]
            res["mass_g"] = d["roi_mass_g"]
        out[label] = res
    return out


def run_patient(name: str, t5: dict) -> dict:
    p = PATIENTS[name]
    paths = [spect_path(name, k) for k in SCANS]
    spects = [sitk.ReadImage(str(x)) for x in paths]
    ct = resample_image(sitk.ReadImage(str(nifti_dir(name) / "ct_scan1.nii.gz")), image_type="ct",
                        reference_image=spects[0])
    masks = {v: [read_mask(mask_path(name, k, v)) for k in SCANS] for v in p.organs + p.lesions}
    results = {}
    for v in p.organs:
        roi = ORGAN_NAMES[v]
        results[roi] = voi_t4(spects, paths, masks[v], ct, roi, True, p.sex, t5[roi])
    no_lesion = []
    for i in range(len(SCANS)):
        m = masks["liver"][i].copy()
        for v in p.lesions:
            m &= ~masks[v][i]
        no_lesion.append(m)
    results["Liver (no lesions)"] = voi_t4(spects, paths, no_lesion, ct, "Liver", True, p.sex, t5["Liver (no lesions)"])
    for v in p.lesions:
        results[v.capitalize()] = voi_t4(spects, paths, masks[v], ct, v.capitalize(), False, p.sex, t5[v.capitalize()])
    # total kidney: mass-weighted mean of the two kidneys (masses of T5, CT of time point 1)
    kl, kr = results["Kidney_L"], results["Kidney_R"]
    ml, mr = t5["Kidney_L"]["mass_g"], t5["Kidney_R"]["mass_g"]
    total = {}
    for label in METHODS:
        a, b = kl[label], kr[label]
        if "failed" in a or "failed" in b:
            total[label] = {"failed": a.get("failed") or b.get("failed")}
            continue
        total[label] = {
            "tia_mbq_h": a["tia_mbq_h"] + b["tia_mbq_h"],
            "tia_vs_t5_map_pct": 100.0 * ((a["tia_mbq_h"] + b["tia_mbq_h"]) / t5["Kidneys (total)"]["tia_mbq_h"] - 1),
            "svalue": (a["svalue"] * ml + b["svalue"] * mr) / (ml + mr),
            "ldm_scaled": (a["ldm_scaled"] * ml + b["ldm_scaled"] * mr) / (ml + mr),
        }
    results["Kidneys (total)"] = total
    return results


def compare(all_results: dict) -> dict:
    """Add the participants' T4 statistics and the radot difference with their median."""
    for name, results in all_results.items():
        stats = challenge_stats("T4", name)
        for voi, per_method in results.items():
            s = stats[CHALLENGE_VOI[voi]]
            per_method["challenge"] = {k: s[k] for k in ("q1", "median", "q3", "qcd_pct")}
            for label in METHODS:
                r = per_method[label]
                for dose in ("svalue", "ldm_scaled"):
                    if r.get(dose) is not None:
                        r[f"{dose}_diff_pct"] = 100.0 * (r[dose] / s["median"] - 1.0)
                        r[f"{dose}_in_iqr"] = bool(s["q1"] <= r[dose] <= s["q3"])
    return all_results


def print_table(all_results: dict) -> None:
    print(f"\n{'Pat.':4} {'VOI':19} {'method':22} {'TIA vs map':>10} {'svalue (Gy)':>15} {'ldm* (Gy)':>15} "
          f"{'T4 median [IQR] (Gy)':>22}")
    for name, results in all_results.items():
        pid = PATIENTS[name].challenge_id
        for voi, per_method in results.items():
            c = per_method["challenge"]
            iqr = f"{c['median']:.2f} [{c['q1']:.2f}-{c['q3']:.2f}]"
            for label in METHODS:
                r = per_method[label]
                if "failed" in r:
                    print(f"{pid:4} {voi:19} {label:22} failed: {r['failed'][:70]}")
                    continue

                def cell(k):
                    if r.get(k) is None:
                        return "-"
                    return f"{r[k]:.2f} ({r[k + '_diff_pct']:+.0f}%){'' if r[k + '_in_iqr'] else ' *'}"

                print(f"{pid:4} {voi:19} {label:22} {r['tia_vs_t5_map_pct']:+9.1f}% {cell('svalue'):>15} "
                      f"{cell('ldm_scaled'):>15} {iqr:>22}")
    print("(TIA vs map: radot T4 TIA vs the TIA of the challenge map in the VOI; ldm*: T5 radot LDM dose x "
          "TIA ratio; * outside the participants' interquartile range)")


def main() -> None:
    configure_logging(quiet=True)
    if not T5_RESULTS.exists():
        sys.exit(f"No {T5_RESULTS}: run t5.py first.")
    t5 = json.loads(T5_RESULTS.read_text())
    all_results = {}
    for name in PATIENTS:
        if not has_t4_nifti(name):
            sys.exit(f"No converted T4 data for {name}: run prepare_data.py first.")
        all_results[name] = run_patient(name, t5[name])
    compare(all_results)
    print_table(all_results)
    RESULTS_DIR.mkdir(exist_ok=True)
    out = RESULTS_DIR / "t4_radot.json"
    out.write_text(json.dumps(all_results, indent=2))
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
