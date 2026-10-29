#!/usr/bin/env python3
"""
SNMMI challenge task 5: absorbed doses from the provided TIA map and VOIs (dose step only).

For each patient and VOI of time point 1:
- svalue: TIA in the VOI x OpenDose self S-value, mass-scaled (organs only, no S-value for lesions)
- ldm: local deposition of the voxel TIA, mass-weighted mean in the VOI
- gate: OpenGATE Monte Carlo dose-rate map with the TIA map as source (run_gate.py), used with
  radot --dose-rate-map: D = map x TIA / simulated activity (only when the map exists)
Masses from the CT of time point 1 resampled to the TIA grid (density = 1 + HU/1000).

Two liver variants: the liver VOI as given (it contains the liver lesions) and the liver
without the lesion VOIs. Total kidney: mass-weighted mean of the left and right kidney doses.

Compared with the participants' T5 doses (challenge_stats.json). Writes results/t5_radot.json
and prints a table. Run prepare_data.py first.
"""

import json
import sys
from pathlib import Path

import numpy as np
import SimpleITK as sitk

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(1, str(Path(__file__).resolve().parents[1] / "opengate"))
from cases import (  # noqa: E402
    CHALLENGE_VOI,
    ORGAN_NAMES,
    PATIENTS,
    RESULTS_DIR,
    challenge_stats,
    has_nifti,
    nifti_dir,
)

from run_gate import gate_config, gate_dose_rate_map, has_gate  # noqa: E402

from radot.dose import compute_roi_dose  # noqa: E402
from radot.log import configure_logging  # noqa: E402
from radot.resampling import resample_image  # noqa: E402

ISOTOPE = "Lu177"
METHODS = ("svalue", "ldm", "gate")


def read_mask(path: Path) -> np.ndarray:
    return sitk.GetArrayFromImage(sitk.ReadImage(str(path))) > 0


def voi_dose(tia, tia_path, ct, mask, roi_name, phantom, methods, gate=None):
    """Doses of one VOI (binary mask) for each method; None when the method does not apply.
    gate: (dose-rate map, simulated activity MBq, output folder) for the 'gate' method."""
    labels = sitk.GetImageFromArray(mask.astype(np.uint8))
    labels.CopyInformation(tia)
    out = {}
    for method in methods:
        if method == "gate":
            if gate is None:
                out["gate"] = None
                continue
            from benchmarks import gate_roi_uncertainty

            rate_map, activity_mbq, gate_out = gate
            res = compute_roi_dose(
                image=tia,
                image_path=tia_path,
                labels=labels,
                labels_map={1: roi_name},
                ct_image=ct,
                method="svalue",
                isotope=ISOTOPE,
                dose_rate_map=rate_map,
                dose_rate_map_activity_mbq=activity_mbq,
            )
            out["gate"] = res[roi_name]["mean_dose_gy"]
            out["gate_unc_pct"] = 100.0 * gate_roi_uncertainty(gate_out, mask, ISOTOPE)
            continue
        res = compute_roi_dose(
            image=tia,
            image_path=tia_path,
            labels=labels,
            labels_map={1: roi_name},
            ct_image=ct,
            method=method,
            isotope=ISOTOPE,
            method_params={"phantom": phantom},
        )
        if method == "ldm":
            res = res[0]
        r = res.get(roi_name)
        out[method] = None if r is None else r["mean_dose_gy"]
        if r is not None:
            out["mass_g"] = r["roi_mass_g"]
            out["volume_ml"] = r["roi_volume_ml"]
            out["tia_mbq_h"] = r["tia_mbq_h"]
    return out


def run_patient(name: str) -> dict:
    patient = PATIENTS[name]
    d = nifti_dir(name)
    tia_path = d / "tia.nii.gz"
    tia = sitk.ReadImage(str(tia_path))
    ct = resample_image(sitk.ReadImage(str(d / "ct_scan1.nii.gz")), image_type="ct", reference_image=tia)

    gate = None
    if has_gate(name):
        # The GATE map is on the TIA grid: copy the TIA geometry (MHD header rounding)
        rate_map = sitk.ReadImage(str(gate_dose_rate_map(name)))
        rate_map.CopyInformation(tia)
        gate = (rate_map, gate_config(name)["activity_bq"] / 1e6, gate_dose_rate_map(name).parent)

    masks = {v: read_mask(d / "masks" / f"{v}.nii.gz") for v in patient.organs + patient.lesions}
    any_lesion = np.zeros_like(masks[patient.organs[0]])
    for v in patient.lesions:
        any_lesion |= masks[v]

    results = {}
    for v in patient.organs:
        results[ORGAN_NAMES[v]] = voi_dose(tia, tia_path, ct, masks[v], ORGAN_NAMES[v], patient.sex, METHODS, gate)
    results["Liver (no lesions)"] = voi_dose(
        tia, tia_path, ct, masks["liver"] & ~any_lesion, "Liver", patient.sex, METHODS, gate
    )
    for v in patient.lesions:
        # No S-value for a lesion: LDM and GATE only
        results[v.capitalize()] = voi_dose(tia, tia_path, ct, masks[v], v.capitalize(), patient.sex, ["ldm", "gate"], gate)
        results[v.capitalize()]["svalue"] = None

    kl, kr = results["Kidney_L"], results["Kidney_R"]
    m = kl["mass_g"] + kr["mass_g"]
    results["Kidneys (total)"] = {
        k: None if kl[k] is None else (kl[k] * kl["mass_g"] + kr[k] * kr["mass_g"]) / m for k in METHODS
    } | {
        "mass_g": m,
        "volume_ml": kl["volume_ml"] + kr["volume_ml"],
        "tia_mbq_h": kl["tia_mbq_h"] + kr["tia_mbq_h"],
    }
    return results


def compare(all_results: dict) -> dict:
    """Add the participants' T5 statistics and the radot difference with their median."""
    for name, results in all_results.items():
        stats = challenge_stats("T5", name)
        for voi, r in results.items():
            s = stats[CHALLENGE_VOI[voi]]
            r["challenge"] = {k: s[k] for k in ("q1", "median", "q3", "qcd_pct")}
            for method in METHODS:
                if r[method] is not None:
                    r[f"{method}_diff_pct"] = 100.0 * (r[method] / s["median"] - 1.0)
                    r[f"{method}_in_iqr"] = bool(s["q1"] <= r[method] <= s["q3"])
    return all_results


def print_table(all_results: dict) -> None:
    def cell(r, method):
        if r[method] is None:
            return "-"
        flag = "" if r[f"{method}_in_iqr"] else " *"
        return f"{r[method]:.2f} ({r[f'{method}_diff_pct']:+.0f}%){flag}"

    print(
        f"\n{'Pat.':4} {'VOI':19} {'Mass (g)':>8} {'TIA (MBq.h)':>11} "
        f"{'svalue (Gy)':>15} {'ldm (Gy)':>15} {'gate (Gy)':>15} {'T5 median [IQR] (Gy)':>24} {'QCD':>5}"
    )
    for name, results in all_results.items():
        pid = PATIENTS[name].challenge_id
        for voi, r in results.items():
            c = r["challenge"]
            iqr = f"{c['median']:.2f} [{c['q1']:.2f}-{c['q3']:.2f}]"
            print(
                f"{pid:4} {voi:19} {r['mass_g']:8.1f} {r['tia_mbq_h']:11.1f} "
                f"{cell(r, 'svalue'):>15} {cell(r, 'ldm'):>15} {cell(r, 'gate'):>15} {iqr:>24} {c['qcd_pct']:4.0f}%"
            )
    print("(difference with the participants' median; * outside the participants' interquartile range)")


def main() -> None:
    configure_logging(quiet=True)
    all_results = {}
    for name in PATIENTS:
        if not has_nifti(name):
            sys.exit(f"No converted data for {name}: run prepare_data.py first.")
        all_results[name] = run_patient(name)
    compare(all_results)
    print_table(all_results)
    RESULTS_DIR.mkdir(exist_ok=True)
    out = RESULTS_DIR / "t5_radot.json"
    out.write_text(json.dumps(all_results, indent=2))
    print(f"\nSaved {out}")


if __name__ == "__main__":
    main()
