#!/usr/bin/env python3
"""
Measure the uncertainty of the organ mean OpenGATE dose rate with independent simulations (batch method).

OpenGATE gives a history-by-history uncertainty per voxel, not for an organ: the organ uncertainty
is estimated from the voxel maps (benchmarks.gate_roi_uncertainty). To check this estimate, the
reference simulation is repeated n times with different random seeds; the relative standard
deviation of the n organ means is the uncertainty of one simulation.

    python validation/opengate/batch_uncertainty.py --run      # run the missing batches, then analyse
    python validation/opengate/batch_uncertainty.py            # analyse the existing batches

Simulations go to p142/data/batches/run_NN/ (about 1 min each for Lu177 on 8 threads); the result
is written to p142/p142_<radionuclide>_batch_uncertainty.json (versioned, used by the tests).
"""

import argparse
import json
from pathlib import Path

import numpy as np
import SimpleITK as sitk

from benchmarks import (
    BENCH_DIR,
    CROSS_TARGETS,
    DATA_DIR,
    cross_dose,
    download_data,
    gate_output_dir,
    gate_roi_uncertainty,
    load_benchmark,
    organ_mask,
    run_opengate,
)

BATCH_DIR = DATA_DIR / "batches"


def batch_file(isotope: str) -> Path:
    return BENCH_DIR / f"p142_{isotope.lower()}_batch_uncertainty.json"


def run_batches(isotope: str, n: int) -> None:
    for i in range(1, n + 1):
        sim_dir = BATCH_DIR / f"run_{i:02d}"
        if (gate_output_dir(isotope, sim_dir) / "output_dose.mhd").exists():
            continue
        print(f"batch {i}/{n} ...")
        run_opengate(isotope, DATA_DIR, sim_dir)


def analyse(isotope: str) -> dict:
    bench = load_benchmark(isotope)
    outputs = [gate_output_dir(isotope, d) for d in sorted(BATCH_DIR.glob("run_*"))]
    outputs = [o for o in outputs if (o / "output_dose.mhd").exists()]
    n = len(outputs)
    if n < 2:
        raise RuntimeError(f"need at least 2 batches in {BATCH_DIR}, found {n}")

    result = {
        "description": f"Uncertainty of the organ mean OpenGATE dose rate, measured with {n} independent "
        f"simulations of p142_dose_rate_{isotope.lower()}.json (batch method, see batch_uncertainty.py)",
        "radionuclide": isotope,
        "n_batches": n,
        "organs": {},
    }
    # The source organ (self-dose), then the cross-dose targets
    organs = [(bench["roi_name"], "", bench["opengate"])]
    organs += [(t, t, cross_dose(bench, t)["opengate"]) for t in CROSS_TARGETS]
    for organ, target, ref in organs:
        mask = organ_mask(DATA_DIR, target)
        rates, estimates = [], []
        for out in outputs:
            dose = sitk.ReadImage(str(out / "output_dose.mhd"))
            rates.append(float(sitk.GetArrayViewFromImage(dose)[mask].mean()))
            estimates.append(gate_roi_uncertainty(out, mask, isotope))
        rates = np.array(rates)
        measured = float(rates.std(ddof=1) / rates.mean())
        mean = float(rates.mean())
        result["organs"][organ] = {
            "mean_dose_rate_gy_s": rates.tolist(),
            # uncertainty of one simulation, and the uncertainty of this std (n samples)
            "measured_uncertainty_fraction": measured,
            "measured_uncertainty_precision_fraction": measured / np.sqrt(2 * (n - 1)),
            # estimate from the voxel maps (gate_roi_uncertainty), mean over the batches
            "estimated_uncertainty_fraction": float(np.mean(estimates)),
            # mean of the batches vs the reference simulation (benchmark JSON)
            "batch_mean_vs_reference_fraction": (mean - ref["mean_dose_rate_gy_s"]) / ref["mean_dose_rate_gy_s"],
            "batch_mean_uncertainty_fraction": measured / np.sqrt(n),
        }
    versions = []
    for out in outputs:
        f = out / "opengate_version.json"
        versions.append(json.loads(f.read_text()) if f.exists() else None)
    result["opengate"] = versions
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--isotope", default="Lu177")
    parser.add_argument("-n", type=int, default=20, help="number of batches (default 20)")
    parser.add_argument("--run", action="store_true", help="run the missing batches first")
    args = parser.parse_args()

    download_data()
    if args.run:
        run_batches(args.isotope, args.n)
    result = analyse(args.isotope)
    with open(batch_file(args.isotope), "w") as f:
        json.dump(result, f, indent=2)
    print(f"{result['n_batches']} batches, {args.isotope}")
    for organ, r in result["organs"].items():
        print(
            f"  {organ:9s} measured {r['measured_uncertainty_fraction']:.3%} "
            f"(+/- {r['measured_uncertainty_precision_fraction']:.3%}), "
            f"estimated {r['estimated_uncertainty_fraction']:.3%}; "
            f"batch mean vs reference {r['batch_mean_vs_reference_fraction']:+.3%} "
            f"(+/- {r['batch_mean_uncertainty_fraction']:.3%})"
        )
    print(f"Saved to {batch_file(args.isotope)}")


if __name__ == "__main__":
    main()
