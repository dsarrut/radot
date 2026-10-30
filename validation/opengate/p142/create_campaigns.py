#!/usr/bin/env python3
"""
Build, and optionally run, the opengate_jobs campaigns of a p142 dose-rate simulation (VRT: electron
and photon runs). Run in p142/data/ (inputs downloaded with download_data.py).

Build only, e.g. to run the jobs on a cluster; writes campaigns/<radionuclide>/<part>/simulation.json,
with the input paths relative to it:

    python ../create_campaigns.py --radionuclide Y90 --activity 30e6
    opengate_jobs split campaigns/y90/e -j 100
    opengate_jobs run campaigns/y90/e --backend slurm --backend-options options.json
    opengate_jobs merge campaigns/y90/e
    (same for campaigns/y90/gamma_tle)

Build and run locally in parallel jobs, then assemble the output folder (merged images, e/ and
gamma_tle/ parts, opengate_version.json, reference.json):

    python ../create_campaigns.py --radionuclide Y90 --activity 30e6 --run -j 12 -o new/output_y90_vrt

run_references.zsh does it for all the VRT references. The configuration is
p142_dose_rate_<radionuclide>.json, except the activity (--activity, in Bq: more activity = more
statistics). reference.json records the simulated activity: create_benchmark_json.py gives every
result per MBq of it.

Electron run: bremsstrahlung included. opengate_dose_rate sets a 1 m production cut for all particles in
the electron run, so no bremsstrahlung photon is created (its energy is deposited along the electron
track). Here the photon cut is 2 mm (electrons and positrons still 1 m) and the dose actor is a
TLEDoseActor: electrons are scored as with the DoseActor, the bremsstrahlung photons with the
track-length estimator. To remove once opengate_dose_rate does it.
"""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

from box import Box
from opengate.contrib.dose.doserate import create_simulation, merge_vrt_dose_rate
from opengate.utility import g4_units

P142_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(P142_DIR.parent))
from benchmarks import opengate_version
PARTS = {"e": "e-", "gamma_tle": "gamma_tle"}


def include_bremsstrahlung(sim):
    """Electron run: create the bremsstrahlung photons and score them with TLE."""
    mm, m = g4_units.mm, g4_units.m
    pm = sim.physics_manager
    # in a region, an "all" cut overrides the per-particle cuts: unset it
    pm.set_production_cut("ct", "all", None)
    pm.set_production_cut("ct", "gamma", 2 * mm)
    pm.set_production_cut("ct", "electron", 1 * m)
    pm.set_production_cut("ct", "positron", 1 * m)

    old = sim.actor_manager.get_actor("dose")
    sim.actor_manager.remove_actor("dose")
    dose = sim.add_actor("TLEDoseActor", "dose")
    dose.output_filename = "output.mhd"
    for k in ("attached_to", "size", "spacing", "translation", "output_coordinate_system", "hit_type"):
        setattr(dose, k, getattr(old, k))
    dose.dose.active = True
    dose.dose_squared.active = True
    dose.dose_uncertainty.active = True


def create_campaigns(radionuclide: str, activity_bq: float, campaigns_dir: Path, threads: int = 1) -> dict:
    """Write the electron and photon campaigns; returns {part: campaign folder} and the e-factor."""
    with open(P142_DIR / f"p142_dose_rate_{radionuclide.lower()}.json") as f:
        config = json.load(f)
    if config.get("mode", "vrt") != "vrt":
        raise SystemExit(f"{radionuclide}: mode {config['mode']}, only VRT configurations are split here")
    e_factor = float(config.get("e_factor", 10.0))

    campaigns = {}
    for part, mode in PARTS.items():
        campaign = Path(campaigns_dir) / radionuclide.lower() / part
        if campaign.exists():
            shutil.rmtree(campaign)
        campaign.mkdir(parents=True)
        param = Box(config)
        param.mode = mode
        param.number_of_threads = threads
        param.visu = False
        param.activity_bq = int(activity_bq / e_factor if part == "e" else activity_bq)
        param.output_folder = str(campaign / "output")
        sim = create_simulation(param)
        if part == "e":
            include_bremsstrahlung(sim)
        sim.to_json_file(campaign)
        print(f"{campaign}/simulation.json  ({mode}, {param.activity_bq:.3g} Bq)")
        campaigns[part] = campaign

    reference = {
        "radionuclide": radionuclide,
        "activity_bq": float(activity_bq),
        "e_factor": e_factor,
        "bremsstrahlung": True,
        "opengate": opengate_version(),
    }
    folder = Path(campaigns_dir) / radionuclide.lower()
    (folder / "reference.json").write_text(json.dumps(reference, indent=2))
    (folder / "opengate_version.json").write_text(json.dumps(reference["opengate"], indent=2))
    print(f"{folder}/reference.json: {reference}")
    return campaigns, e_factor


def run_campaigns(campaigns: dict, e_factor: float, jobs: int, output: Path) -> Path:
    """Run each campaign in parallel jobs, merge both parts into output (as opengate_dose_rate does)."""
    for campaign in campaigns.values():
        # the opengate_jobs command, without relying on the installed entry points
        subprocess.run(
            [sys.executable, "-c", "from opengate.bin.opengate_jobs import go; go()",
             "run", str(campaign), "-j", str(jobs)],
            check=True,
        )
    output = Path(output)
    merge_vrt_dose_rate(campaigns["e"] / "output", campaigns["gamma_tle"] / "output", output, e_factor=e_factor)
    for part, campaign in campaigns.items():
        (output / part).mkdir(parents=True, exist_ok=True)
        for name in ("output_dose.mhd", "output_dose.raw", "output_dose_uncertainty.mhd",
                     "output_dose_uncertainty.raw", "stats.txt"):
            shutil.copy2(campaign / "output" / name, output / part / name)
    folder = campaigns["e"].parent
    for name in ("reference.json", "opengate_version.json"):
        shutil.copy2(folder / name, output / name)
    print(f"Output: {output}")
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--radionuclide", required=True, help="e.g. Y90, Lu177, I131")
    parser.add_argument("--activity", type=float, required=True, help="total simulated activity (Bq)")
    parser.add_argument("--threads", type=int, default=1, help="threads per job (default 1)")
    parser.add_argument("--campaigns", default="campaigns", help="campaigns folder (default campaigns/)")
    parser.add_argument("--run", action="store_true", help="run the campaigns locally, then merge")
    parser.add_argument("-j", "--jobs", type=int, default=12, help="parallel jobs with --run (default 12)")
    parser.add_argument("-o", "--output", help="output folder with --run (default new/output_<radionuclide>_vrt)")
    args = parser.parse_args()

    campaigns, e_factor = create_campaigns(args.radionuclide, args.activity, Path(args.campaigns), args.threads)
    if args.run:
        output = args.output or f"new/output_{args.radionuclide.lower()}_vrt"
        run_campaigns(campaigns, e_factor, args.jobs, Path(output))


if __name__ == "__main__":
    main()
