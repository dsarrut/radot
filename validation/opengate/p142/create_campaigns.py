#!/usr/bin/env python3
"""
Build the opengate_jobs campaigns of a p142 dose-rate simulation (VRT: electron and photon runs),
to run them on a cluster instead of with opengate_dose_rate.

Run in p142/data/ (inputs downloaded with download_data.py). For each part, writes
campaigns/<radionuclide>/<part>/simulation.json, with the input paths relative to it:

    python ../create_campaigns.py --radionuclide Y90 --activity 30e6
    opengate_jobs split campaigns/y90/e -j 100
    opengate_jobs run campaigns/y90/e --backend slurm --backend-options options.json
    opengate_jobs merge campaigns/y90/e
    (same for campaigns/y90/gamma_tle)
    opengate_dose_rate --merge-only campaigns/y90/e/output campaigns/y90/gamma_tle/output \\
        --e-factor 1 -o output_y90_vrt
    cp campaigns/y90/opengate_version.json output_y90_vrt/

The opengate version (and git commit) is written to campaigns/<radionuclide>/opengate_version.json;
create_benchmark_json.py copies it into the benchmark JSON when it is in the output folder.

The configuration is p142_dose_rate_<radionuclide>.json, except the activity (--activity, in Bq:
more activity = more statistics; the benchmark JSON gives every result per MBq of the simulated activity).

Electron run: bremsstrahlung included. opengate_dose_rate sets a 1 m production cut for all particles in
the electron run, so no bremsstrahlung photon is created (its energy is deposited along the electron
track). Here the photon cut is 2 mm (electrons and positrons still 1 m) and the dose actor is a
TLEDoseActor: electrons are scored as with the DoseActor, the bremsstrahlung photons with the
track-length estimator. To remove once opengate_dose_rate does it.
"""

import argparse
import json
import sys
from pathlib import Path

from box import Box
from opengate.contrib.dose.doserate import create_simulation
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


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--radionuclide", required=True, help="e.g. Y90, Lu177, I131")
    parser.add_argument("--activity", type=float, required=True, help="total simulated activity (Bq)")
    parser.add_argument("--threads", type=int, default=1, help="threads per job (default 1)")
    parser.add_argument("--output", default="campaigns", help="campaigns folder (default campaigns/)")
    args = parser.parse_args()

    with open(P142_DIR / f"p142_dose_rate_{args.radionuclide.lower()}.json") as f:
        config = json.load(f)
    if config.get("mode", "vrt") != "vrt":
        raise SystemExit(f"{args.radionuclide}: mode {config['mode']}, only VRT configurations are split here")
    e_factor = float(config.get("e_factor", 10.0))

    for part, mode in PARTS.items():
        campaign = Path(args.output) / args.radionuclide.lower() / part
        campaign.mkdir(parents=True, exist_ok=True)
        param = Box(config)
        param.mode = mode
        param.number_of_threads = args.threads
        param.visu = False
        param.activity_bq = int(args.activity / e_factor if part == "e" else args.activity)
        param.output_folder = str(campaign / "output")
        sim = create_simulation(param)
        if part == "e":
            include_bremsstrahlung(sim)
        sim.to_json_file(campaign)
        print(f"{campaign}/simulation.json  ({mode}, {param.activity_bq:.3g} Bq)")
    version_file = Path(args.output) / args.radionuclide.lower() / "opengate_version.json"
    with open(version_file, "w") as f:
        json.dump(opengate_version(), f, indent=2)
    print(f"{version_file}: {opengate_version()}")


if __name__ == "__main__":
    main()
