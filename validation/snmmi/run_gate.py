#!/usr/bin/env python3
"""
OpenGATE Monte Carlo dose rate for the SNMMI challenge, with the TIA map as source.

For each patient: the CT of time point 1 resampled to the TIA grid (the dose grid of
opengate_dose_rate is the source grid, and the CT is aligned by image centre: identical grids
avoid any shift), the TIA map as relative source, and an opengate_dose_rate configuration
(Lu177, VRT mode). The output dose-rate map (Gy/s for activity_bq distributed like the TIA)
is then used by t5.py with radot --dose-rate-map: D = map x TIA / activity.

Outputs in DATA_DIR/gate/<patient>/ (derived from the challenge data: never committed).
Needs opengate (not a radot dependency). Run prepare_data.py first.

Usage:
    python validation/snmmi/run_gate.py                    # both patients
    python validation/snmmi/run_gate.py patient_6 --activity 1e7 --threads 16
"""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import SimpleITK as sitk

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cases import DATA_DIR, PATIENTS, has_nifti, nifti_dir  # noqa: E402

from radot.resampling import resample_image  # noqa: E402

GATE_DIR = DATA_DIR / "gate"
OUTPUT_FOLDER = "output_tia_vrt"
DEFAULT_ACTIVITY_BQ = 1e7


def gate_dir(patient: str) -> Path:
    return GATE_DIR / patient


def gate_dose_rate_map(patient: str) -> Path:
    return gate_dir(patient) / OUTPUT_FOLDER / "output_dose.mhd"


def gate_config(patient: str) -> dict:
    return json.loads((gate_dir(patient) / "dose_rate_tia.json").read_text())


def has_gate(patient: str) -> bool:
    return gate_dose_rate_map(patient).exists()


def prepare_inputs(patient: str, activity_bq: float, threads: int) -> Path:
    """CT on the TIA grid, TIA source, HU tables and the opengate_dose_rate configuration."""
    import opengate

    d = gate_dir(patient)
    d.mkdir(parents=True, exist_ok=True)
    tia = sitk.ReadImage(str(nifti_dir(patient) / "tia.nii.gz"))
    ct = sitk.ReadImage(str(nifti_dir(patient) / "ct_scan1.nii.gz"))
    ct = resample_image(ct, image_type="ct", reference_image=tia)
    sitk.WriteImage(sitk.Cast(ct, sitk.sitkInt16), str(d / "ct_tia_grid.nii.gz"))
    sitk.WriteImage(sitk.Cast(tia, sitk.sitkFloat32), str(d / "tia_source.nii.gz"))
    tables = Path(opengate.__file__).parent / "data"
    for name in ("Schneider2000MaterialsTable.txt", "Schneider2000DensitiesTable.txt"):
        shutil.copy(tables / name, d / name)
    config = {
        "ct_image": "ct_tia_grid.nii.gz",
        "table_mat": "Schneider2000MaterialsTable.txt",
        "table_density": "Schneider2000DensitiesTable.txt",
        "density_tolerance_gcm3": 0.2,
        "activity_image": "tia_source.nii.gz",
        "radionuclide": "Lu177",
        "activity_bq": activity_bq,
        "number_of_threads": threads,
        "visu": False,
        "verbose": False,
        "mode": "vrt",
        "e_factor": 1,
        "output_folder": OUTPUT_FOLDER,
    }
    path = d / "dose_rate_tia.json"
    path.write_text(json.dumps(config, indent=2))
    return path


def run(patient: str, activity_bq: float, threads: int) -> Path:
    config = prepare_inputs(patient, activity_bq, threads)
    print(f"{patient}: opengate_dose_rate {config} ({activity_bq:.1e} Bq, {threads} threads)")
    subprocess.run(["opengate_dose_rate", config.name], cwd=config.parent, check=True)
    # The output grid is the source grid; keep the TIA geometry exactly (origin, direction)
    out = gate_dose_rate_map(patient)
    rate = sitk.ReadImage(str(out))
    tia = sitk.ReadImage(str(gate_dir(patient) / "tia_source.nii.gz"))
    if rate.GetSize() != tia.GetSize():
        raise RuntimeError(f"{out}: size {rate.GetSize()} differs from the TIA {tia.GetSize()}")
    # Sidecar read by radot --dose-rate-map: unit and simulated activity
    out.with_suffix(".json").write_text(
        json.dumps({"modality": "dose", "unit": "gy/s", "source_activity_mbq": activity_bq / 1e6}, indent=2)
    )
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.strip().splitlines()[0])
    parser.add_argument("patients", nargs="*", default=list(PATIENTS))
    parser.add_argument("--activity", type=float, default=DEFAULT_ACTIVITY_BQ, help="Simulated activity (Bq)")
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()
    for name in args.patients:
        if not has_nifti(name):
            sys.exit(f"No converted data for {name}: run prepare_data.py first.")
        print(f"Wrote {run(name, args.activity, args.threads)}")


if __name__ == "__main__":
    main()
