"""
Comparison of radot organ doses with OpenGATE Monte Carlo benchmarks.

Single source of truth for the validation: the cases and acceptance ranges below
are used by both test_opengate.py and report.py (docs/7_validation.md).

A benchmark JSON (see p142/p142_*_opengate_benchmark.json) gives the mean
OpenGATE dose rate at t=0 in a source organ (self-dose) and in one target organ
(cross-dose), for a uniform activity in the source organ. The OpenGATE dose rate
is integrated with the same mono-exponential clearance (Teff) as radot, so the
comparison tests the absorbed dose per decay, not the kinetics.

The organ is modelled as a single voxel with the benchmark volume, mass and
total activity; this gives the same organ-level results as the full images
(within ~0.5 %).
"""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import SimpleITK as sitk

from radot.db import get_reference_mass, get_s_value, resolve_organ_name
from radot.dose import compute_roi_dose, dose_hanscheid2018, dose_madsen2018

BENCH_DIR = Path(__file__).resolve().parent / "p142"

BENCHMARK_FILES = {
    "Lu177": "p142_lu177_opengate_benchmark.json",
    "Y90": "p142_y90_opengate_benchmark.json",
    "I131": "p142_i131_opengate_benchmark.json",
    "Ac225": "p142_ac225_opengate_benchmark.json",
}

METHOD_NAMES = {
    "madsen2018": "Madsen 2018",
    "hanscheid2017": "Hänscheid 2017",
    "hanscheid2018": "Hänscheid 2018",
    "ldm": "LDM",
}


@dataclass(frozen=True)
class Case:
    isotope: str
    method: str
    # Acceptance range for |radot - GATE| / GATE, in %; None = shown but not checked
    accept: Optional[Tuple[float, float]]
    note: str = ""
    cascade: bool = False
    cross: bool = False

    @property
    def label(self) -> str:
        label = METHOD_NAMES[self.method]
        return label + " (cascade)" if self.cascade else label

    def __str__(self) -> str:
        kind = "cross" if self.cross else "self"
        cascade = "-cascade" if self.cascade else ""
        return f"{self.isotope}-{self.method}{cascade}-{kind}"


H2018_NOTE = "t = T<sub>eff</sub> / 2: same as Madsen 2018 by construction"
AC225_PARENT_NOTE = "parent only; GATE tracks the full decay chain"

CASES: List[Case] = [
    # Self-dose (liver <- liver)
    Case("Lu177", "madsen2018", (0, 5)),
    Case("Lu177", "hanscheid2017", (0, 15), "assumes density 1 g/mL and a fixed 0.125 factor"),
    Case("Lu177", "hanscheid2018", (0, 5), H2018_NOTE),
    Case("Lu177", "ldm", (0, 1)),
    Case("Y90", "madsen2018", (0, 10)),
    Case("Y90", "ldm", (0, 5), "edge effect: β range ~11 mm"),
    Case("I131", "madsen2018", (0, 5)),
    Case("I131", "ldm", (15, 35), "LDM ignores the 364 keV photon self-dose (expected)"),
    Case("Ac225", "madsen2018", None, AC225_PARENT_NOTE),
    Case("Ac225", "ldm", None, AC225_PARENT_NOTE),
    Case("Ac225", "ldm", (0, 5), "with daughters in secular equilibrium", cascade=True),
    # Cross-dose (right kidney <- liver): photon dose, sensitive to phantom vs patient geometry
    Case("Lu177", "madsen2018", (0, 50), cross=True),
    Case("Lu177", "hanscheid2018", (0, 50), cross=True),
]


def load_benchmark(isotope: str, bench_dir: Path = BENCH_DIR) -> dict:
    with open(Path(bench_dir) / BENCHMARK_FILES[isotope]) as f:
        return json.load(f)


def gate_integrated_dose_gy(dose_rate_gy_s: float, teff_h: float) -> float:
    """Integrate a dose rate at t=0 with mono-exponential clearance of half-life teff_h."""
    k = np.log(2) / (teff_h * 3600.0)
    return dose_rate_gy_s / k


def organ_proxy(bench: dict, time_from_injection_h: float):
    """Single-voxel activity, label and CT images with the organ volume, mass and activity."""
    roi_volume_ml = float(bench["roi_volume_ml"])
    roi_mass_g = float(bench["roi_mass_g"])
    s = (roi_volume_ml * 1000.0) ** (1.0 / 3.0)

    activity = sitk.Image([1, 1, 1], sitk.sitkFloat32)
    activity.SetSpacing((s, s, s))
    activity.SetMetaData("unit", "Bq")
    activity.SetMetaData("time_from_injection", str(time_from_injection_h))
    activity[0, 0, 0] = float(bench["total_activity_bq"])

    labels = sitk.Image([1, 1, 1], sitk.sitkUInt16)
    labels.SetSpacing((s, s, s))
    labels[0, 0, 0] = 1

    ct = sitk.Image([1, 1, 1], sitk.sitkInt16)
    ct.SetSpacing((s, s, s))
    ct[0, 0, 0] = int(round((roi_mass_g / roi_volume_ml - 1.0) * 1000.0))
    return activity, labels, ct


def radot_self_dose_gy(bench: dict, isotope: str, method: str, cascade: bool = False) -> float:
    teff_h = float(bench["teff_h"])
    roi = bench["roi_name"]
    phantom = bench.get("phantom", "male")

    if method == "ldm":
        activity, labels, ct = organ_proxy(bench, 0.0)
        results, _ = compute_roi_dose(
            image=activity,
            labels=labels,
            method="ldm",
            isotope=isotope,
            ct_image=ct,
            labels_map={1: roi},
            dose_rate=True,
            cascade=cascade,
        )
        rate_gy_s = results[roi]["mean_dose_rate_gyh"] / 3600.0
        return gate_integrated_dose_gy(rate_gy_s, teff_h)

    if method == "hanscheid2017":
        activity, labels, _ = organ_proxy(bench, float(bench["time_from_injection_h"]))
        results = compute_roi_dose(
            image=activity,
            labels=labels,
            method=method,
            isotope=isotope,
            labels_map={1: roi},
            teff_params={roi: teff_h},
        )
        return results[roi]["mean_dose_gy"]

    if method == "hanscheid2018":
        # Teff = 2 t: acquisition at t = Teff / 2
        activity, labels, ct = organ_proxy(bench, teff_h / 2.0)
        results = compute_roi_dose(
            image=activity,
            labels=labels,
            method=method,
            isotope=isotope,
            ct_image=ct,
            labels_map={1: roi},
            method_params={"phantom": phantom},
        )
        return results[roi]["mean_dose_gy"]

    activity, labels, ct = organ_proxy(bench, float(bench["time_from_injection_h"]))
    results = compute_roi_dose(
        image=activity,
        labels=labels,
        method=method,
        isotope=isotope,
        ct_image=ct,
        labels_map={1: roi},
        teff_params={roi: teff_h},
        method_params={"phantom": phantom},
    )
    return results[roi]["mean_dose_gy"]


def radot_cross_dose_gy(bench: dict, isotope: str, method: str) -> float:
    cross = bench["cross_dose"]
    teff_h = float(bench["teff_h"])
    phantom = bench.get("phantom", "male")
    od_source = resolve_organ_name(cross.get("source_roi_name", bench["roi_name"]), phantom)
    od_target = resolve_organ_name(cross["target_roi_name"], phantom)

    s_cross = get_s_value(phantom, isotope, od_source, od_target)
    mass_scaling = float(cross["target_mass_g"]) / get_reference_mass(phantom, od_target)

    # Two voxels: the source (all the activity) and the target
    activity = np.array([float(bench["total_activity_bq"]), 0.0])
    source_mask = np.array([1, 0], dtype=np.uint8)
    target_mask = np.array([0, 1], dtype=np.uint8)

    if method == "hanscheid2018":
        return dose_hanscheid2018(
            spect_bq_a=activity,
            roi_a=target_mask,
            time_from_injection_h=teff_h / 2.0,
            s_value=s_cross,
            mass_scaling=mass_scaling,
            source_roi_a=source_mask,
        )
    return dose_madsen2018(
        spect_bq_a=activity,
        roi_a=target_mask,
        time_from_injection_h=float(bench["time_from_injection_h"]),
        s_value=s_cross,
        mass_scaling=mass_scaling,
        effective_time_h=teff_h,
        source_roi_a=source_mask,
    )


def evaluate(case: Case, bench_dir: Path = BENCH_DIR) -> dict:
    """Compute radot and OpenGATE doses for one case. Doses are in mGy per MBq."""
    bench = load_benchmark(case.isotope, bench_dir)
    mbq = float(bench["total_activity_bq"]) / 1e6
    teff_h = float(bench["teff_h"])

    if case.cross:
        gate = bench["cross_dose"]["opengate"]
        radot_gy = radot_cross_dose_gy(bench, case.isotope, case.method)
    else:
        gate = bench["opengate"]
        radot_gy = radot_self_dose_gy(bench, case.isotope, case.method, case.cascade)
    gate_gy = gate_integrated_dose_gy(float(gate["mean_dose_rate_gy_s"]), teff_h)

    diff_pct = (radot_gy - gate_gy) / gate_gy * 100.0
    passed = None
    if case.accept is not None:
        lo, hi = case.accept
        passed = lo <= abs(diff_pct) < hi

    return {
        "case": case,
        "bench": bench,
        "teff_h": teff_h,
        "radot": radot_gy * 1000.0 / mbq,
        "gate": gate_gy * 1000.0 / mbq,
        "gate_unc": float(gate["uncertainty_estimate_fraction"]) * 100.0,
        "diff": diff_pct,
        "passed": passed,
    }
