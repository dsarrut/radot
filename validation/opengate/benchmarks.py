"""
Comparison of radot organ doses with OpenGATE Monte Carlo benchmarks.

Single source of truth for the validation: the cases and acceptance ranges below
are used by both test_opengate.py and report.py (docs/7_validation_opengate.md).

A benchmark JSON (see p142/p142_*_opengate_benchmark.json) gives the mean
OpenGATE dose rate at t=0 in a source organ (self-dose) and in one target organ
(cross-dose), for a uniform activity in the source organ. The OpenGATE dose rate
is integrated with the same mono-exponential clearance (Teff) as radot, so the
comparison tests the absorbed dose per decay, not the kinetics.

Three levels:
- default (data_dir=None): the organ is modelled as a single voxel with the
  benchmark volume, mass and total activity, and the OpenGATE values are read from
  the benchmark JSON. No image is needed; on p142 the differences with GATE are
  the same as with the full images within 0.02 percentage points.
- images (data_dir given): radot runs on the real CT, liver and organ labels, and
  the OpenGATE dose rate is averaged over the same masks in the OpenGATE images.
  The data are downloaded from Zenodo with download_data().
- simulate (gate_dir given): as images, but with OpenGATE outputs from a new
  simulation (run_opengate) instead of the published ones.
"""

import hashlib
import importlib.metadata
import json
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.parse
import urllib.request
import zipfile
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import SimpleITK as sitk

from radot.db import get_reference_mass, get_s_value, resolve_organ_name
from radot.dose import (
    compute_roi_dose,
    dose_hanscheid2018,
    dose_madsen2018,
    roi_dose_uncertainty_estimate,
)

BENCH_DIR = Path(__file__).resolve().parent / "p142"
DATA_DIR = BENCH_DIR / "data"
SIM_DIR = DATA_DIR / "simulated"

# Images and OpenGATE outputs, Zenodo doi:10.5281/zenodo.23696023 (version 1.2, doi:10.5281/zenodo.23809698)
DATA_URL = "https://zenodo.org/records/23809698/files/radot_opengate_p142_v1.2.zip?download=1"
DATA_SHA256 = "9408522b321147c7aa4a8cd246b50caab4379a2d4f35184fd6bba7d7860ede59"
DATA_TOP_FOLDER = "radot_opengate_p142_v1.2"
CT_FILE = "ct_4mm_cropped.nii.gz"
LIVER_FILE = "liver_4mm_cropped.nii.gz"
LABELS_FILE = "labels_all_4mm_cropped.nii.gz"
LABELS_NAMES_FILE = "labels_all_4mm_cropped.json"

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
    # Cross-dose target organ (key of CROSS_TARGETS); "" = self-dose in the source organ
    target: str = ""

    @property
    def cross(self) -> bool:
        return bool(self.target)

    @property
    def label(self) -> str:
        label = METHOD_NAMES[self.method]
        return label + " (cascade)" if self.cascade else label

    def __str__(self) -> str:
        kind = self.target if self.cross else "self"
        cascade = "-cascade" if self.cascade else ""
        return f"{self.isotope}-{self.method}{cascade}-{kind}"


# Cross-dose targets: patient label (TG-263) -> OpenDose target organ(s). Whole organs entirely
# inside the cropped CT. The stomach and the heart are segmented whole in the patient: OpenDose
# wall + contents, with the mass-weighted S-value.
CROSS_TARGETS = {
    "Kidney_R": ["Right kidney C+M+P"],
    "Kidney_L": ["Left kidney C+M+P"],
    "Spleen": ["Spleen"],
    "Pancreas": ["Pancreas"],
    "Adrenal_R": ["Adrenal, right"],
    "Adrenal_L": ["Adrenal, left"],
    "Stomach": ["Stomach wall", "Stomach contents"],
    "Heart": ["Heart wall", "Heart contents (blood)"],
}

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
    # Cross-dose (organ <- liver): photon dose, sensitive to phantom vs patient anatomy.
    # Only the right kidney is checked; the other organs are shown.
    Case("Lu177", "madsen2018", (0, 50), target="Kidney_R"),
    Case("Lu177", "hanscheid2018", (0, 50), H2018_NOTE, target="Kidney_R"),
    *[Case("Lu177", "madsen2018", None, target=t) for t in CROSS_TARGETS if t != "Kidney_R"],
    *[Case("I131", "madsen2018", None, target=t) for t in CROSS_TARGETS],
    *[Case("Y90", "madsen2018", None, target=t) for t in CROSS_TARGETS],
]


def load_benchmark(isotope: str, bench_dir: Path = BENCH_DIR) -> dict:
    with open(Path(bench_dir) / BENCHMARK_FILES[isotope]) as f:
        return json.load(f)


# Beyond this distance from the source, the electrons do not reach the organ (Y90 maximum beta range
# ~11 mm): the dose of the electron run there is bremsstrahlung only (see gate_roi_uncertainty)
ELECTRON_REACH_MM = 10.0


def cross_dose(bench: dict, target: str) -> dict:
    """Entry of the benchmark cross_doses list for this target organ."""
    for entry in bench["cross_doses"]:
        if entry["target_roi_name"] == target:
            return entry
    raise KeyError(f"no cross-dose to {target} in the {bench['radionuclide']} benchmark")


def is_distant(bench: dict, target: str) -> bool:
    """The organ is beyond the electron reach of the source (always False for the source itself)."""
    return bool(target) and cross_dose(bench, target)["min_distance_to_source_mm"] > ELECTRON_REACH_MM


def gate_output_dir(isotope: str, data_dir: Path, bench_dir: Path = BENCH_DIR) -> Path:
    """OpenGATE output folder of a radionuclide, as given in its opengate_dose_rate config."""
    with open(Path(bench_dir) / f"p142_dose_rate_{isotope.lower()}.json") as f:
        return Path(data_dir) / json.load(f)["output_folder"]


def has_data(data_dir: Path = DATA_DIR, bench_dir: Path = BENCH_DIR) -> bool:
    data_dir = Path(data_dir)
    files = [data_dir / f for f in (CT_FILE, LIVER_FILE, LABELS_FILE, LABELS_NAMES_FILE)]
    for isotope in BENCHMARK_FILES:
        out = gate_output_dir(isotope, data_dir, bench_dir)
        folders = [out]
        with open(Path(bench_dir) / f"p142_dose_rate_{isotope.lower()}.json") as f:
            if json.load(f)["mode"] == "vrt":
                # since version 1.1: the electron and photon parts of the VRT simulations;
                # since version 1.2 (bremsstrahlung): reference.json
                folders += [out / "e", out / "gamma_tle"]
                files.append(out / "reference.json")
        files += [d / name for d in folders for name in ("output_dose.mhd", "output_dose_uncertainty.mhd")]
    return all(f.exists() for f in files)


def download_data(data_dir: Path = DATA_DIR, force: bool = False) -> Path:
    """Download the images and OpenGATE outputs from Zenodo into data_dir (if needed)."""
    data_dir = Path(data_dir)
    if has_data(data_dir) and not force:
        return data_dir
    with tempfile.TemporaryDirectory() as tmp:
        zip_path = Path(tmp) / "data.zip"
        urllib.request.urlretrieve(DATA_URL, zip_path)
        sha256 = hashlib.sha256(zip_path.read_bytes()).hexdigest()
        if sha256 != DATA_SHA256:
            raise RuntimeError(f"Checksum mismatch for {DATA_URL}: {sha256}")
        with zipfile.ZipFile(zip_path) as z:
            for info in z.infolist():
                if info.is_dir():
                    continue
                rel = Path(info.filename).relative_to(DATA_TOP_FOLDER)
                if ".." in rel.parts:
                    raise RuntimeError(f"Unexpected path in archive: {info.filename}")
                dest = data_dir / rel
                dest.parent.mkdir(parents=True, exist_ok=True)
                with z.open(info) as src, open(dest, "wb") as dst:
                    shutil.copyfileobj(src, dst)
    return data_dir


@lru_cache(maxsize=None)
def _read_images(data_dir: str):
    data_dir = Path(data_dir)
    liver = sitk.ReadImage(str(data_dir / LIVER_FILE))
    ct = sitk.ReadImage(str(data_dir / CT_FILE))
    labels = sitk.ReadImage(str(data_dir / LABELS_FILE))
    with open(data_dir / LABELS_NAMES_FILE) as f:
        label_names = {name: int(i) for i, name in json.load(f).items()}
    return liver, ct, labels, label_names


def organ_images(bench: dict, time_from_injection_h: float, data_dir: Path):
    """Real activity (uniform in the liver), label and CT images."""
    liver, ct, _, _ = _read_images(str(data_dir))
    mask = liver > 0
    num_voxels = int(sitk.GetArrayViewFromImage(mask).sum())
    activity = sitk.Cast(mask, sitk.sitkFloat32) * (
        float(bench["total_activity_bq"]) / num_voxels
    )
    activity.SetMetaData("unit", "Bq")
    activity.SetMetaData("time_from_injection", str(time_from_injection_h))
    labels = sitk.Cast(mask, sitk.sitkUInt16)
    return activity, labels, ct


def organ_mask(data_dir: Path, target: str = "") -> np.ndarray:
    """Mask of a target organ (patient label name), or of the source (liver activity) by default."""
    liver, _, labels, label_names = _read_images(str(data_dir))
    if not target:
        return sitk.GetArrayFromImage(liver) > 0
    return sitk.GetArrayFromImage(labels) == label_names[target]


def opengate_version() -> dict:
    """Installed opengate version, and its git commit when installed from a checkout."""
    # Use the install record, not an import: under pytest, "opengate" can resolve
    # to this folder (validation/opengate) instead of the OpenGATE package.
    dist = importlib.metadata.distribution("opengate")
    info = {"version": dist.version}
    direct_url = json.loads(dist.read_text("direct_url.json") or "{}")
    url = direct_url.get("url", "")
    if url.startswith("file://"):
        repo = Path(urllib.parse.urlparse(url).path)
        result = subprocess.run(
            ["git", "-C", str(repo), "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            info["commit"] = result.stdout.strip()
            status = subprocess.run(
                ["git", "-C", str(repo), "status", "--porcelain", "--untracked-files=no"],
                capture_output=True,
                text=True,
            )
            info["dirty"] = bool(status.stdout.strip())
        elif (repo / ".git").exists():
            info["commit_error"] = result.stderr.strip()
    return info


def run_opengate(
    isotope: str,
    data_dir: Path = DATA_DIR,
    sim_dir: Path = SIM_DIR,
    bench_dir: Path = BENCH_DIR,
) -> Path:
    """
    Run the p142 OpenGATE simulation of this radionuclide, with the activity of its configuration.
    Inputs are read from data_dir; outputs go to sim_dir/<output_folder>, with the log and
    opengate_version.json. Returns the output folder.
    VRT radionuclides: p142/create_campaigns.py, as for the references (bremsstrahlung included,
    parallel jobs); analog (Ac225): opengate_dose_rate.
    """
    config = (Path(bench_dir) / f"p142_dose_rate_{isotope.lower()}.json").resolve()
    with open(config) as f:
        cfg = json.load(f)
    out = gate_output_dir(isotope, sim_dir, bench_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    if cfg.get("mode", "vrt") == "vrt":
        # a separate process: under pytest, "import opengate" resolves to this folder
        cmd = [
            sys.executable, str(Path(bench_dir) / "create_campaigns.py"),
            "--radionuclide", isotope, "--activity", str(float(cfg["activity_bq"])),
            "--campaigns", str(Path(sim_dir).resolve() / "campaigns"),
            "--run", "-j", str(os.cpu_count() or 4), "-o", str(out),
        ]
    else:
        exe = shutil.which("opengate_dose_rate")
        if exe is None:
            raise RuntimeError("opengate_dose_rate not found: install opengate (pip install opengate)")
        cmd = [exe, str(config), "--output_folder", str(out)]
    with open(out / "opengate.log", "w") as log:
        subprocess.run(cmd, cwd=data_dir, stdout=log, stderr=subprocess.STDOUT, check=True)
    with open(out / "opengate_version.json", "w") as f:
        json.dump(opengate_version(), f, indent=2)
    return out


def simulated_activity_bq(out: Path, isotope: str, bench_dir: Path = BENCH_DIR) -> float:
    """Activity simulated in an OpenGATE output folder: reference.json (create_campaigns.py), else the configuration."""
    reference = Path(out) / "reference.json"
    if reference.exists():
        return float(json.loads(reference.read_text())["activity_bq"])
    with open(Path(bench_dir) / f"p142_dose_rate_{isotope.lower()}.json") as f:
        return float(json.load(f)["activity_bq"])


def gate_from_images(
    isotope: str,
    mask: np.ndarray,
    data_dir: Path,
    gate_dir: Optional[Path] = None,
    distant: bool = False,
    activity_bq: Optional[float] = None,
) -> Tuple[float, float]:
    """
    Mean OpenGATE dose rate (Gy/s) in the mask, and its relative uncertainty.
    The OpenGATE outputs are read from gate_dir (default: data_dir).
    distant: the organ is beyond the electron reach of the source (see gate_roi_uncertainty).
    activity_bq: scale the dose rate to this activity (e.g. the one of the benchmark JSON);
    default: the activity of the simulation.
    """
    out = gate_output_dir(isotope, gate_dir if gate_dir is not None else data_dir)
    dose = sitk.ReadImage(str(out / "output_dose.mhd"))
    rate = float(sitk.GetArrayViewFromImage(dose)[mask].mean())
    if activity_bq is not None:
        rate *= activity_bq / simulated_activity_bq(out, isotope)
    return rate, gate_roi_uncertainty(out, mask, isotope, distant)


def gate_roi_uncertainty(out: Path, mask: np.ndarray, isotope: str, distant: bool = False) -> float:
    """
    Relative uncertainty of the mean OpenGATE dose rate in the mask.

    In VRT mode, the electrons (e/) and the photons (gamma_tle/) are two independent
    simulations whose dose maps add up to the merged one: each part is estimated with
    its own correlation model (mean beta range, ROI mean chord) and the variances are
    added. The electron run also contains the bremsstrahlung photons (opengate_dose_rate option):
    in an organ beyond the electron reach (distant=True, see ELECTRON_REACH_MM), its dose is
    bremsstrahlung only, so the photon model is used there too; in an organ touching the
    source, the electrons crossing the boundary dominate (electron model). Otherwise (analog mode, or parts missing), the merged map is
    used with the electron model, which underestimates the uncertainty of a photon dose.
    Checked on p142 Lu177 against 20 independent simulations (batch_uncertainty.py).
    """
    out = Path(out)
    parts = [(out / "e", distant), (out / "gamma_tle", True)]
    if not all((folder / "output_dose_uncertainty.mhd").exists() for folder, _ in parts):
        dose = sitk.ReadImage(str(out / "output_dose.mhd"))
        unc = sitk.ReadImage(str(out / "output_dose_uncertainty.mhd"))
        return roi_dose_uncertainty_estimate(dose, unc, mask, isotope=isotope)
    variance, total = 0.0, 0.0
    for folder, photons in parts:
        dose = sitk.ReadImage(str(folder / "output_dose.mhd"))
        unc = sitk.ReadImage(str(folder / "output_dose_uncertainty.mhd"))
        mean = float(sitk.GetArrayViewFromImage(dose)[mask].mean())
        total += mean
        if mean > 0:
            u = roi_dose_uncertainty_estimate(dose, unc, mask, isotope=isotope, photons=photons)
            variance += (u * mean) ** 2
    if total <= 0:
        # no dose at all (e.g. Y90 in a distant organ): no meaningful value
        return 1.0
    return float(np.sqrt(variance) / total)


def gate_integrated_dose_gy(dose_rate_gy_s: float, teff_h: float) -> float:
    """Integrate a dose rate at t=0 with mono-exponential clearance of half-life teff_h."""
    k = np.log(2) / (teff_h * 3600.0)
    return dose_rate_gy_s / k


def organ_single_voxel(bench: dict, time_from_injection_h: float):
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


def radot_self_dose_gy(
    bench: dict,
    isotope: str,
    method: str,
    cascade: bool = False,
    data_dir: Optional[Path] = None,
) -> float:
    def organ_inputs(bench, time_from_injection_h):
        if data_dir is None:
            return organ_single_voxel(bench, time_from_injection_h)
        return organ_images(bench, time_from_injection_h, data_dir)

    teff_h = float(bench["teff_h"])
    roi = bench["roi_name"]
    phantom = bench.get("phantom", "male")

    if method == "ldm":
        activity, labels, ct = organ_inputs(bench, 0.0)
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
        activity, labels, _ = organ_inputs(bench, float(bench["time_from_injection_h"]))
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
        activity, labels, ct = organ_inputs(bench, teff_h / 2.0)
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

    activity, labels, ct = organ_inputs(bench, float(bench["time_from_injection_h"]))
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


def radot_cross_dose_gy(
    bench: dict, isotope: str, method: str, target: str, data_dir: Optional[Path] = None
) -> float:
    cross = cross_dose(bench, target)
    teff_h = float(bench["teff_h"])
    phantom = bench.get("phantom", "male")
    od_source = resolve_organ_name(cross.get("source_roi_name", bench["roi_name"]), phantom)

    # Several OpenDose targets (e.g. wall + contents): mass-weighted S-value
    od_masses = {t: get_reference_mass(phantom, t) for t in CROSS_TARGETS[target]}
    ref_mass_g = sum(od_masses.values())
    s_cross = sum(m * get_s_value(phantom, isotope, od_source, t) for t, m in od_masses.items()) / ref_mass_g

    if data_dir is None:
        # Two voxels: the source (all the activity) and the target
        target_mass_g = float(cross["target_mass_g"])
        activity = np.array([float(bench["total_activity_bq"]), 0.0])
        source_mask = np.array([1, 0], dtype=np.uint8)
        target_mask = np.array([0, 1], dtype=np.uint8)
    else:
        source, target_arr = organ_mask(data_dir), organ_mask(data_dir, target)
        activity_img, _, ct = organ_images(bench, 0.0, data_dir)
        activity = sitk.GetArrayFromImage(activity_img)
        source_mask = source.astype(np.uint8)
        target_mask = target_arr.astype(np.uint8)
        spacing = ct.GetSpacing()
        voxel_volume_ml = spacing[0] * spacing[1] * spacing[2] / 1000.0
        ct_arr = sitk.GetArrayViewFromImage(ct)
        density = np.maximum(0.0, 1.0 + ct_arr[target_arr] / 1000.0)
        target_mass_g = float(density.sum() * voxel_volume_ml)

    mass_scaling = target_mass_g / ref_mass_g

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


def evaluate(
    case: Case,
    bench_dir: Path = BENCH_DIR,
    data_dir: Optional[Path] = None,
    gate_dir: Optional[Path] = None,
) -> dict:
    """
    Compute radot and OpenGATE doses for one case. Doses are in mGy per MBq.
    With data_dir, use the real images and the OpenGATE images (see download_data);
    with gate_dir, read the OpenGATE images from there (see run_opengate).
    """
    bench = load_benchmark(case.isotope, bench_dir)
    mbq = float(bench["total_activity_bq"]) / 1e6
    teff_h = float(bench["teff_h"])

    if case.cross:
        radot_gy = radot_cross_dose_gy(bench, case.isotope, case.method, case.target, data_dir)
        gate = cross_dose(bench, case.target)["opengate"]
    else:
        radot_gy = radot_self_dose_gy(
            bench, case.isotope, case.method, case.cascade, data_dir
        )
        gate = bench["opengate"]
    gate_rate = float(gate["mean_dose_rate_gy_s"])
    gate_unc = float(gate["uncertainty_estimate_fraction"])
    if data_dir is not None:
        gate_rate, gate_unc = gate_from_images(
            case.isotope, organ_mask(data_dir, case.target), data_dir, gate_dir,
            is_distant(bench, case.target), activity_bq=float(bench["total_activity_bq"]),
        )
    gate_gy = gate_integrated_dose_gy(gate_rate, teff_h)

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
        "gate_unc": gate_unc * 100.0,
        "diff": diff_pct,
        "passed": passed,
    }
