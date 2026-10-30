import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from benchmarks import (
    BENCH_DIR,
    BENCHMARK_FILES,
    CASES,
    CROSS_TARGETS,
    Case,
    cross_dose,
    is_distant,
    download_data,
    evaluate,
    gate_from_images,
    gate_integrated_dose_gy,
    gate_output_dir,
    load_benchmark,
    organ_images,
    organ_mask,
    run_opengate,
    SIM_DIR,
)

CHECKED_CASES = [c for c in CASES if c.accept is not None]
# The source organ ("": self-dose), then the cross-dose targets
ORGANS = [""] + list(CROSS_TARGETS)


def stored_gate(bench: dict, organ: str) -> dict:
    return cross_dose(bench, organ)["opengate"] if organ else bench["opengate"]


@pytest.fixture(scope="session")
def data_dir(request):
    """None for the default level; the downloaded images for --opengate=images."""
    if request.config.getoption("--opengate") == "json":
        return None
    return download_data()


@pytest.fixture(scope="session")
def gate_dir(request, data_dir):
    """None (published OpenGATE outputs), or new simulations for --opengate=simulate."""
    if request.config.getoption("--opengate") != "simulate":
        return None
    for isotope in BENCHMARK_FILES:
        run_opengate(isotope, data_dir, SIM_DIR)
    return SIM_DIR


@pytest.mark.parametrize("case", CHECKED_CASES, ids=str)
def test_radot_vs_opengate(case, data_dir, gate_dir):
    row = evaluate(case, data_dir=data_dir, gate_dir=gate_dir)
    lo, hi = case.accept
    assert row["radot"] > 0.0 and row["gate"] > 0.0
    assert row["passed"], (
        f"{case}: radot {row['radot']:.4g} vs GATE {row['gate']:.4g} mGy/MBq, "
        f"difference {row['diff']:+.2f} % (expected |diff| in [{lo}, {hi}) %)"
    )


def test_ac225_parent_only_and_cascade(data_dir, gate_dir):
    # Parent-only Madsen and LDM agree (both ignore the daughters) ...
    madsen = evaluate(Case("Ac225", "madsen2018", None), data_dir=data_dir, gate_dir=gate_dir)
    ldm = evaluate(Case("Ac225", "ldm", None), data_dir=data_dir, gate_dir=gate_dir)
    assert abs(madsen["radot"] - ldm["radot"]) / ldm["radot"] < 0.01

    # ... and GATE (full decay chain) is ~4.7 times higher than the parent alone
    ratio = ldm["gate"] / ldm["radot"]
    assert 4.5 < ratio < 5.0, f"GATE / LDM parent-only ratio: {ratio:.3f}"


@pytest.mark.parametrize("isotope", list(BENCHMARK_FILES))
def test_benchmark_json_matches_opengate_images(isotope, data_dir):
    # The stored OpenGATE values are those of the published images
    if data_dir is None:
        pytest.skip("needs --opengate=images")
    bench = load_benchmark(isotope)
    for organ in ORGANS:
        rate, _ = gate_from_images(isotope, organ_mask(data_dir, organ), data_dir)
        assert rate == pytest.approx(stored_gate(bench, organ)["mean_dose_rate_gy_s"], rel=1e-12), organ


@pytest.mark.parametrize("isotope", list(BENCHMARK_FILES))
def test_radot_with_opengate_dose_rate_map(isotope, data_dir):
    # radot --dose-rate-map with the OpenGATE map gives back the GATE dose, D = Ddot * Teff / ln2,
    # for the liver and the cross-dose organs (mass-weighted mean in radot, plain mean in GATE)
    if data_dir is None:
        pytest.skip("needs --opengate=images")
    import numpy as np
    import SimpleITK as sitk
    from radot.dose import compute_roi_dose

    bench = load_benchmark(isotope)
    activity, _, ct = organ_images(bench, 0.0, data_dir)
    rate_map = sitk.ReadImage(str(gate_output_dir(isotope, data_dir) / "output_dose.mhd"))
    rate_arr = sitk.GetArrayViewFromImage(rate_map)
    density = np.clip(1.0 + sitk.GetArrayViewFromImage(ct) / 1000.0, 1e-3, None)
    for organ in ORGANS:
        name = organ or "Liver"
        mask = organ_mask(data_dir, organ)
        labels = sitk.GetImageFromArray(mask.astype(np.uint16))
        labels.CopyInformation(activity)
        res = compute_roi_dose(
            image=activity,
            labels=labels,
            labels_map={1: name},
            ct_image=ct,
            method="madsen2018",
            isotope=isotope,
            method_params={name: {"teff": bench["teff_h"]}},
            dose_rate_map=rate_map,
            dose_rate_map_activity_mbq=bench["total_activity_bq"] / 1e6,
        )[name]
        weighted = float(np.sum(rate_arr[mask] * density[mask]) / np.sum(density[mask]))
        assert res["mean_dose_gy"] == pytest.approx(gate_integrated_dose_gy(weighted, bench["teff_h"]), rel=1e-6)
        stored = stored_gate(bench, organ)
        # the stored GATE value is a plain voxel mean: close to the mass-weighted one where the dose is
        # smooth (liver, organs beyond the electron reach), not in an organ touching the liver, whose
        # electron dose is in the border voxels (Y90: stomach -12 % with its gas, heart -2 %)
        smooth = not organ or is_distant(bench, organ)
        if stored["uncertainty_estimate_fraction"] < 0.1 and smooth:
            gate = gate_integrated_dose_gy(stored["mean_dose_rate_gy_s"], bench["teff_h"])
            assert res["mean_dose_gy"] == pytest.approx(gate, rel=0.02), name


@pytest.mark.parametrize("organ", ORGANS, ids=lambda o: o or "Liver")
@pytest.mark.parametrize("isotope", list(BENCHMARK_FILES))
def test_new_simulation_matches_stored(isotope, organ, data_dir, gate_dir):
    # A new Monte Carlo run agrees with the stored one within 3 combined standard uncertainties
    if gate_dir is None:
        pytest.skip("needs --opengate=simulate")
    bench = load_benchmark(isotope)
    stored = stored_gate(bench, organ)
    name = organ or "Liver"
    if stored["uncertainty_estimate_fraction"] >= 0.5:
        # e.g. Y90 cross-doses: almost no photons, the stored value is noise
        pytest.skip(f"{isotope} {name}: stored uncertainty {stored['uncertainty_estimate_fraction']:.0%}, not meaningful")
    # the new simulation may have another activity: compared at the benchmark activity
    rate, unc = gate_from_images(
        isotope, organ_mask(data_dir, organ), data_dir, gate_dir, is_distant(bench, organ),
        activity_bq=float(bench["total_activity_bq"]),
    )
    ref = stored["mean_dose_rate_gy_s"]
    sigma = (unc**2 + stored["uncertainty_estimate_fraction"] ** 2) ** 0.5
    diff = abs(rate - ref) / ref
    assert diff < 3 * sigma, f"{isotope} {name}: new {rate:.4e} vs stored {ref:.4e} Gy/s, {diff:.2%} > 3 x {sigma:.2%}"


@pytest.mark.parametrize("organ", ORGANS, ids=lambda o: o or "Liver")
def test_uncertainty_estimate_matches_batches(organ):
    # The organ uncertainty estimated from the voxel maps (benchmark JSON) agrees with the one
    # measured with 20 independent simulations (batch_uncertainty.py), within 2 sigma of the measurement
    bench = load_benchmark("Lu177")
    batches = json.loads((BENCH_DIR / "p142_lu177_batch_uncertainty.json").read_text())
    measured = batches["organs"][organ or bench["roi_name"]]
    # the uncertainty scales as 1/sqrt(activity): the estimate at the activity of one batch
    scale = np.sqrt(float(bench["total_activity_bq"]) / batches["activity_bq_per_batch"])
    estimate = stored_gate(bench, organ)["uncertainty_estimate_fraction"] * scale
    m, precision = measured["measured_uncertainty_fraction"], measured["measured_uncertainty_precision_fraction"]
    assert abs(estimate - m) < 2 * precision, (
        f"Lu177 {organ or 'Liver'}: estimated {estimate:.3%}, measured {m:.3%} +/- {precision:.3%}"
    )


def test_validation_page_is_up_to_date():
    # docs/7_validation_opengate.md is generated: python validation/opengate/report.py
    report = Path(__file__).resolve().parent / "report.py"
    result = subprocess.run(
        [sys.executable, str(report), "--check"], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
