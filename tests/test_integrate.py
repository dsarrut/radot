"""radot integrate, ROI level (fit first): models, error models, AICc selection, ROI curves, CLI."""
import json
import logging
import math

import numpy as np
import pytest
import SimpleITK as sitk
from click.testing import CliRunner

from radot.cli import main
from radot.db import get_physical_half_life_h
from radot.integrate import (
    TAC_MODELS,
    akaike_weights,
    fit_tac,
    integrate_rois,
    integrate_tac,
    parse_error_model,
    read_integrate,
    roi_time_activity_curves,
)

LP = math.log(2) / get_physical_half_life_h("Lu177")
T6 = np.array([2.0, 4.0, 24.0, 48.0, 96.0, 168.0])
T7 = np.array([1.0, 4.0, 24.0, 48.0, 96.0, 168.0, 240.0])

# (model, parameters, analytic TIA in MBq.h)
TRUTH = {
    "f1": ([80.0], 80.0 / LP),
    "f2": ([100.0, 0.02], 100.0 / (0.02 + LP)),
    "f3": ([60.0, 20.0, 0.05], 60.0 / (0.05 + LP) + 20.0 / LP),
    "f4": ([60.0, 20.0, 0.2, 0.01], 60.0 / (0.2 + LP) + 20.0 / (0.01 + LP)),
    "f4a": ([50.0, 0.01, 0.3], 50.0 / (0.01 + LP) - 50.0 / (0.3 + LP)),
    "f5": ([40.0, 30.0, 10.0, 0.3, 0.02], 40.0 / (0.3 + LP) + 30.0 / (0.02 + LP) + 10.0 / LP),
    "f6": ([50.0, 30.0, 5.0, 0.5, 0.05, 0.01], 50.0 / (0.5 + LP) + 30.0 / (0.05 + LP) + 5.0 / (0.01 + LP)),
    "f6a": ([40.0, 20.0, 0.01, 0.05, 1.0],
            40.0 / (0.01 + LP) + 20.0 / (0.05 + LP) - 60.0 / (1.0 + LP)),
}


def _curve(model, p, t):
    m = TAC_MODELS[model]
    amps, rates = p[: m.n_amplitudes], p[m.n_amplitudes:]
    y = np.zeros_like(t)
    for coefs, r in m.terms:
        y += np.dot(coefs, amps) * np.exp(-((rates[r] if r is not None else 0.0) + LP) * t)
    return y


# ----------------------------------------------------------------------------------------------
# models
# ----------------------------------------------------------------------------------------------


@pytest.mark.parametrize("model", list(TRUTH))
def test_exact_curve_recovers_the_integral(model):
    p, tia = TRUTH[model]
    y = _curve(model, p, T7)
    res = fit_tac(model, T7, y, LP, parse_error_model("fsd:0.05"))
    assert res.tia_mbq_h == pytest.approx(tia, rel=1e-5)


def test_noisy_curve_se_and_cv():
    p, tia = TRUTH["f2"]
    y = _curve("f2", p, T6) * (1 + 0.05 * np.random.default_rng(3).standard_normal(len(T6)))
    res = fit_tac("f2", T6, y, LP, parse_error_model("proportional"))
    assert res.tia_se_mbq_h > 0 and abs(res.tia_mbq_h - tia) < 4 * res.tia_se_mbq_h
    assert set(res.params_cv) == {"A1", "lambda1"} and res.aicc is not None
    assert res.n_estimated == 3  # K + 1: the variance scale is estimated


def test_number_of_points_rules():
    y = _curve("f2", TRUTH["f2"][0], T6[:3])
    # J = K + 1 = 3 <= N = 3: fitted, but no SE (J + 1 > N) and no AICc (J + 2 > N)
    res = fit_tac("f2", T6[:3], y, LP, parse_error_model("proportional"))
    assert res.params_se is None and res.aicc is None
    assert any("SE not computable" in i for i in res.issues)
    # absolute weighting: J = K = 2, SE computable
    assert fit_tac("f2", T6[:3], y, LP, parse_error_model("fsd:0.1")).params_se is not None
    with pytest.raises(ValueError, match="J <= N"):
        fit_tac("f4a", T6[:3], y, LP, parse_error_model("proportional"))
    with pytest.raises(ValueError, match="J \\+ 2 <= N"):
        integrate_tac(T6[:5], _curve("f2", TRUTH["f2"][0], T6[:5]), LP, "auto", "proportional",
                      candidates=["f2", "f4a"], select="best")


def test_rate_at_bound_reported():
    y = _curve("f1", [80.0], T6)
    res = fit_tac("f2", T6, y, LP, parse_error_model("fsd:0.05"))
    assert res.at_bound == ["lambda1"] and any("bound" in i for i in res.issues)
    assert res.teff_terminal_h == pytest.approx(get_physical_half_life_h("Lu177"), rel=1e-6)


def test_non_physical_fit_is_an_error():
    # negative late values (noise around 0): f3 fits a negative retained fraction, the curve goes < 0
    y = np.array([50.0, 40.0, 10.0, 2.0, -1.0, -2.0])
    with pytest.raises(ValueError, match="non-physical"):
        fit_tac("f3", T6, y, LP, parse_error_model("sd:1"))


@pytest.mark.parametrize("value, err", [(None, "Missing error model"), ("", "Missing error model"),
                                        ("relative", "Unknown error model"), ("fsd:0", "X > 0"),
                                        ("sd:abc", "X > 0")])
def test_error_model_no_default(value, err):
    with pytest.raises(ValueError, match=err):
        parse_error_model(value)


def test_error_model_kinds():
    assert parse_error_model("Proportional").relative
    em = parse_error_model("fsd:0.1")
    assert not em.relative and em.theta(np.array([10.0]))[0] == pytest.approx(1.0)
    assert parse_error_model("sd:2").theta(np.array([10.0, 50.0])).tolist() == [2.0, 2.0]
    with pytest.raises(ValueError, match="> 0 at every time point"):
        parse_error_model("proportional").theta(np.array([1.0, 0.0]))


# ----------------------------------------------------------------------------------------------
# model selection
# ----------------------------------------------------------------------------------------------


def test_auto_best_and_average():
    p, tia = TRUTH["f4a"]
    y = _curve("f4a", p, T7) * (1 + 0.03 * np.random.default_rng(0).standard_normal(len(T7)))
    best = integrate_tac(T7, y, LP, "auto", "proportional", candidates=["f2", "f4a"], select="best")
    avg = integrate_tac(T7, y, LP, "auto", "proportional", candidates=["f2", "f4a"], select="average")
    weights = {c: f["weight"] for c, f in avg["fits"].items() if f["weight"] is not None}
    assert sum(weights.values()) == pytest.approx(1.0)
    assert best["model"] == max(weights, key=weights.get)
    used = {c: w for c, w in weights.items() if w > 0.01}
    expected = sum(w * avg["fits"][c]["tia_mbq_h"] for c, w in used.items()) / sum(used.values())
    assert avg["model"] == "average" and avg["tia_mbq_h"] == pytest.approx(expected)


def test_akaike_weights():
    w = akaike_weights([10.0, 12.0])
    assert w.sum() == pytest.approx(1.0) and w[0] / w[1] == pytest.approx(math.e)


def test_auto_needs_candidates_and_select():
    y = _curve("f2", TRUTH["f2"][0], T7)
    with pytest.raises(ValueError, match="candidate"):
        integrate_tac(T7, y, LP, "auto", "proportional", select="best")
    with pytest.raises(ValueError, match="select"):
        integrate_tac(T7, y, LP, "auto", "proportional", candidates=["f2", "f3"])
    with pytest.raises(ValueError, match="only used with model 'auto'"):
        integrate_tac(T7, y, LP, "f2", "proportional", select="best")
    with pytest.raises(ValueError, match="Unknown model"):
        integrate_tac(T7, y, LP, "biexp", "proportional")


# ----------------------------------------------------------------------------------------------
# ROI level, images
# ----------------------------------------------------------------------------------------------

LIVER = (slice(2, 5), slice(2, 5), slice(2, 5))
SPLEEN = (slice(6, 8), slice(6, 8), slice(6, 8))
TIMES = [4.0, 24.0, 96.0, 168.0]


def _write(path, arr, sidecar):
    img = sitk.GetImageFromArray(arr)
    img.SetSpacing((2.0, 2.0, 2.0))
    sitk.WriteImage(img, str(path))
    if sidecar is not None:
        path.with_name(path.name.split(".")[0] + ".json").write_text(json.dumps(sidecar))
    return path


def _study(tmp_path, decay_correction="start", spleen_at=(0, 1, 2, 3), drop=None, unit="bqml", sex=None):
    """Liver: f2 (A = 100 MBq, lambda1 = 0.02/h); spleen: f2 (20 MBq, 0.05/h). bqml, 2 mm voxels."""
    lab = np.zeros((10, 10, 10), np.uint8)
    lab[LIVER] = 1
    lab[SPLEEN] = 2
    labels = _write(tmp_path / "labels.nii.gz", lab, {"1": "Liver", "2": "Spleen"})
    images = []
    for i, t in enumerate(TIMES):
        arr = np.zeros((10, 10, 10), np.float32)
        liver = 100.0 * math.exp(-(0.02 + LP) * t)
        spleen = 20.0 * math.exp(-(0.05 + LP) * t)
        if decay_correction == "admin":
            liver, spleen = liver * math.exp(LP * t), spleen * math.exp(LP * t)
        arr[LIVER] = liver * 1e6 / (27 * 0.008)  # Bq/mL, 27 voxels of 0.008 mL
        if i in spleen_at:
            arr[SPLEEN] = spleen * 1e6 / (8 * 0.008)
        sc = {"modality": "spect", "unit": unit, "radionuclide": "Lu177", "time_from_injection": t,
              "decay_correction": decay_correction}
        if sex:
            sc["sex"] = sex
        if drop:
            sc.pop(drop)
        images.append(_write(tmp_path / f"spect_tp{i + 1}.nii.gz", arr, sc))
    if spleen_at != (0, 1, 2, 3):  # the spleen label is absent at the other time points
        lbls = []
        for i in range(len(TIMES)):
            li = lab.copy()
            if i not in spleen_at:
                li[SPLEEN] = 0
            lbls.append(_write(tmp_path / f"labels_tp{i + 1}.nii.gz", li, {"1": "Liver", "2": "Spleen"}))
        return images, lbls
    return images, labels


@pytest.mark.parametrize("dc", ["start", "admin"])
def test_roi_tia_from_images(tmp_path, dc):
    images, labels = _study(tmp_path, decay_correction=dc)
    doc = integrate_rois(images, labels, model="f2", error_model="proportional")
    assert doc["radot_integrate"] == 1 and doc["radionuclide"] == "Lu177"
    assert doc["decay_correction"] == [dc] * 4
    liver = doc["rois"]["Liver"]
    assert liver["activities_mbq"][0] == pytest.approx(100.0 * math.exp(-(0.02 + LP) * 4.0), rel=1e-5)
    assert liver["tia_mbq_h"] == pytest.approx(100.0 / (0.02 + LP), rel=1e-5)
    assert liver["tia_mbq_s"] == pytest.approx(liver["tia_mbq_h"] * 3600)
    assert doc["rois"]["Spleen"]["tia_mbq_h"] == pytest.approx(20.0 / (0.05 + LP), rel=1e-5)


def test_roi_sorted_by_time(tmp_path):
    images, labels = _study(tmp_path)
    tac = roi_time_activity_curves(images[::-1], labels)
    assert tac["times_h"] == TIMES and tac["image_names"][0] == "spect_tp1.nii.gz"


# AGENTS rule 5: values only in the sidecar work (above); missing everywhere is an explicit error
@pytest.mark.parametrize("key, err", [("decay_correction", "Decay correction"),
                                      ("time_from_injection", "Time point 1"),
                                      ("radionuclide", "Radionuclide not given"),
                                      ("unit", "unit None")])
def test_roi_missing_metadata_is_an_error(tmp_path, key, err):
    images, labels = _study(tmp_path, drop=key)
    with pytest.raises(ValueError, match=err):
        integrate_rois(images, labels, model="f2", error_model="proportional")


def test_isotope_option_vs_sidecar(tmp_path):
    images, labels = _study(tmp_path, drop="radionuclide")
    doc = integrate_rois(images, labels, model="f2", error_model="proportional", isotope="Lu177")
    assert doc["radionuclide"] == "Lu177"
    images, labels = _study(tmp_path)
    with pytest.raises(ValueError, match="Radionuclide mismatch"):
        integrate_rois(images, labels, model="f2", error_model="proportional", isotope="Y90")


def test_roi_unit_not_activity_is_an_error(tmp_path):
    images, labels = _study(tmp_path, unit="suv")
    with pytest.raises(ValueError, match="bq or bqml"):
        integrate_rois(images, labels, model="f2", error_model="proportional")


def test_roi_absent_at_a_time_point(tmp_path, caplog):
    images, labels = _study(tmp_path, spleen_at=(0, 1, 2))
    with caplog.at_level(logging.WARNING, logger="radot"):
        doc = integrate_rois(images, labels, model="f2", error_model="proportional")
    assert "Spleen" not in doc["rois"] and "Liver" in doc["rois"]
    assert "ROI 'Spleen' has no voxel at time point(s) 4" in caplog.text
    with pytest.raises(ValueError, match="needs it at every time point"):
        integrate_rois(images, labels, model="f2", error_model="proportional", rois=["Spleen"])


def test_model_choice_checked_once_for_all_rois(tmp_path):
    images, labels = _study(tmp_path)
    # the same for every ROI: one error, not one skipped ROI each (f6: J = 7 > N = 4)
    with pytest.raises(ValueError, match="Model f6: J = 7 .* N = 4"):
        integrate_rois(images, labels, model="f6", error_model="proportional")
    with pytest.raises(ValueError, match="Candidate model\\(s\\) f2, f4a"):
        integrate_rois(images, labels, model="auto", candidates=["f2", "f4a"], select="best",
                       error_model="proportional")


def test_roi_fit_failure_skipped_in_loop_error_when_requested(tmp_path, caplog, monkeypatch):
    from radot.integrate import roi as roi_module

    images, labels = _study(tmp_path)
    real = roi_module.integrate_curve

    def failing_for_spleen(times, acts, *a, **k):
        if acts[0] < 50:  # the spleen starts at ~16 MBq, the liver at ~91 MBq
            raise ValueError("Model f2: non-physical fit")
        return real(times, acts, *a, **k)

    monkeypatch.setattr(roi_module, "integrate_curve", failing_for_spleen)
    with caplog.at_level(logging.WARNING, logger="radot"):
        doc = integrate_rois(images, labels, model="f2", error_model="proportional")
    assert list(doc["rois"]) == ["Liver"] and "ROI 'Spleen' skipped: Model f2" in caplog.text
    with pytest.raises(ValueError, match="ROI 'Spleen': Model f2"):
        integrate_rois(images, labels, model="f2", error_model="proportional", rois=["Spleen"])


# ----------------------------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------------------------


def test_cli_integrate(tmp_path):
    images, labels = _study(tmp_path)
    out = tmp_path / "tia.json"
    args = ["integrate", *map(str, images), "-l", str(labels), "--model", "f2",
            "--error-model", "proportional", "-o", str(out)]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, result.output
    assert "Liver" in result.output
    doc = read_integrate(out)
    assert doc["rois"]["Liver"]["tia_mbq_h"] == pytest.approx(100.0 / (0.02 + LP), rel=1e-5)
    assert doc["error_model"] == "proportional" and doc["model"] == "f2" and doc["sex"] is None


def test_patient_sex_recorded(tmp_path):
    images, labels = _study(tmp_path, sex="female")
    assert integrate_rois(images, labels, model="f2", error_model="proportional")["sex"] == "female"


@pytest.mark.parametrize("drop, err", [("--error-model", "Missing error model"),
                                       ("--model", "Missing option '--model'"),
                                       ("-l", "Missing option '-l'")])
def test_cli_integrate_no_default(tmp_path, drop, err):
    images, labels = _study(tmp_path)
    opts = {"-l": str(labels), "--model": "f2", "--error-model": "proportional"}
    opts.pop(drop)
    args = ["integrate", *map(str, images)] + [x for kv in opts.items() for x in kv]
    result = CliRunner().invoke(main, args)
    assert result.exit_code != 0 and err in result.output


def test_cli_list_models():
    result = CliRunner().invoke(main, ["integrate", "--list-models"])
    assert result.exit_code == 0 and all(m in result.output for m in TAC_MODELS)


def test_read_integrate_rejects_other_files(tmp_path):
    p = tmp_path / "x.json"
    p.write_text(json.dumps({"radot_results": 1, "rois": {}}))
    with pytest.raises(ValueError, match="not a radot integrate file"):
        read_integrate(p)
