"""radot integrate without a fitted function: trapezoid (head / body / tail) and single time point models."""
import math

import numpy as np
import pytest
import SimpleITK as sitk
from click.testing import CliRunner

from radot.cli import main
from radot.db import get_physical_half_life_h
from radot.dose import compute_roi_dose, compute_roi_dose_tia_table
from radot.integrate import (
    check_integration_options,
    integrate_rois,
    integrate_stp,
    integrate_trapezoid,
    read_integrate,
)

from test_integrate import LP, TIMES, _study

T_PHYS = get_physical_half_life_h("Lu177")


# ----------------------------------------------------------------------------------------------
# trapezoid
# ----------------------------------------------------------------------------------------------


@pytest.mark.parametrize("head, expected", [("zero", 0.0), ("ramp", 10.0), ("const", 20.0)])
def test_trapezoid_parts(head, expected):
    r = integrate_trapezoid([2.0, 4.0], [10.0, 6.0], LP, head, "phys")
    parts = r["fits"]["trapezoid"]["parts_mbq_h"]
    assert parts == pytest.approx({"head": expected, "body": 16.0, "tail": 6.0 / LP})
    assert r["tia_mbq_h"] == pytest.approx(expected + 16.0 + 6.0 / LP)
    assert r["tia_se_mbq_h"] is None and r["model"] == "trapezoid"


def test_trapezoid_exact_parts_on_a_mono_exponential():
    lam = 0.05 + LP
    t = np.array([4.0, 24.0, 48.0, 96.0])
    a = 100.0 * np.exp(-lam * t)
    r = integrate_trapezoid(t, a, LP, "extrap", "fit", tail_points=3)
    parts = r["fits"]["trapezoid"]["parts_mbq_h"]
    # head: the mono-exponential through the first two points, back to 0: exact here
    assert parts["head"] == pytest.approx(100.0 * (1 - math.exp(-lam * 4.0)) / lam, rel=1e-9)
    assert parts["tail"] == pytest.approx(a[-1] / lam, rel=1e-9)
    assert r["fits"]["trapezoid"]["teff_terminal_h"] == pytest.approx(math.log(2) / lam)


def test_trapezoid_tail_fit_slower_than_physical_is_an_error():
    with pytest.raises(ValueError, match="slower than the physical decay.*use tail 'phys'"):
        integrate_trapezoid([24.0, 96.0, 168.0], [10.0, 10.0, 9.9], LP, "zero", "fit", tail_points=2)


def test_trapezoid_extrap_needs_positive_values():
    with pytest.raises(ValueError, match="first two values > 0"):
        integrate_trapezoid([2.0, 4.0], [0.0, 6.0], LP, "extrap", "phys")


# ----------------------------------------------------------------------------------------------
# single time point
# ----------------------------------------------------------------------------------------------


def test_stp_models():
    t, a = 24.0, 50.0
    r = integrate_stp(t, a, LP, "physical")
    assert r["tia_mbq_h"] == pytest.approx(a * math.exp(LP * t) / LP)  # exact for a pure physical decay
    assert r["fits"]["physical"]["teff_h"] == pytest.approx(T_PHYS)
    assert integrate_stp(t, a, LP, "hanscheid2018", 48.0, "default")["tia_mbq_h"] == pytest.approx(a * 48 / math.log(2))
    lam = math.log(2) / 60.0
    r = integrate_stp(t, a, LP, "madsen2018", 60.0, "given")
    assert r["tia_mbq_h"] == pytest.approx(a * math.exp(lam * t) / lam)
    assert r["fits"]["madsen2018"]["teff_source"] == "given"


# ----------------------------------------------------------------------------------------------
# options: no default, an option the model does not use is an error
# ----------------------------------------------------------------------------------------------


@pytest.mark.parametrize("model, n, kwargs, err", [
    ("trapezoid", 3, {"tail": "phys"}, "needs head"),
    ("trapezoid", 3, {"head": "ramp"}, "needs tail"),
    ("trapezoid", 3, {"head": "ramp", "tail": "fit"}, "needs tail_points"),
    ("trapezoid", 3, {"head": "ramp", "tail": "fit", "tail_points": 4}, "between 2 and"),
    ("trapezoid", 3, {"head": "ramp", "tail": "phys", "tail_points": 2}, "only used with tail 'fit'"),
    ("trapezoid", 3, {"head": "ramp", "tail": "phys", "error_model": "proportional"}, "Not used.*error_model"),
    ("trapezoid", 1, {"head": "ramp", "tail": "phys"}, "at least 2 time points"),
    ("f2", 4, {"error_model": "proportional", "head": "ramp"}, "Not used.*head"),
    ("f2", 4, {"error_model": "proportional", "teff_given": True}, "Not used.*teff"),
    ("f2", 4, {}, "Missing error model"),
    ("f2", 1, {"error_model": "fsd:0.1"}, "single time point model"),
    ("hanscheid2018", 2, {}, "give one image"),
    ("physical", 1, {"teff_given": True}, "Not used.*teff"),
    ("physical", 1, {"error_model": "fsd:0.1"}, "Not used.*error_model"),
    ("biexp", 3, {}, "Unknown model"),
])
def test_options(model, n, kwargs, err):
    with pytest.raises(ValueError, match=err):
        check_integration_options(model, n, **kwargs)


# ----------------------------------------------------------------------------------------------
# ROI level
# ----------------------------------------------------------------------------------------------


def test_roi_trapezoid(tmp_path):
    images, labels = _study(tmp_path)
    doc = integrate_rois(images, labels, model="trapezoid", head="ramp", tail="phys")
    liver = doc["rois"]["Liver"]
    a = np.array(liver["activities_mbq"])
    t = np.array(TIMES)
    expected = a[0] * t[0] / 2 + np.sum((a[1:] + a[:-1]) * np.diff(t) / 2) + a[-1] / LP
    assert liver["tia_mbq_h"] == pytest.approx(expected)
    assert doc["head"] == "ramp" and doc["tail"] == "phys" and doc["error_model"] is None


def test_roi_madsen_needs_teff(tmp_path):
    images, labels = _study(tmp_path)
    with pytest.raises(ValueError, match="teff\\) must be provided for ROI"):
        integrate_rois(images[:1], labels, model="madsen2018")
    doc = integrate_rois(images[:1], labels, model="madsen2018", teff={"Liver": 30.0, "Spleen": 10.0})
    assert doc["rois"]["Liver"]["fits"]["madsen2018"]["teff_h"] == 30.0 and doc["teff_h"]["Spleen"] == 10.0


@pytest.mark.parametrize("method, teff", [("hanscheid2018", None), ("hanscheid2018", {"Liver": 30.0}),
                                          ("madsen2018", {"Liver": 30.0, "Spleen": 12.0})])
def test_stp_integrate_then_svalue_equals_dose_method(tmp_path, method, teff):
    """integrate --model M + dose -m svalue = dose -m M (same masses: labels, no CT)."""
    images, labels = _study(tmp_path, sex="male")
    doc = integrate_rois(images[1:2], labels, model=method, teff=teff)
    via_table = compute_roi_dose_tia_table(doc, labels)
    params = {r: {"teff": v} for r, v in (teff or {}).items()}
    direct = compute_roi_dose(image=sitk.ReadImage(str(images[1])), image_path=images[1], method=method,
                              isotope="Lu177", labels=sitk.ReadImage(str(labels)), labels_path=labels,
                              method_params={**params, "phantom": "male"})
    for roi in ("Liver", "Spleen"):
        # rel 1e-6: the dose path converts bqml -> bq in the float32 image, integrate in float64
        assert via_table[roi]["mean_dose_gy"] == pytest.approx(direct[roi]["mean_dose_gy"], rel=1e-6)


def test_cli_trapezoid_and_stp(tmp_path):
    images, labels = _study(tmp_path)
    out = tmp_path / "tia.json"
    r = CliRunner().invoke(main, ["integrate", *map(str, images), "-l", str(labels), "--model", "trapezoid",
                                  "--head", "extrap", "--tail", "fit", "--tail-points", "2", "-o", str(out)])
    assert r.exit_code == 0, r.output
    assert read_integrate(out)["rois"]["Liver"]["model"] == "trapezoid"
    r = CliRunner().invoke(main, ["integrate", str(images[1]), "-l", str(labels), "--model", "madsen2018",
                                  "--teff", "Liver:30,Spleen:12", "-o", str(out)])
    assert r.exit_code == 0, r.output
    assert read_integrate(out)["teff_h"] == {"Liver": 30.0, "Spleen": 12.0}
    r = CliRunner().invoke(main, ["integrate", str(images[1]), "-l", str(labels), "--model", "madsen2018",
                                  "--teff", "Liver=30"])
    assert r.exit_code != 0 and "ROI:hours" in r.output
