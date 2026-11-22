"""No CT: the masses use water density 1 g/mL, said once per run with a WARNING (stats.density_map_g_ml),
on every dose path; never with a CT."""
import logging

import numpy as np
import pytest
import SimpleITK as sitk

from radot.dose import compute_roi_dose, compute_roi_dose_mtp, compute_roi_dose_tia_table
from radot.integrate import integrate_rois
from radot.stats import WATER_DENSITY_WARNING, density_map_g_ml

from test_dose_input_errors import _mtp_study
from test_dose_rate_map import svalue_map
from test_dose_tia import make_image, make_labels
from test_integrate import _study

COMMON = dict(isotope="Lu177", labels_map={1: "Left kidney"})


def _ct(hu=0.0):
    ct = sitk.Image([5, 5, 5], sitk.sitkFloat32) + hu
    ct.SetSpacing((10.0, 10.0, 10.0))
    return ct


def _count(caplog):
    return sum(WATER_DENSITY_WARNING in r.message and r.levelno == logging.WARNING for r in caplog.records)


def test_density_map():
    assert density_map_g_ml(np.array([0.0, 1000.0])).tolist() == [1.0, 2.0]
    assert density_map_g_ml(None, (2, 3), log=False).tolist() == [[1.0] * 3] * 2


RUNS = {
    "svalue": lambda ct: compute_roi_dose(image=make_image(1e12, "bqs"), labels=make_labels(), method="svalue",
                                          ct_image=ct, **COMMON),
    "hanscheid2018": lambda ct: compute_roi_dose(image=make_image(1e9, "bq", 24.0), labels=make_labels(),
                                                 method="hanscheid2018", ct_image=ct, **COMMON),
    "ldm": lambda ct: compute_roi_dose(image=make_image(1e9, "bq", 24.0), labels=make_labels(), method="ldm",
                                       ct_image=ct, **COMMON),
    "map": lambda ct: compute_roi_dose(image=make_image(1e9, "bq", 24.0), labels=make_labels(),
                                       method="hanscheid2018", dose_rate_map=svalue_map(),
                                       dose_rate_map_activity_mbq=1.0, ct_image=ct, **COMMON),
}


@pytest.mark.parametrize("path", list(RUNS))
def test_single_time_point_paths(path, caplog):
    with caplog.at_level(logging.WARNING, logger="radot"):
        no_ct = RUNS[path](None)
    assert _count(caplog) == 1
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="radot"):
        water_ct = RUNS[path](_ct(0.0))
    assert _count(caplog) == 0
    if isinstance(no_ct, tuple):  # ldm: (ROI table, dose image)
        no_ct, water_ct = no_ct[0], water_ct[0]
    assert no_ct["Left kidney"]["roi_mass_g"] == pytest.approx(water_ct["Left kidney"]["roi_mass_g"])
    assert no_ct["Left kidney"]["mean_dose_gy"] == pytest.approx(water_ct["Left kidney"]["mean_dose_gy"])


def test_mtp(tmp_path, caplog):
    spects, labels = _mtp_study(tmp_path)
    with caplog.at_level(logging.WARNING, logger="radot"):
        compute_roi_dose_mtp(spects, labels, method="jackson2013", isotope="Lu177", labels_path=labels)
    assert _count(caplog) == 1  # once, not once per time point


def test_tia_table(tmp_path, caplog):
    images, labels = _study(tmp_path, sex="male")
    doc = integrate_rois(images, labels, model="f2", error_model="proportional")
    with caplog.at_level(logging.WARNING, logger="radot"):
        compute_roi_dose_tia_table(doc, labels)
    assert _count(caplog) == 1
