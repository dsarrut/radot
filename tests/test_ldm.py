import pytest
import numpy as np
import SimpleITK as sitk
import json
from click.testing import CliRunner
from radot.cli import main
from radot.db import get_non_penetrating_energy_constant_mev
from radot.dose import compute_roi_dose


def test_get_non_penetrating_energy_constant_mev():
    # Lu-177 non-penetrating energy
    lu_const = get_non_penetrating_energy_constant_mev("Lu177")
    assert lu_const > 0.1  # Typically around 0.14 MeV

    # Y-90 non-penetrating energy
    y_const = get_non_penetrating_energy_constant_mev("Y90")
    assert y_const > 0.8  # Typically around 0.93 MeV


def test_decay_chains_and_cascade_energy():
    from radot.db import get_decay_chain, get_cascade_energy_constant_mev

    # 1. Non-cascade isotope returns single Delta_local
    lu_tot, lu_contrib = get_cascade_energy_constant_mev("Lu177", cascade=True)
    assert lu_tot == pytest.approx(get_non_penetrating_energy_constant_mev("Lu177"))
    assert "Lu177" in lu_contrib

    # 2. Ac-225 decay chain lookup
    chain = get_decay_chain("Ac225")
    assert chain is not None
    assert chain["parent"] == "Ac225"
    daughters = [d["nuclide"] for d in chain["daughters"]]
    assert "Fr221" in daughters
    assert "Bi213" in daughters
    assert "Pb209" in daughters

    # 3. Ac-225 cascade Delta_local in secular equilibrium (~28.7 MeV vs 5.92 MeV parent alone)
    ac_parent, _ = get_cascade_energy_constant_mev("Ac225", cascade=False)
    ac_cascade, ac_contrib = get_cascade_energy_constant_mev("Ac225", cascade=True)
    assert ac_parent == pytest.approx(5.9168, rel=1e-3)
    assert ac_cascade == pytest.approx(28.704, rel=1e-2)
    assert ac_cascade / ac_parent == pytest.approx(4.85, rel=1e-2)
    assert len(ac_contrib) >= 6

    # 4. Ra-223 and Pb-212 cascade constants
    ra_cascade, _ = get_cascade_energy_constant_mev("Ra223", cascade=True)
    assert ra_cascade > 25.0  # 4 alphas cascade (~28 MeV)

    pb_cascade, _ = get_cascade_energy_constant_mev("Pb212", cascade=True)
    assert pb_cascade > 7.0  # 1 alpha cascade (~8 MeV)


def test_ldm_core_calculation_without_ct():
    # Setup test SPECT image (2x2x2) with voxel values of 1e6 Bq
    image = sitk.Image([2, 2, 2], sitk.sitkFloat32)
    # Voxel spacing: 2mm isotropic
    image.SetSpacing([2.0, 2.0, 2.0])
    image_arr = np.ones((2, 2, 2), dtype=np.float32) * 1e6

    # Write array to SITK image
    sitk_image = sitk.GetImageFromArray(image_arr)
    sitk_image.SetSpacing([2.0, 2.0, 2.0])
    sitk_image.SetMetaData("time_from_injection", "24.0")
    sitk_image.SetMetaData("unit", "Bq")

    # Run LDM without CT (assumes water density = 1.0 g/cm3)
    dose_img = compute_roi_dose(
        image=sitk_image,
        labels=None,
        method="ldm",
        isotope="Lu177",
        ct_image=None,
        verbose=True,
    )

    assert isinstance(dose_img, sitk.Image)
    assert dose_img.GetSize() == (2, 2, 2)
    assert dose_img.GetSpacing() == (2.0, 2.0, 2.0)

    # Let's inspect the calculated voxel values
    dose_arr = sitk.GetArrayFromImage(dose_img)
    # Dose value should be > 0.0 Gy since voxel has activity
    assert np.all(dose_arr > 0.0)


def test_ldm_core_calculation_with_ct():
    # Setup activity image
    image = sitk.Image([2, 2, 2], sitk.sitkFloat32)
    image.SetSpacing([2.0, 2.0, 2.0])
    image_arr = np.ones((2, 2, 2), dtype=np.float32) * 1e6
    sitk_image = sitk.GetImageFromArray(image_arr)
    sitk_image.SetSpacing([2.0, 2.0, 2.0])
    sitk_image.SetMetaData("time_from_injection", "24.0")
    sitk_image.SetMetaData("unit", "Bq")

    # Setup CT image with density corresponding to bone (HU = 1000)
    ct_arr = np.ones((2, 2, 2), dtype=np.int16) * 1000
    ct_image = sitk.GetImageFromArray(ct_arr)
    ct_image.SetSpacing([2.0, 2.0, 2.0])

    # Run LDM with CT
    dose_img = compute_roi_dose(
        image=sitk_image,
        labels=None,
        method="ldm",
        isotope="Lu177",
        ct_image=ct_image,
        verbose=True,
    )

    dose_arr = sitk.GetArrayFromImage(dose_img)
    assert np.all(dose_arr > 0.0)


def test_ldm_cli_execution(tmp_path):
    runner = CliRunner()

    image_path = tmp_path / "spect.nii.gz"
    ct_path = tmp_path / "ct.nii.gz"
    output_path = tmp_path / "dose.nii.gz"

    # Save SPECT image
    img = sitk.Image([3, 3, 3], sitk.sitkFloat32)
    img.SetSpacing([2.0, 2.0, 2.0])
    sitk.WriteImage(img, str(image_path))

    with open(tmp_path / "spect.json", "w") as f:
        json.dump({"unit": "Bq", "time_from_injection": "24.0"}, f)

    # Save CT image
    ct = sitk.Image([3, 3, 3], sitk.sitkFloat32)
    ct.SetSpacing([2.0, 2.0, 2.0])
    sitk.WriteImage(ct, str(ct_path))

    # Run LDM CLI command
    result = runner.invoke(
        main,
        [
            "dose",
            str(image_path),
            "-m",
            "ldm",
            "--isotope",
            "Lu177",
            "--ct",
            str(ct_path),
            "-o",
            str(output_path),
        ],
    )

    assert result.exit_code == 0
    assert output_path.exists()

    json_path = output_path.with_suffix("").with_suffix(".json")
    assert json_path.exists()
    with open(json_path, "r") as f:
        meta_data = json.load(f)
    assert meta_data["unit"] == "gy"
    assert meta_data["modality"] == "dose"


def test_ldm_core_calculation_with_labels_no_output():
    # Setup test SPECT image (2x2x2) with voxel values of 1e6 Bq
    image = sitk.Image([2, 2, 2], sitk.sitkFloat32)
    image.SetSpacing([2.0, 2.0, 2.0])
    image_arr = np.ones((2, 2, 2), dtype=np.float32) * 1e6
    sitk_image = sitk.GetImageFromArray(image_arr)
    sitk_image.SetSpacing([2.0, 2.0, 2.0])
    sitk_image.SetMetaData("time_from_injection", "24.0")
    sitk_image.SetMetaData("unit", "Bq")

    # Setup labels image
    labels = sitk.Image([2, 2, 2], sitk.sitkUInt16)
    labels.SetSpacing([2.0, 2.0, 2.0])
    labels[1, 1, 1] = 1

    # Run LDM with labels
    results = compute_roi_dose(
        image=sitk_image,
        labels=labels,
        method="ldm",
        isotope="Lu177",
        labels_map={1: "Liver"},
        verbose=True,
    )

    assert isinstance(results, tuple)
    results_dict, dose_img = results
    assert "Liver" in results_dict
    assert results_dict["Liver"]["mean_dose_gy"] > 0.0
    assert results_dict["Liver"]["roi_volume_ml"] > 0.0
    assert isinstance(dose_img, sitk.Image)


def test_ldm_cli_with_labels_no_output(tmp_path):
    runner = CliRunner()

    image_path = tmp_path / "spect.nii.gz"
    labels_path = tmp_path / "labels.nii.gz"

    # Save SPECT image
    img = sitk.Image([3, 3, 3], sitk.sitkFloat32)
    img.SetSpacing([2.0, 2.0, 2.0])
    sitk.WriteImage(img, str(image_path))

    with open(tmp_path / "spect.json", "w") as f:
        json.dump({"unit": "Bq", "time_from_injection": "24.0"}, f)

    # Save labels image
    labels = sitk.Image([3, 3, 3], sitk.sitkUInt16)
    labels.SetSpacing([2.0, 2.0, 2.0])
    labels[1, 1, 1] = 1
    sitk.WriteImage(labels, str(labels_path))

    with open(tmp_path / "labels.json", "w") as f:
        json.dump({"1": "Liver"}, f)

    # Run LDM CLI command without output path option -o
    result = runner.invoke(
        main,
        [
            "dose",
            str(image_path),
            "-l",
            str(labels_path),
            "-m",
            "ldm",
            "--isotope",
            "Lu177",
        ],
    )

    assert result.exit_code == 0
    assert "Liver" in result.output


@pytest.mark.parametrize(
    "output, image_name, table_name",
    [("dose.nii.gz", "dose.nii.gz", "dose_rois.json"), ("results.json", "results_dose.nii.gz", "results.json")],
)
def test_ldm_with_labels_writes_image_and_table(tmp_path, output, image_name, table_name):
    # voxel method with ROIs: the dose image and the ROI table are both written, whatever -o is
    image_path = tmp_path / "spect.nii.gz"
    labels_path = tmp_path / "labels.nii.gz"
    img = sitk.Image([3, 3, 3], sitk.sitkFloat32) + 1000.0
    img.SetSpacing([2.0, 2.0, 2.0])
    sitk.WriteImage(img, str(image_path))
    (tmp_path / "spect.json").write_text(json.dumps({"unit": "bq", "time_from_injection": "24.0"}))
    labels = sitk.Image([3, 3, 3], sitk.sitkUInt16)
    labels.SetSpacing([2.0, 2.0, 2.0])
    labels[1, 1, 1] = 1
    sitk.WriteImage(labels, str(labels_path))
    (tmp_path / "labels.json").write_text(json.dumps({"1": "Liver"}))

    out = tmp_path / output
    result = CliRunner().invoke(
        main, ["dose", str(image_path), "-l", str(labels_path), "-m", "ldm", "--isotope", "Lu177", "-o", str(out)]
    )
    assert result.exit_code == 0, result.output
    assert "Liver" in result.output
    assert sitk.ReadImage(str(tmp_path / image_name)).GetSize() == (3, 3, 3)
    sidecar = json.loads((tmp_path / image_name.replace(".nii.gz", ".json")).read_text())
    assert sidecar["unit"] == "gy" and sidecar["modality"] == "dose"
    table = json.loads((tmp_path / table_name).read_text())
    assert table["Liver"]["mean_dose_gy"] > 0


def test_ldm_hu_threshold():
    # Activity image (2x2x2) with activity = 1e6 Bq
    image = sitk.Image([2, 2, 2], sitk.sitkFloat32)
    image.SetSpacing([2.0, 2.0, 2.0])
    image_arr = np.ones((2, 2, 2), dtype=np.float32) * 1e6
    sitk_image = sitk.GetImageFromArray(image_arr)
    sitk_image.SetSpacing([2.0, 2.0, 2.0])
    sitk_image.SetMetaData("time_from_injection", "24.0")
    sitk_image.SetMetaData("unit", "Bq")

    # CT image where one voxel has HU = -900, another has HU = 0
    ct_arr = np.zeros((2, 2, 2), dtype=np.int16)
    ct_arr[0, 0, 0] = -900
    ct_image = sitk.GetImageFromArray(ct_arr)
    ct_image.SetSpacing([2.0, 2.0, 2.0])

    # Run LDM with hu_threshold = -800
    dose_img = compute_roi_dose(
        image=sitk_image,
        labels=None,
        method="ldm",
        isotope="Lu177",
        ct_image=ct_image,
        hu_threshold=-800.0,
        verbose=True,
    )

    dose_arr = sitk.GetArrayFromImage(dose_img)
    # The voxel at (0, 0, 0) should be zeroed out
    assert dose_arr[0, 0, 0] == 0.0
    # The other voxels (HU = 0 >= -800) should have dose > 0.0
    assert dose_arr[1, 1, 1] > 0.0


def test_ldm_dose_rate(tmp_path):
    runner = CliRunner()

    image_path = tmp_path / "spect.nii.gz"
    labels_path = tmp_path / "labels.nii.gz"
    output_path = tmp_path / "dose_rate.nii.gz"

    # Save SPECT image
    img = sitk.Image([3, 3, 3], sitk.sitkFloat32)
    img.SetSpacing([2.0, 2.0, 2.0])
    sitk.WriteImage(img, str(image_path))

    with open(tmp_path / "spect.json", "w") as f:
        json.dump({"unit": "Bq", "time_from_injection": "24.0"}, f)

    # Save labels image
    labels = sitk.Image([3, 3, 3], sitk.sitkUInt16)
    labels.SetSpacing([2.0, 2.0, 2.0])
    labels[1, 1, 1] = 1
    sitk.WriteImage(labels, str(labels_path))

    with open(tmp_path / "labels.json", "w") as f:
        json.dump({"1": "Liver"}, f)

    # Run LDM CLI command with --dose-rate
    result = runner.invoke(
        main,
        [
            "dose",
            str(image_path),
            "-l",
            str(labels_path),
            "-m",
            "ldm",
            "--isotope",
            "Lu177",
            "-o",
            str(output_path),
            "--dose-rate",
        ],
    )

    assert result.exit_code == 0
    assert "Mean Dose Rate (Gy/h)" in result.output
    assert output_path.exists()

    json_path = output_path.with_suffix("").with_suffix(".json")
    assert json_path.exists()
    with open(json_path, "r") as f:
        meta_data = json.load(f)
    assert meta_data["unit"] == "gy/h"
    assert meta_data["modality"] == "dose"


def _fake_chain(monkeypatch, daughters):
    import radot.db as db

    # each entry is the alpha daughter of the previous one (Z - 2, A - 4), short-lived after the parent
    chain = {"parent": "Ac225", "half_life_s": 857088.0, "half_life_str": "9.92 d", "daughters": [
        {"nuclide": n, "yield": y, "z": 89 - 2 * i, "a": 225 - 4 * i, "half_life_str": "1 s",
         "half_life_s": 857088.0 if i == 0 else 1.0, "is_stable": False}
        for i, (n, y) in enumerate(daughters)]}
    monkeypatch.setattr(db, "_decay_chains_cache", {"chains": {"Ac225": chain}})


def test_cascade_daughter_without_icrp107_data_raises(monkeypatch):
    """A daughter without ICRP 107 data is an error, not silently dropped from the cascade sum."""
    from radot.db import get_cascade_energy_constant_mev

    _fake_chain(monkeypatch, [("Ac225", 1.0), ("Xx999", 0.02)])
    with pytest.raises(ValueError, match="Xx999.*0.02.*--no-cascade"):
        get_cascade_energy_constant_mev("Ac225", cascade=True)
    # without the cascade, the chain is not used
    tot, _ = get_cascade_energy_constant_mev("Ac225", cascade=False)
    assert tot == pytest.approx(get_non_penetrating_energy_constant_mev("Ac225"))


def test_cascade_daughter_missing_in_radionuclide_info(monkeypatch):
    """radot radionuclide info still shows the chain, with the reason the total is not available."""
    _fake_chain(monkeypatch, [("Ac225", 1.0), ("Fr221", 1.0), ("Xx999", 0.02)])
    result = CliRunner().invoke(main, ["radionuclide", "info", "Ac225"])
    assert result.exit_code == 0, result.output
    assert "Xx999" in result.output and "N/A" in result.output
    assert "not available" in result.output


def test_decay_chains_unreadable_raises(monkeypatch, tmp_path):
    """An unreadable decay_chains.json is an error: the cascade must not become single-nuclide."""
    import importlib.resources
    import radot.db as db

    (tmp_path / "data").mkdir()
    (tmp_path / "data" / "decay_chains.json").write_text("{ broken")
    monkeypatch.setattr(db, "_decay_chains_cache", None)
    monkeypatch.setattr(importlib.resources, "files", lambda pkg: tmp_path)
    with pytest.raises(RuntimeError, match="decay_chains.json"):
        db.get_cascade_energy_constant_mev("Ac225", cascade=True)
    assert db._decay_chains_cache is None  # the failure is not cached


# Cascade: only the daughters in secular equilibrium with the parent (T1/2 << T1/2 parent)

@pytest.mark.parametrize("parent, excluded", [
    ("Ac225", set()),
    ("Pb212", set()),
    ("Ra223", set()),
    ("At211", {"Bi207"}),  # Po211 (0.5 s) kept, Bi207 (31.6 y) excluded
    ("Th227", {"Ra223", "Rn219", "Po215", "Pb211", "Bi211", "Tl207", "Po211"}),  # Ra223 (11.4 d) and below
    ("Tb149", {"Gd149", "Eu145", "Eu149", "Sm145", "Pm145"}),
])
def test_cascade_members(parent, excluded):
    from radot.db import cascade_members, get_decay_chain

    kept, out = cascade_members(get_decay_chain(parent))
    assert {d["nuclide"] for d in out} == excluded
    assert kept[0]["nuclide"] == parent


def test_cascade_energy_excludes_long_lived_daughters(caplog):
    from radot.db import get_cascade_energy_constant_mev

    # Ac225: unchanged (all daughters short-lived), the OpenGATE validation case
    ac, _ = get_cascade_energy_constant_mev("Ac225", cascade=True)
    assert ac == pytest.approx(28.704, rel=1e-3)
    # Tb149: all the daughters live longer than T1/2 / 10: the parent alone (was 0.857 MeV, +12 %)
    tb, contrib = get_cascade_energy_constant_mev("Tb149", cascade=True)
    assert tb == pytest.approx(get_non_penetrating_energy_constant_mev("Tb149"))
    assert set(contrib) == {"Tb149"}
    assert "Cascade of Tb149: not in secular equilibrium" in caplog.text and "Gd149" in caplog.text
    # At211: the Po211 alphas are kept, Bi207 is not
    _, contrib = get_cascade_energy_constant_mev("At211", cascade=True)
    assert set(contrib) == {"At211", "Po211"}


def test_radionuclide_info_marks_excluded():
    result = CliRunner().invoke(main, ["radionuclide", "info", "Tb149"])
    assert result.exit_code == 0, result.output
    assert "excluded" in result.output
    assert "Excluded (not in secular equilibrium" in result.output and "Gd149" in result.output
    assert "Cascade Factor: 1.00x" in result.output
