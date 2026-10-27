#!/usr/bin/env python3
"""
Read the OpenGATE dose-rate outputs (in data/) and write p142_<radionuclide>_opengate_benchmark.json
for every p142_dose_rate_<radionuclide>.json configuration found next to this script.
"""
import json
from pathlib import Path
import numpy as np
import SimpleITK as sitk
from radot.db import get_physical_half_life_h
from radot.dose import (
    roi_dose_uncertainty_bounds,
    roi_dose_uncertainty_estimate,
)


def process_single_benchmark(config_path: Path, p142_dir: Path, data_dir: Path) -> bool:
    with open(config_path) as f:
        config = json.load(f)

    radionuclide = config.get("radionuclide", "Lu177")
    rad_lower = radionuclide.lower()
    total_activity_bq = float(config.get("activity_bq", 1e6))
    mode = config.get("mode", "vrt")

    # Determine default Teff: 45h for Lu177, otherwise physical half life (e.g. 64.1h for Y90)
    if "teff_h" in config:
        teff_h = float(config["teff_h"])
    elif rad_lower == "lu177":
        teff_h = 45.0
    else:
        teff_h = get_physical_half_life_h(radionuclide)

    time_from_injection_h = float(config.get("time_from_injection_h", 0.0))

    # Determine candidate output folder for OpenGate
    candidate_folders = []
    if "output_folder" in config and config["output_folder"]:
        candidate_folders.append(data_dir / config["output_folder"])
    else:
        candidate_folders.append(data_dir / f"output_{rad_lower}_{mode}")
        candidate_folders.append(data_dir / f"output_{rad_lower}")
        # Default output_vrt only if Lu177 or only 1 config
        if rad_lower == "lu177":
            candidate_folders.append(data_dir / f"output_{mode}")

    gate_dose_path = None
    gate_unc_path = None
    for folder in candidate_folders:
        dose_file = folder / "output_dose.mhd"
        if dose_file.exists():
            gate_dose_path = dose_file
            unc_file = folder / "output_dose_uncertainty.mhd"
            if unc_file.exists():
                gate_unc_path = unc_file
            break

    if not gate_dose_path:
        print(
            f"[-] Skipping {radionuclide} ({config_path.name}): simulation output not found in any of:"
        )
        for c in candidate_folders:
            print(f"      {c}")
        return False

    liver_path = data_dir / "liver_4mm_cropped.nii.gz"
    ct_path = data_dir / "ct_4mm_cropped.nii.gz"
    if not liver_path.exists() or not ct_path.exists():
        raise FileNotFoundError(f"Liver or CT images missing in {data_dir}.")

    # 1. Load liver mask and CT
    liver_img = sitk.ReadImage(str(liver_path))
    mask = sitk.GetArrayViewFromImage(liver_img) > 0
    num_voxels = int(mask.sum())

    spacing = liver_img.GetSpacing()
    voxel_volume_ml = (spacing[0] * spacing[1] * spacing[2]) / 1000.0
    roi_volume_ml = float(num_voxels * voxel_volume_ml)

    ct_img = sitk.ReadImage(str(ct_path))
    ct_arr = sitk.GetArrayViewFromImage(ct_img)
    density = np.maximum(0.0, 1.0 + (ct_arr[mask] / 1000.0))
    roi_mass_g = float(np.sum(density) * voxel_volume_ml)

    # 2. Read OpenGate simulated dose rate map (Gy/s) and uncertainty
    gate_img = sitk.ReadImage(str(gate_dose_path))
    gate_arr = sitk.GetArrayViewFromImage(gate_img)
    gate_mean_dose_rate_gy_s = float(gate_arr[mask].mean())

    opengate_dict = {
        "mean_dose_rate_gy_s": gate_mean_dose_rate_gy_s,
        "simulation_mode": mode,
    }
    if gate_unc_path and gate_unc_path.exists():
        unc_img = sitk.ReadImage(str(gate_unc_path))
        unc_min, unc_max = roi_dose_uncertainty_bounds(gate_img, unc_img, mask)
        unc_est = roi_dose_uncertainty_estimate(
            gate_img, unc_img, mask, isotope=radionuclide
        )
        opengate_dict["uncertainty_min_fraction"] = unc_min
        opengate_dict["uncertainty_max_fraction"] = unc_max
        opengate_dict["uncertainty_estimate_fraction"] = unc_est
        opengate_dict["mean_voxel_uncertainty_fraction"] = unc_max

    benchmark_data = {
        "description": f"Patient 142 cropped liver {radionuclide} OpenGate VRT simulation benchmark",
        "patient": "p142",
        "roi_name": "Liver",
        "radionuclide": radionuclide,
        "phantom": "male",
        "total_activity_bq": total_activity_bq,
        "teff_h": teff_h,
        "time_from_injection_h": time_from_injection_h,
        "num_voxels": num_voxels,
        "roi_volume_ml": roi_volume_ml,
        "roi_mass_g": roi_mass_g,
        "opengate": opengate_dict,
    }

    # 3. Read cross-dose to Right Kidney (Kidney_R, label ID 6) if labels are present
    labels_crop_path = data_dir / "labels_all_4mm_cropped.nii.gz"
    if not labels_crop_path.exists():
        labels_full_path = data_dir / "labels_all_4mm.nii.gz"
        if labels_full_path.exists():
            full_lbl = sitk.ReadImage(str(labels_full_path))
            lbl_crop = full_lbl[:, :, 282:366]
            sitk.WriteImage(lbl_crop, str(labels_crop_path))

    if labels_crop_path.exists():
        lbl_crop = sitk.ReadImage(str(labels_crop_path))
        lbl_crop_arr = sitk.GetArrayViewFromImage(lbl_crop)
        mask_rk = lbl_crop_arr == 6
        num_vox_rk = int(mask_rk.sum())

        rk_vol_ml = float(num_vox_rk * voxel_volume_ml)
        density_rk = np.maximum(0.0, 1.0 + (ct_arr[mask_rk] / 1000.0))
        rk_mass_g = float(np.sum(density_rk) * voxel_volume_ml)

        gate_rk_mean_rate = float(gate_arr[mask_rk].mean())
        rk_opengate = {
            "mean_dose_rate_gy_s": gate_rk_mean_rate,
        }
        if gate_unc_path and gate_unc_path.exists():
            unc_img = sitk.ReadImage(str(gate_unc_path))
            rk_unc_min, rk_unc_max = roi_dose_uncertainty_bounds(
                gate_img, unc_img, mask_rk
            )
            rk_unc_est = roi_dose_uncertainty_estimate(
                gate_img, unc_img, mask_rk, isotope=radionuclide
            )
            rk_opengate["uncertainty_min_fraction"] = rk_unc_min
            rk_opengate["uncertainty_max_fraction"] = rk_unc_max
            rk_opengate["uncertainty_estimate_fraction"] = rk_unc_est

        benchmark_data["cross_dose"] = {
            "target_roi_name": "Kidney_R",
            "source_roi_name": "Liver",
            "target_num_voxels": num_vox_rk,
            "target_volume_ml": rk_vol_ml,
            "target_mass_g": rk_mass_g,
            "opengate": rk_opengate,
        }

    output_json_path = p142_dir / f"p142_{rad_lower}_opengate_benchmark.json"
    with open(output_json_path, "w") as f:
        json.dump(benchmark_data, f, indent=2)

    print(f"[+] Benchmark created for {radionuclide}:")
    print(f"    Source config:       {config_path.name}")
    print(f"    Simulation output:   {gate_dose_path.parent.name}")
    print(f"    Total activity:      {total_activity_bq:.2e} Bq")
    print(f"    Teff:                {teff_h:.1f} h")
    print(f"    Liver dose rate:     {gate_mean_dose_rate_gy_s:.6e} Gy/s")
    print(f"    Saved to:            {output_json_path}")
    return True


def main():
    p142_dir = Path(__file__).resolve().parent
    data_dir = p142_dir / "data"
    config_files = sorted(p142_dir.glob("p142_dose_rate_*.json"))
    if not config_files:
        print("No p142_dose_rate_*.json configurations found.")
        return

    print(f"Found {len(config_files)} configuration(s) in {p142_dir}:")
    for cfg in config_files:
        print(f"  - {cfg.name}")
    print()

    count = 0
    for cfg in config_files:
        if process_single_benchmark(cfg, p142_dir, data_dir):
            count += 1
            print()

    print(f"Done. Processed {count}/{len(config_files)} benchmark(s).")


if __name__ == "__main__":
    main()
