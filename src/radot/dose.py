import json
import logging
from typing import List, Optional, Dict, Any
import numpy as np
import SimpleITK as sitk
from pathlib import Path
from datetime import datetime

logger = logging.getLogger("radot")
from radot.db import get_physical_half_life_h
from radot.meta import get_time_from_injection_h


def dose_hanscheid2017(
    spect_bq_a: np.ndarray,
    roi_a: np.ndarray,
    time_from_injection_h: float,
    volume_voxel_m_l: float,
    effective_time_h: float,
) -> float:
    """
    Input img and ROI must be numpy arrays
    Input SPECT img must be in Bq
    """
    # compute mean activity concentration in the ROI
    v = spect_bq_a[roi_a == 1] / volume_voxel_m_l
    concentration = np.mean(v) / 1e6

    # compute dose
    dose = (
        0.125
        * concentration
        * np.power(2, time_from_injection_h / effective_time_h)
        * effective_time_h
    )

    return float(dose)


def dose_hanscheid2018(
    spect_bq_a: np.ndarray,
    roi_a: np.ndarray,
    time_from_injection_h: float,
    s_value: float,
    mass_scaling: float,
) -> float:
    """
    Input image and ROI must be numpy arrays
    - spect must be in Bq (not concentration)
    - acquisition time in hours
    - output is in Gray
    """
    # compute mean activity in the ROI, in MBq
    v = spect_bq_a[roi_a == 1]
    At = np.sum(v) / 1e6

    # S is in (mGy/MBq/s), so we get dose in mGy
    dose = (
        At
        * (2 * time_from_injection_h * 3600.0)
        / np.log(2)
        * s_value
        / mass_scaling
        / 1000.0
    )

    return float(dose)


def strategy_hanscheid2017(
    image_arr: np.ndarray,
    mask: np.ndarray,
    voxel_volume_ml: float,
    roi_name: str,
    time_from_injection_h: float,
    isotope: str,
    roi_mass_g: float,
    roi_volume_ml: float,
    method_params: Optional[Dict[str, Any]],
    verbose: bool,
    ct_arr: Optional[np.ndarray] = None,
) -> tuple[float, float]:
    """Strategy implementation for Hanscheid 2017 method."""
    # Teff parameter resolution
    teff = None
    if method_params and roi_name in method_params:
        teff = method_params[roi_name].get("teff")

    if teff is None:
        try:
            teff = get_physical_half_life_h(isotope)
        except ValueError as e:
            raise ValueError(
                f"No effective half-life (teff) provided for ROI '{roi_name}' and {e}"
            )
        logger.warning(
            f"No Teff specified for '{roi_name}'. Using physical half-life of {isotope} ({teff:.2f} h) as fallback."
        )

    dose = dose_hanscheid2017(
        spect_bq_a=image_arr,
        roi_a=mask.astype(np.uint8),
        time_from_injection_h=time_from_injection_h,
        volume_voxel_m_l=voxel_volume_ml,
        effective_time_h=teff,
    )

    if verbose or logger.isEnabledFor(logging.DEBUG):
        voxels_in_roi = np.sum(mask)
        total_activity_bq = np.sum(image_arr[mask])
        total_activity_mbq = total_activity_bq / 1e6
        mean_activity_concentration_bq_ml = np.mean(image_arr[mask] / voxel_volume_ml)
        mean_activity_concentration_mbq_ml = mean_activity_concentration_bq_ml / 1e6
        density_str = (
            f"computed using companion CT, mean density: {roi_mass_g / roi_volume_ml:.3f} g/mL"
            if ct_arr is not None
            else "no companion CT, assuming water density 1.0 g/mL"
        )
        formula_str = f"0.125 * {mean_activity_concentration_mbq_ml:.3e} * 2^({time_from_injection_h:.2f} / {teff:.2f}) * {teff:.2f}"

        log_msg = (
            f"\nDosimetry calculation log for ROI '{roi_name}' using Hanscheid 2017:\n"
            f"  - Radionuclide: {isotope}\n"
            f"  - Scaling factor: 0.125\n"
            f"  - ROI Volume: {roi_volume_ml:.3f} mL ({voxels_in_roi} voxels)\n"
            f"  - ROI Mass: {roi_mass_g:.3f} g\n"
            f"    ({density_str})\n"
            f"  - Total Activity in ROI: {total_activity_mbq:.3e} MBq\n"
            f"  - Mean Activity Concentration: {mean_activity_concentration_mbq_ml:.3e} MBq/mL ({mean_activity_concentration_bq_ml:.3e} Bq/mL)\n"
            f"  - Time from injection: {time_from_injection_h:.3f} h\n"
            f"  - Effective half-life (Teff): {teff:.3f} h\n"
            f"  - Formula: {formula_str}\n"
            f"  - Voxel-mean Absorbed Dose: {dose:.6f} Gy"
        )
        if verbose:
            logger.info(log_msg)
        else:
            logger.debug(log_msg)

    return float(dose), float(teff)


def strategy_hanscheid2018(
    image_arr: np.ndarray,
    mask: np.ndarray,
    voxel_volume_ml: float,
    roi_name: str,
    time_from_injection_h: float,
    isotope: str,
    roi_mass_g: float,
    roi_volume_ml: float,
    method_params: Optional[Dict[str, Any]],
    verbose: bool,
    ct_arr: Optional[np.ndarray] = None,
) -> tuple[float, float]:
    """Strategy implementation for Hanscheid 2018 method."""
    from radot.db import resolve_organ_name, get_s_value, get_reference_mass

    phantom = "male"
    if method_params:
        phantom = method_params.get("phantom", "male")

    # Standardize phantom/organ/rad
    opendose_organ = resolve_organ_name(roi_name, phantom)
    s_value = get_s_value(phantom, isotope, opendose_organ, opendose_organ)
    s_mass = get_reference_mass(phantom, opendose_organ)
    mass_scaling = roi_mass_g / s_mass

    dose = dose_hanscheid2018(
        spect_bq_a=image_arr,
        roi_a=mask.astype(np.uint8),
        time_from_injection_h=time_from_injection_h,
        s_value=s_value,
        mass_scaling=mass_scaling,
    )
    # Teff for Hanscheid 2018 is always assumed as 2 * t
    teff = 2.0 * time_from_injection_h

    if verbose or logger.isEnabledFor(logging.DEBUG):
        voxels_in_roi = np.sum(mask)
        total_activity_bq = np.sum(image_arr[mask])
        total_activity_mbq = total_activity_bq / 1e6
        density_str = (
            f"computed using companion CT, mean density: {roi_mass_g / roi_volume_ml:.3f} g/mL"
            if ct_arr is not None
            else "no companion CT, assuming water density 1.0 g/mL"
        )
        formula_str = f"At * (2 * t * 3600) / ln(2) * S / mass_scaling / 1000"

        log_msg = (
            f"\nDosimetry calculation log for ROI '{roi_name}' using Hanscheid 2018:\n"
            f"  - Radionuclide: {isotope}\n"
            f"  - Phantom model: {phantom}\n"
            f"  - OpenDose Organ Name: {opendose_organ}\n"
            f"  - Reference Organ Mass: {s_mass:.3f} g\n"
            f"  - ROI Volume: {roi_volume_ml:.3f} mL ({voxels_in_roi} voxels)\n"
            f"  - ROI Mass: {roi_mass_g:.3f} g\n"
            f"    ({density_str})\n"
            f"  - Mass scaling factor: {mass_scaling:.4f}\n"
            f"  - S-value (self): {s_value:.3e} mGy/MBq.s\n"
            f"  - Total Activity in ROI: {total_activity_mbq:.3e} MBq\n"
            f"  - Time from injection: {time_from_injection_h:.3f} h\n"
            f"  - Effective half-life (Teff = 2*t): {teff:.3f} h\n"
            f"  - Formula: {formula_str}\n"
            f"  - Voxel-mean Absorbed Dose: {dose:.6f} Gy"
        )
        if verbose:
            logger.info(log_msg)
        else:
            logger.debug(log_msg)

    return float(dose), float(teff)


DOSIMETRY_STRATEGIES = {
    "hanscheid2017": strategy_hanscheid2017,
    "hanscheid2018": strategy_hanscheid2018,
}


def compute_dose_generic(
    image_arr: np.ndarray,
    labels_arr: np.ndarray,
    isotope: str,
    voxel_volume_ml: float,
    time_from_injection_h: float,
    labels_map: Dict[int, str],
    ct_arr: Optional[np.ndarray] = None,
    rois: Optional[List[str]] = None,
    method_params: Optional[Dict[str, Any]] = None,
    verbose: bool = False,
    method: str = "hanscheid2017",
) -> dict:
    """
    Generic loop to calculate dosimetry using strategy functions.
    """
    method_key = method.lower()
    strategy = DOSIMETRY_STRATEGIES.get(method_key)
    if not strategy:
        raise ValueError(f"Unknown dosimetry method: {method}")

    unique_labels = sorted(np.unique(labels_arr))
    results = {}

    for label_id in unique_labels:
        if label_id == 0:
            continue

        roi_name = labels_map.get(label_id, f"unknown_{label_id}")

        # Filter by ROIs if specified
        if rois is not None and roi_name not in rois:
            continue

        mask = labels_arr == label_id
        if not np.any(mask):
            continue

        # Common calculations: volume and mass
        roi_volume_ml = np.sum(mask) * voxel_volume_ml

        if ct_arr is not None:
            density = np.maximum(0.0, 1.0 + (ct_arr[mask] / 1000.0))
            roi_mass_g = np.sum(density) * voxel_volume_ml
        else:
            roi_mass_g = roi_volume_ml

        # Delegate calculation and logging to the strategy
        try:
            dose, teff = strategy(
                image_arr=image_arr,
                mask=mask,
                voxel_volume_ml=voxel_volume_ml,
                roi_name=roi_name,
                time_from_injection_h=time_from_injection_h,
                isotope=isotope,
                roi_mass_g=roi_mass_g,
                roi_volume_ml=roi_volume_ml,
                method_params=method_params,
                verbose=verbose,
                ct_arr=ct_arr,
            )
            results[roi_name] = {
                "mean_dose_gy": float(dose),
                "roi_mass_g": float(roi_mass_g),
                "roi_volume_ml": float(roi_volume_ml),
                "roi_teff_h": float(teff),
            }
        except ValueError as e:
            if "not found in" in str(e) and (
                "S-value matrix" in str(e) or "phantom" in str(e)
            ):
                if rois is None:
                    logger.warning(f"Skipping dosimetry for ROI '{roi_name}': {e}")
                else:
                    raise ValueError(
                        f"Dosimetry Halt: {e}\n"
                        "Note: OpenDose lacks a single 'Skeleton' compartment and instead requires specific cortical/spongiosa targets.\n"
                        "To proceed, you must explicitly provide an external S-value override file via the CLI parameters (e.g., using --custom-s-values override.json) or adjust your input sidecar JSON mapping definitions to explicitly match target OpenDose bone/tissue definitions."
                    )
            else:
                raise e

    return results


def check_dosimetry_inputs(
    method: str,
    isotope: str,
    labels_path: Optional[Path] = None,
    rois: Optional[List[str]] = None,
    method_params: Optional[Dict[str, Any]] = None,
) -> Dict[int, str]:
    """
    Perform fast-fail boundary checks (method, isotope, ROI typos) before loading heavy images.
    Returns the loaded labels map.
    """
    method_key = method.lower()
    if method_key not in DOSIMETRY_STRATEGIES and method_key != "ldm":
        raise ValueError(
            f"Dosimetry method '{method}' is not implemented yet. Available: {list(DOSIMETRY_STRATEGIES.keys()) + ['ldm']}"
        )

    if not isotope:
        raise ValueError("Isotope is required.")

    if method_key == "hanscheid2017" and isotope.lower() != "lu177":
        raise ValueError(
            "The Hanscheid 2017 method (0.125 scaling factor) is specific to Lu177 isotope."
        )

    # Load label mapping
    labels_map = {}
    if labels_path is not None:
        labels_path = Path(labels_path)
        json_path = labels_path.with_suffix("").with_suffix(".json")
        if json_path.exists():
            try:
                with open(json_path, "r") as f:
                    raw_map = json.load(f)
                    labels_map = {int(k): v for k, v in raw_map.items()}
            except Exception:
                pass

    if labels_map:
        norm_to_std = {v.lower(): v for v in labels_map.values()}
        known_rois = sorted(list(labels_map.values()))

        # 1. Validate & Normalize rois in-place
        if rois:
            normalized_rois = []
            for r in rois:
                r_lower = r.lower()
                if r_lower in norm_to_std:
                    normalized_rois.append(norm_to_std[r_lower])
                else:
                    import difflib

                    close_matches = difflib.get_close_matches(
                        r, known_rois, n=5, cutoff=0.4
                    )
                    if close_matches:
                        suggestions = " or ".join(f"'{m}'" for m in close_matches)
                        raise ValueError(
                            f"Requested ROI '{r}' not found. Did you mean {suggestions}?"
                        )
                    else:
                        raise ValueError(
                            f"Requested ROI '{r}' not found. Available ROIs: {known_rois}"
                        )
            rois.clear()
            rois.extend(normalized_rois)

        # 2. Validate & Normalize method_params in-place
        if method_params:
            normalized_params = {}
            for r, val in method_params.items():
                if r.lower() == "phantom":
                    normalized_params[r] = val
                    continue
                r_lower = r.lower()
                if r_lower in norm_to_std:
                    normalized_params[norm_to_std[r_lower]] = val
                else:
                    import difflib

                    close_matches = difflib.get_close_matches(
                        r, known_rois, n=5, cutoff=0.4
                    )
                    if close_matches:
                        suggestions = " or ".join(f"'{m}'" for m in close_matches)
                        raise ValueError(
                            f"Teff ROI '{r}' not found. Did you mean {suggestions}?"
                        )
                    else:
                        raise ValueError(
                            f"Teff ROI '{r}' not found. Available ROIs: {known_rois}"
                        )
            method_params.clear()
            method_params.update(normalized_params)

    return labels_map


def compute_roi_dose(
    image: Any,
    method: str,
    isotope: str,
    labels: Optional[Any] = None,
    rois: Optional[List[str]] = None,
    ct_image: Optional[Any] = None,
    total_activity_mbq: Optional[float] = None,
    method_params: Optional[Dict[str, Any]] = None,
    image_path: Optional[Path] = None,
    labels_path: Optional[Path] = None,
    labels_map: Optional[Dict[int, str]] = None,
    verbose: bool = False,
    ct_path: Optional[Any] = None,
    teff_params: Optional[Dict[str, float]] = None,
    hu_threshold: Optional[float] = -800.0,
) -> Any:
    """
    Computes average/voxel absorbed doses for specified regions of interest.
    """
    from radot.stats import has_same_grid
    from radot.convert import convert_image_units
    from radot.meta import get_metadata_case_insensitive

    # Convert paths/strings to SimpleITK images on-the-fly
    if isinstance(image, (str, Path)):
        if image_path is None:
            image_path = Path(image)
        image = sitk.ReadImage(str(image))

    if labels is not None and isinstance(labels, (str, Path)):
        if labels_path is None:
            labels_path = Path(labels)
        labels = sitk.ReadImage(str(labels))

    if ct_image is None and ct_path is not None:
        ct_image = ct_path

    if isinstance(ct_image, (str, Path)):
        ct_image = sitk.ReadImage(str(ct_image))

    # Convert/merge teff_params into method_params structure
    if method_params is None:
        method_params = {}
    if teff_params:
        for roi, teff in teff_params.items():
            if roi not in method_params:
                method_params[roi] = {}
            method_params[roi]["teff"] = teff

    # 1. Fast-Fail Boundary Checks
    method_key = method.lower()
    if labels_map is None:
        if method_key == "ldm":
            if labels_path is not None:
                labels_map = check_dosimetry_inputs(
                    method=method,
                    isotope=isotope,
                    labels_path=labels_path,
                    rois=rois,
                    method_params=method_params,
                )
            else:
                labels_map = {}
        else:
            labels_map = check_dosimetry_inputs(
                method=method,
                isotope=isotope,
                labels_path=labels_path,
                rois=rois,
                method_params=method_params,
            )

    # 2. Strict Grid Validation (skip labels validation for ldm)
    if method_key != "ldm":
        if not has_same_grid(image, labels):
            raise ValueError(
                "Activity image grid does not match labels grid, resample image first (option --resample)."
            )
        if ct_image is not None:
            if not has_same_grid(ct_image, labels):
                raise ValueError("CT image grid does not match labels grid.")
    else:
        if ct_image is not None:
            if not has_same_grid(image, ct_image):
                raise ValueError("CT image grid does not match activity image grid. Please resample CT first.")

    # 3. Timing Extraction
    time_from_injection_h = get_time_from_injection_h(image, image_path)

    # 4. Calibrated Activity Scaling
    unit = get_metadata_case_insensitive(image, "unit", image_path)
    if unit is None:
        raise ValueError("Image unit is not specified and was not found in metadata.")

    if unit.lower() != "bq":
        logger.info(
            f"Converting image units on-the-fly from {unit} to bq for dose computation..."
        )
        image = convert_image_units(
            image, to_unit="bq", from_unit=unit, image_path=image_path
        )

    # Apply total activity override normalization if requested
    if total_activity_mbq is not None:
        stat_filter = sitk.StatisticsImageFilter()
        stat_filter.Execute(image)
        current_total_bq = stat_filter.GetSum()
        if current_total_bq > 0:
            scale_factor = (total_activity_mbq * 1e6) / current_total_bq
            logger.info(
                f"Normalizing total image activity to {total_activity_mbq} MBq (scale factor: {scale_factor:.6f})"
            )
            image = image * scale_factor

    # 5. Extract 3D data and volumes
    spacing = image.GetSpacing()
    voxel_volume_ml = (spacing[0] * spacing[1] * spacing[2]) / 1000.0

    image_arr = sitk.GetArrayFromImage(image)

    if method_key == "ldm":
        from radot.db import get_non_penetrating_energy_constant_mev
        
        # Get Delta_local in MeV
        delta_local_mev = get_non_penetrating_energy_constant_mev(isotope)
        
        # Convert MeV to Joules: 1 MeV = 1.60218e-13 Joules
        # Since mass will be in grams, we scale by 1000 to convert Joules/g to Gy (Joules/kg):
        # 1.60218e-13 * 1000 = 1.60218e-10
        mev_to_joules_per_g = 1.60218e-10
        
        # Physical decay constant lambda = ln(2) / T_p (in hours)
        t_p_h = get_physical_half_life_h(isotope)
        decay_constant_h = np.log(2) / t_p_h
        
        # Compute Time Integrated Activity (TIA) per voxel in Bq.s:
        # TIA = A_voxel_t * e^(lambda * t) / lambda_s (where lambda_s is in seconds)
        decay_constant_s = decay_constant_h / 3600.0
        tia_arr = image_arr * np.exp(decay_constant_h * time_from_injection_h) / decay_constant_s
        
        # Compute density map (g/cm3 = g/mL) from Hounsfield Units (HU)
        if ct_image is not None:
            ct_arr = sitk.GetArrayFromImage(ct_image)
            density_arr = 1.0 + (ct_arr / 1000.0)
            # Ensure density is positive
            density_arr = np.clip(density_arr, 1e-3, None)
        else:
            density_arr = np.ones_like(image_arr)
            
        # Voxel mass in grams: Volume (mL) * Density (g/mL)
        voxel_mass_g = voxel_volume_ml * density_arr
        
        # Absorbed Dose (Gy) per voxel: TIA (Bq.s) * Delta_local (MeV) * conversions / mass (g)
        dose_arr = (tia_arr * delta_local_mev * mev_to_joules_per_g) / voxel_mass_g
        
        # Apply HU threshold if CT is present (set dose to 0 for voxels with CT HU below threshold)
        if ct_image is not None and hu_threshold is not None:
            dose_arr[ct_arr < hu_threshold] = 0.0
        
        # Create SimpleITK image matching activity image space
        dose_image = sitk.GetImageFromArray(dose_arr.astype(np.float32))
        dose_image.CopyInformation(image)
        dose_image.SetMetaData("unit", "gy")
        dose_image.SetMetaData("modality", "dose")
        
        if labels is None:
            return dose_image

        labels_arr = sitk.GetArrayFromImage(labels)
        results = {}
        for label_val, roi_name in labels_map.items():
            if rois and roi_name not in rois:
                continue
            mask = (labels_arr == label_val)
            if not np.any(mask):
                continue

            roi_mass_g = np.sum(voxel_mass_g[mask])
            mean_dose_gy = np.sum(dose_arr[mask] * voxel_mass_g[mask]) / roi_mass_g if roi_mass_g > 0 else 0.0
            roi_volume_ml = np.sum(mask) * voxel_volume_ml

            results[roi_name] = {
                "mean_dose_gy": float(mean_dose_gy),
                "roi_mass_g": float(roi_mass_g),
                "roi_volume_ml": float(roi_volume_ml),
                "roi_teff_h": float(t_p_h),
            }
        return results, dose_image

    labels_arr = sitk.GetArrayFromImage(labels)
    ct_arr = sitk.GetArrayFromImage(ct_image) if ct_image is not None else None

    # Call the generic dosimetry pipeline
    return compute_dose_generic(
        image_arr=image_arr,
        labels_arr=labels_arr,
        isotope=isotope,
        voxel_volume_ml=voxel_volume_ml,
        time_from_injection_h=time_from_injection_h,
        labels_map=labels_map,
        ct_arr=ct_arr,
        rois=rois,
        method_params=method_params,
        verbose=verbose,
        method=method,
    )


def compute_roi_dose_file(
    image_path: Path,
    labels_path: Optional[Path],
    method: str,
    isotope: str,
    ct_path: Optional[Path] = None,
    total_activity_mbq: Optional[float] = None,
    rois: Optional[List[str]] = None,
    method_params: Optional[Dict[str, Any]] = None,
    output_path: Optional[Path] = None,
    resample_target: Optional[str] = None,
    verbose: bool = False,
    hu_threshold: Optional[float] = -800.0,
) -> None:
    image_path = Path(image_path)
    if labels_path is not None:
        labels_path = Path(labels_path)
    if ct_path is not None:
        ct_path = Path(ct_path)
    if output_path is not None:
        output_path = Path(output_path)

    # 1. Fast boundary/input checks before reading heavy images
    if method_params is None:
        method_params = {}

    method_key = method.lower()
    labels_map = {}
    if method_key != "ldm":
        if not labels_path:
            raise ValueError("Labels path is required for S-value dosimetry methods.")
        labels_map = check_dosimetry_inputs(
            method=method,
            isotope=isotope,
            labels_path=labels_path,
            rois=rois,
            method_params=method_params,
        )
    else:
        if labels_path is not None:
            labels_map = check_dosimetry_inputs(
                method=method,
                isotope=isotope,
                labels_path=labels_path,
                rois=rois,
                method_params=method_params,
            )
        else:
            check_dosimetry_inputs(
                method=method,
                isotope=isotope,
                labels_path=None,
                rois=None,
                method_params=None,
            )

    if resample_target == "like_ct" and not ct_path:
        raise ValueError(
            "Cannot resample like CT image because companion CT path was not provided."
        )

    if not image_path.exists():
        raise FileNotFoundError(f"Activity image not found: {image_path}")
    if labels_path is not None and not labels_path.exists():
        raise FileNotFoundError(f"Labels image not found: {labels_path}")

    logger.info(f"Reading activity image: {image_path}")
    image = sitk.ReadImage(str(image_path))

    labels = None
    if labels_path is not None:
        logger.info(f"Reading labels image: {labels_path}")
        labels = sitk.ReadImage(str(labels_path))

    ct_image = None
    if ct_path:
        if not ct_path.exists():
            raise FileNotFoundError(f"CT image not found: {ct_path}")
        logger.info(f"Reading companion CT image: {ct_path}")
        ct_image = sitk.ReadImage(str(ct_path))

    # Perform resampling if target is specified
    if resample_target:
        logger.info(f"Resampling images to align grids: target is {resample_target}")
        from radot.resampling import resample_image
        from radot.meta import get_metadata_case_insensitive

        modality = (
            get_metadata_case_insensitive(image, "modality", image_path) or "spect"
        )
        unit = get_metadata_case_insensitive(image, "unit", image_path)
        preserve_total = unit is not None and unit.lower() == "bq"

        if resample_target == "like_activity":
            if labels is not None:
                labels = resample_image(labels, image_type="labels", reference_image=image)
            if ct_image is not None:
                ct_image = resample_image(
                    ct_image,
                    image_type="ct",
                    preserve_total=False,
                    gaussian=True,
                    reference_image=image,
                )
        elif resample_target == "like_labels":
            if labels is None:
                raise ValueError("Cannot resample like labels because labels image is not provided.")
            image = resample_image(
                image,
                image_type=modality.lower(),
                preserve_total=preserve_total,
                gaussian=True,
                reference_image=labels,
            )
            if ct_image is not None:
                ct_image = resample_image(
                    ct_image,
                    image_type="ct",
                    preserve_total=False,
                    gaussian=True,
                    reference_image=labels,
                )
        elif resample_target == "like_ct":
            if ct_image is None:
                raise ValueError("Cannot resample like CT because CT image is not provided.")
            image = resample_image(
                image,
                image_type=modality.lower(),
                preserve_total=preserve_total,
                gaussian=True,
                reference_image=ct_image,
            )
            if labels is not None:
                labels = resample_image(
                    labels, image_type="labels", reference_image=ct_image
                )
        else:
            raise ValueError(f"Unknown resample target: {resample_target}")

    results = compute_roi_dose(
        image=image,
        labels=labels,
        method=method,
        isotope=isotope,
        rois=rois,
        ct_image=ct_image,
        total_activity_mbq=total_activity_mbq,
        image_path=image_path,
        labels_path=labels_path,
        verbose=verbose,
        method_params=method_params,
        hu_threshold=hu_threshold,
    )

    if method_key == "ldm":
        from radot.meta import read_metadata, write_metadata
        if labels_path is not None:
            results_dict, dose_image = results
            if output_path:
                if output_path.name.endswith((".nii", ".nii.gz", ".nrrd", ".mhd")):
                    logger.info(f"Saving dose image to: {output_path}")
                    output_path.parent.mkdir(parents=True, exist_ok=True)
                    sitk.WriteImage(dose_image, str(output_path))
                    try:
                        orig_meta = read_metadata(image_path).get("metadata", {})
                    except Exception:
                        orig_meta = {}
                    write_metadata(
                        image_path=output_path,
                        modality="dose",
                        unit="gy",
                        weight=orig_meta.get("weight"),
                        activity=orig_meta.get("activity"),
                        injection_date=orig_meta.get("injection_date"),
                        acquisition_date=orig_meta.get("acquisition_date"),
                    )
                else:
                    logger.info(f"Saving dose results to: {output_path}")
                    output_path.parent.mkdir(parents=True, exist_ok=True)
                    with open(output_path, "w") as f:
                        json.dump(results_dict, f, indent=4)
            if not output_path or output_path.name.endswith((".nii", ".nii.gz", ".nrrd", ".mhd")):
                pretty_print_dose(results_dict)
        else:
            if not output_path:
                output_path = Path("dose.nii.gz")
            logger.info(f"Saving dose image to: {output_path}")
            output_path.parent.mkdir(parents=True, exist_ok=True)
            sitk.WriteImage(results, str(output_path))
            try:
                orig_meta = read_metadata(image_path).get("metadata", {})
            except Exception:
                orig_meta = {}
            write_metadata(
                image_path=output_path,
                modality="dose",
                unit="gy",
                weight=orig_meta.get("weight"),
                activity=orig_meta.get("activity"),
                injection_date=orig_meta.get("injection_date"),
                acquisition_date=orig_meta.get("acquisition_date"),
            )
        return

    # Print results to stdout
    if output_path:
        logger.info(f"Saving dose results to: {output_path}")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(results, f, indent=4)
    else:
        # Pretty print results as a table
        pretty_print_dose(results)


def pretty_print_dose(results: dict) -> None:
    """Pretty print dosimetry results as a table."""
    import sys

    headers = ["ROI", "Mean Dose (Gy)", "Mass (g)", "Volume (mL)", "Teff (h)"]
    rows = []
    for roi, data in results.items():
        rows.append(
            [
                roi,
                f"{data['mean_dose_gy']:.3f}",
                f"{data['roi_mass_g']:.3f}",
                f"{data['roi_volume_ml']:.3f}",
                f"{data['roi_teff_h']:.1f}",
            ]
        )

    if not rows:
        sys.stdout.write("No dose computed.\n")
        return

    # Compute column widths
    col_widths = [len(h) for h in headers]
    for row in rows:
        for i, val in enumerate(row):
            col_widths[i] = max(col_widths[i], len(val))

    def format_row(row_data):
        parts = []
        for i, val in enumerate(row_data):
            width = col_widths[i]
            if i == 0:  # left aligned
                parts.append(val.ljust(width))
            else:  # right aligned
                parts.append(val.rjust(width))
        return " | ".join(parts)

    header_line = format_row(headers)
    separator = "-+-".join("-" * w for w in col_widths)

    sys.stdout.write("\n" + header_line + "\n")
    sys.stdout.write(separator + "\n")
    for row in rows:
        sys.stdout.write(format_row(row) + "\n")
    sys.stdout.write("\n")
