from pathlib import Path
from typing import Optional, Dict, Union
import logging
import json
from datetime import datetime

logger = logging.getLogger("radot")
import SimpleITK as sitk


def calculate_lbm(
    weight: float,
    height: Optional[float] = None,
    sex: Optional[str] = None,
    method: str = "james",
    labels: Optional[sitk.Image] = None,
    labels_map: Optional[Dict[int, str]] = None,
) -> float:
    """
    Calculate Lean Body Mass (LBM) in kg using different strategies.
    """
    method = method.lower()
    if method == "segmentation":
        if labels is None:
            raise ValueError(
                "Labels image must be provided for segmentation-based LBM calculation."
            )
        if weight is None or weight <= 0:
            raise ValueError("Weight must be positive.")
        if not labels_map:
            raise ValueError(
                "A labels_map mapping label IDs to names is required to identify fat/adipose tissue classes."
            )

        spacing = labels.GetSpacing()
        voxel_volume_ml = (spacing[0] * spacing[1] * spacing[2]) / 1000.0
        labels_arr = sitk.GetArrayViewFromImage(labels)

        # Identify label IDs matching fat or adipose tissue
        fat_label_ids = []
        for label_id, name in labels_map.items():
            name_lower = name.lower().strip()
            # Split name by common word separators
            tokens = name_lower.replace("_", " ").replace("-", " ").split()
            if (
                "fat" in name_lower
                or "adipose" in name_lower
                or "lipoma" in name_lower
                or any(t in ("sat", "vat", "tat") for t in tokens)
            ):
                fat_label_ids.append(label_id)

        if not fat_label_ids:
            logger.warning(
                "No fat or adipose tissue classes found in the provided labels_map."
            )
            return weight  # No fat tissue means LBM = total weight

        # Sum voxels of fat
        fat_voxels = 0
        for fid in fat_label_ids:
            fat_voxels += (labels_arr == fid).sum()

        fat_volume_ml = fat_voxels * voxel_volume_ml
        # Fat Mass in kg: Volume (mL) * Density (0.92 g/mL) * 1e-3 (kg/g)
        fat_mass_kg = fat_volume_ml * 0.92 * 1e-3

        # LBM = Weight - Fat Mass
        lbm = weight - fat_mass_kg
        return float(max(0.01, min(weight, lbm)))

    # Analytical formulas require height and sex
    if height is None or height <= 0:
        raise ValueError(
            f"A positive height (in cm) is required for {method} LBM calculation."
        )
    if sex is None or sex.lower().strip() not in ("male", "female"):
        raise ValueError(
            f"Sex ('male' or 'female') is required for {method} LBM calculation."
        )

    norm_sex = sex.lower().strip()
    if method == "james":
        if norm_sex == "male":
            return 1.10 * weight - 128.0 * ((weight / height) ** 2)
        else:
            return 1.07 * weight - 148.0 * ((weight / height) ** 2)
    elif method == "janmahasatian":
        bmi = weight / ((height / 100.0) ** 2)
        if norm_sex == "male":
            return (9.27 * 1000.0 * weight) / (6.68 * 1000.0 + 216.0 * bmi)
        else:
            return (9.27 * 1000.0 * weight) / (8.78 * 1000.0 + 244.0 * bmi)
    else:
        raise ValueError(f"Unsupported LBM method: {method}")


def compute_decay_correction(
    injection_date: Optional[Union[str, datetime]] = None,
    acquisition_date: Optional[Union[str, datetime]] = None,
    radionuclide: Optional[str] = None,
    decay_factor: Optional[float] = None,
) -> tuple[float, Optional[float]]:
    """
    Computes physical radioactive decay correction factor:
        factor = 2.0 ** (delay_minutes / half_life_minutes)
    such that:
        A_scan = A_inj / factor
    Returns (factor, delay_minutes). If dates or radionuclide cannot be parsed, returns (1.0, None).
    """

    def _parse_dt(val):
        if isinstance(val, datetime):
            return val
        if not val:
            return None
        val_str = str(val).strip()
        for fmt in (
            "%Y-%m-%d %H:%M:%S",
            "%Y-%m-%d %H:%M:%S.%f",
            "%Y-%m-%dT%H:%M:%S",
            "%Y-%m-%dT%H:%M:%S.%f",
            "%Y/%m/%d, %H:%M:%S",
            "%Y/%m/%d, %H:%M:%S.%f",
            "%Y/%m/%d %H:%M:%S",
            "%d/%m/%Y %H:%M:%S",
            "%Y%m%d%H%M%S",
            "%Y-%m-%d",
        ):
            try:
                return datetime.strptime(val_str, fmt)
            except ValueError:
                continue
        return None

    delay_minutes = None
    if injection_date and acquisition_date:
        t_inj = _parse_dt(injection_date)
        t_acq = _parse_dt(acquisition_date)
        if t_inj is not None and t_acq is not None:
            delay_minutes = (t_acq - t_inj).total_seconds() / 60.0

    if decay_factor is not None and float(decay_factor) > 0:
        return float(decay_factor), (
            float(delay_minutes)
            if delay_minutes is not None and delay_minutes >= 0
            else None
        )

    if not injection_date or not acquisition_date or not radionuclide:
        return 1.0, None

    tracer_clean = str(radionuclide).strip().lower()
    tracer_map = {
        "psma": "Ga-68",
        "ga68": "Ga-68",
        "ga-68": "Ga-68",
        "fdg": "F-18",
        "f18": "F-18",
        "f-18": "F-18",
        "dotatate": "Ga-68",
    }
    isotope = tracer_map.get(tracer_clean, str(radionuclide).strip())

    try:
        from radot.db import get_physical_half_life_h

        half_life_minutes = get_physical_half_life_h(isotope) * 60.0
    except Exception as e:
        logger.warning(
            f"Could not resolve half-life for radionuclide '{radionuclide}': {e}"
        )
        return 1.0, None

    if delay_minutes is None:
        t_inj = _parse_dt(injection_date)
        t_acq = _parse_dt(acquisition_date)
        if t_inj is None or t_acq is None:
            return 1.0, None
        delay_minutes = (t_acq - t_inj).total_seconds() / 60.0

    if delay_minutes < 0:
        logger.warning(
            f"Negative scan delay ({delay_minutes:.1f} min): acquisition ({t_acq}) before injection ({t_inj})."
        )
        return 1.0, None

    factor = 2.0 ** (delay_minutes / half_life_minutes)
    return float(factor), float(delay_minutes)


def convert_image_units(
    image: sitk.Image,
    to_unit: str,
    from_unit: Optional[str] = None,
    weight: Optional[float] = None,
    activity: Optional[float] = None,
    height: Optional[float] = None,
    sex: Optional[str] = None,
    lbm_method: str = "james",
    labels: Optional[sitk.Image] = None,
    labels_map: Optional[dict] = None,
    image_path: Optional[Path] = None,
    radionuclide: Optional[str] = None,
    injection_date: Optional[Union[str, datetime]] = None,
    acquisition_date: Optional[Union[str, datetime]] = None,
    decay_factor: Optional[float] = None,
) -> sitk.Image:
    """
    Convert image voxel values between units (bq, bqml, suv, sul).
    Supports automatic decay correction when converting bqml to suv/sul.
    """
    from radot.meta import get_metadata_case_insensitive

    to_unit = to_unit.lower()

    # 1. Resolve from_unit
    if from_unit is None:
        from_unit = get_metadata_case_insensitive(image, "unit", image_path)
        if from_unit is None:
            raise ValueError(
                "Source unit is not specified and was not found in the image metadata."
            )
    else:
        from_unit = from_unit.lower()

    if from_unit == to_unit:
        return image

    # 2. Resolve weight and activity if SUV/SUL is involved
    if any(u in (from_unit, to_unit) for u in ("suv", "sul")):
        if weight is None:
            w_str = get_metadata_case_insensitive(image, "weight", image_path)
            if w_str is not None:
                try:
                    weight = float(w_str)
                except ValueError:
                    raise ValueError(
                        f"Invalid weight format in image metadata: {w_str}"
                    )
            else:
                raise ValueError(
                    "A positive --weight (in kg) is required for SUV/SUL conversion (not found in image metadata)."
                )
        if weight <= 0:
            raise ValueError("Weight must be positive.")

        if activity is None:
            act_str = get_metadata_case_insensitive(image, "activity", image_path)
            if act_str is not None:
                try:
                    activity = float(act_str)
                except ValueError:
                    raise ValueError(
                        f"Invalid activity format in image metadata: {act_str}"
                    )
            else:
                raise ValueError(
                    "A positive --activity (in MBq) is required for SUV/SUL conversion (not found in image metadata)."
                )
        if activity <= 0:
            raise ValueError("Activity must be positive.")

    # 3. Resolve decay correction parameters if converting bqml to SUV/SUL
    if decay_factor is None:
        if radionuclide is None:
            radionuclide = get_metadata_case_insensitive(
                image, "radionuclide", image_path
            ) or get_metadata_case_insensitive(image, "isotope", image_path)
        if injection_date is None:
            injection_date = get_metadata_case_insensitive(
                image, "injection_date", image_path
            )
        if acquisition_date is None:
            acquisition_date = get_metadata_case_insensitive(
                image, "acquisition_date", image_path
            )

    decay_corr_factor, _ = compute_decay_correction(
        injection_date=injection_date,
        acquisition_date=acquisition_date,
        radionuclide=radionuclide,
        decay_factor=decay_factor,
    )

    effective_activity = activity
    if (
        from_unit in ("bq", "bqml")
        and to_unit in ("suv", "sul")
        and activity is not None
        and decay_corr_factor > 0
    ):
        effective_activity = activity / decay_corr_factor
    elif (
        from_unit in ("suv", "sul")
        and to_unit in ("bq", "bqml")
        and activity is not None
        and decay_corr_factor > 0
    ):
        effective_activity = activity / decay_corr_factor

    # 4. Resolve height and sex if SUL is involved (and not using segmentation-based LBM)
    lbm = None
    if any(u in (from_unit, to_unit) for u in ("sul",)):
        lbm_method_key = lbm_method.lower()
        if lbm_method_key != "segmentation":
            if height is None:
                h_str = get_metadata_case_insensitive(image, "height", image_path)
                if h_str is not None:
                    try:
                        height = float(h_str)
                    except ValueError:
                        raise ValueError(
                            f"Invalid height format in image metadata: {h_str}"
                        )
                else:
                    raise ValueError(
                        f"A positive --height (in cm) is required for {lbm_method} SUL conversion (not found in image metadata)."
                    )
            if height <= 0:
                raise ValueError("Height must be positive.")

            if sex is None:
                sex_str = get_metadata_case_insensitive(image, "sex", image_path)
                if sex_str is not None:
                    sex = sex_str.lower().strip()
                else:
                    raise ValueError(
                        f"Patient sex ('male' or 'female') is required for {lbm_method} SUL conversion (not found in image metadata)."
                    )
            if sex.lower().strip() not in ("male", "female"):
                raise ValueError(
                    f"Invalid sex value: {sex}. Expected 'male' or 'female'."
                )

        lbm = calculate_lbm(
            weight=weight,
            height=height,
            sex=sex,
            method=lbm_method,
            labels=labels,
            labels_map=labels_map,
        )

    # Calculate voxel volume in mL (1 mL = 1000 mm^3)
    spacing = image.GetSpacing()
    voxel_volume_ml = (spacing[0] * spacing[1] * spacing[2]) / 1000.0

    # Step 1: Convert input to bqml (Becquerels per mL)
    converted = sitk.Cast(image, sitk.sitkFloat32)

    if from_unit == "bq":
        converted = converted / voxel_volume_ml
    elif from_unit == "suv":
        activity_bq = effective_activity * 1_000_000.0
        weight_g = weight * 1000.0
        suv_factor = activity_bq / weight_g
        converted = converted * suv_factor
    elif from_unit == "sul":
        activity_bq = effective_activity * 1_000_000.0
        lbm_g = lbm * 1000.0
        sul_factor = activity_bq / lbm_g
        converted = converted * sul_factor
    elif from_unit != "bqml":
        raise ValueError(f"Unsupported source unit: {from_unit}")

    # Step 2: Convert bqml to target unit
    if to_unit == "bq":
        converted = converted * voxel_volume_ml
    elif to_unit == "suv":
        activity_bq = effective_activity * 1_000_000.0
        weight_g = weight * 1000.0
        suv_factor = activity_bq / weight_g
        converted = converted / suv_factor
    elif to_unit == "sul":
        activity_bq = effective_activity * 1_000_000.0
        lbm_g = lbm * 1000.0
        sul_factor = activity_bq / lbm_g
        converted = converted / sul_factor
    elif to_unit != "bqml":
        raise ValueError(f"Unsupported target unit: {to_unit}")

    return sitk.Cast(converted, sitk.sitkFloat32)


def convert_image_units_file(
    image_path: Path,
    to_unit: str,
    from_unit: Optional[str] = None,
    output_path: Optional[Path] = None,
    weight: Optional[float] = None,
    activity: Optional[float] = None,
    height: Optional[float] = None,
    sex: Optional[str] = None,
    lbm_method: str = "james",
    labels_path: Optional[Path] = None,
    radionuclide: Optional[str] = None,
    injection_date: Optional[Union[str, datetime]] = None,
    acquisition_date: Optional[Union[str, datetime]] = None,
    decay_factor: Optional[float] = None,
) -> None:
    """
    Load image from path, convert its units, and save to output path.
    """
    if not image_path.exists():
        raise FileNotFoundError(f"Input image not found: {image_path}")

    to_unit = to_unit.lower()

    if output_path is None:
        suffix = "".join(image_path.suffixes)
        base = image_path.name[: -len(suffix)] if suffix else image_path.name
        output_path = image_path.parent / f"{base}_{to_unit}{suffix}"

    logger.info(f"Reading image from: {image_path}")
    image = sitk.ReadImage(str(image_path))

    # Resolve labels and labels_map if method is segmentation
    labels = None
    labels_map = None
    if lbm_method.lower() == "segmentation" and labels_path is not None:
        labels_path = Path(labels_path)
        logger.info(f"Reading labels image for LBM: {labels_path}")
        labels = sitk.ReadImage(str(labels_path))

        json_path = labels_path.with_suffix("").with_suffix(".json")
        if json_path.exists():
            try:
                with open(json_path, "r") as f:
                    raw_map = json.load(f)
                    labels_map = {int(k): v for k, v in raw_map.items()}
            except Exception as e:
                logger.warning(
                    f"Failed to read labels map sidecar JSON {json_path}: {e}"
                )

    # Resolve from_unit if None for print info
    resolved_from_unit = from_unit
    if resolved_from_unit is None:
        from radot.meta import get_metadata_case_insensitive

        resolved_from_unit = get_metadata_case_insensitive(image, "unit", image_path)
        if resolved_from_unit is None:
            raise ValueError(
                "Source unit is not specified and was not found in the image metadata."
            )

    logger.info(f"Converting units from {resolved_from_unit} to {to_unit}...")
    converted_image = convert_image_units(
        image=image,
        from_unit=from_unit,
        to_unit=to_unit,
        weight=weight,
        activity=activity,
        height=height,
        sex=sex,
        lbm_method=lbm_method,
        labels=labels,
        labels_map=labels_map,
        image_path=image_path,
        radionuclide=radionuclide,
        injection_date=injection_date,
        acquisition_date=acquisition_date,
        decay_factor=decay_factor,
    )

    logger.info(f"Saving converted image to: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(converted_image, str(output_path))

    # Write updated metadata to output sidecar JSON
    from radot.meta import (
        get_sidecar_path,
        get_metadata_case_insensitive,
        compute_file_sha256,
    )

    modality = get_metadata_case_insensitive(image, "modality", image_path)
    weight_val = get_metadata_case_insensitive(image, "weight", image_path)
    activity_val = get_metadata_case_insensitive(image, "activity", image_path)
    inj_val = get_metadata_case_insensitive(image, "injection_date", image_path)
    acq_val = get_metadata_case_insensitive(image, "acquisition_date", image_path)
    height_val = get_metadata_case_insensitive(image, "height", image_path)
    sex_val = get_metadata_case_insensitive(image, "sex", image_path)
    radionuclide_val = get_metadata_case_insensitive(
        image, "radionuclide", image_path
    ) or get_metadata_case_insensitive(image, "isotope", image_path)

    new_meta = {}
    if modality:
        new_meta["modality"] = modality
    new_meta["unit"] = to_unit

    # Merge existing values with CLI parameters
    final_weight = (
        weight if weight is not None else (float(weight_val) if weight_val else None)
    )
    if final_weight:
        new_meta["weight"] = final_weight

    final_activity = (
        activity
        if activity is not None
        else (float(activity_val) if activity_val else None)
    )
    if final_activity:
        new_meta["activity"] = final_activity

    final_injection_date = injection_date if injection_date is not None else inj_val
    if final_injection_date:
        new_meta["injection_date"] = (
            final_injection_date.strftime("%Y-%m-%d %H:%M:%S")
            if isinstance(final_injection_date, datetime)
            else str(final_injection_date)
        )

    final_acquisition_date = (
        acquisition_date if acquisition_date is not None else acq_val
    )
    if final_acquisition_date:
        new_meta["acquisition_date"] = (
            final_acquisition_date.strftime("%Y-%m-%d %H:%M:%S")
            if isinstance(final_acquisition_date, datetime)
            else str(final_acquisition_date)
        )

    final_height = (
        height if height is not None else (float(height_val) if height_val else None)
    )
    if final_height:
        new_meta["height"] = final_height

    final_sex = sex if sex is not None else (sex_val if sex_val else None)
    if final_sex:
        new_meta["sex"] = final_sex

    final_radionuclide = radionuclide if radionuclide is not None else radionuclide_val
    if final_radionuclide:
        new_meta["radionuclide"] = str(final_radionuclide)

    # Compute decay correction parameters for metadata
    decay_corr_factor, decay_delay_mins = compute_decay_correction(
        injection_date=final_injection_date,
        acquisition_date=final_acquisition_date,
        radionuclide=final_radionuclide,
        decay_factor=decay_factor,
    )

    if (
        decay_factor is not None
        or decay_delay_mins is not None
        or (decay_corr_factor != 1.0 and to_unit in ("suv", "sul"))
    ):
        new_meta["decay_correction_factor"] = float(decay_corr_factor)
        if decay_delay_mins is not None:
            new_meta["decay_delay_minutes"] = float(decay_delay_mins)
        new_meta["decay_correction_reference"] = "injection_time"

    sha = compute_file_sha256(output_path)
    if sha:
        new_meta["sha256"] = sha

    out_sidecar = get_sidecar_path(output_path)
    with open(out_sidecar, "w") as f:
        json.dump(new_meta, f, indent=4)

    logger.info("Unit conversion completed successfully.")
