"""compute_roi_dose_file: dosimetry from image files (config, resampling, output)."""

import json
import logging
from typing import List, Optional, Dict, Any, Union
import SimpleITK as sitk
from pathlib import Path

from radot.dose.compute import check_dosimetry_inputs, compute_roi_dose
from radot.dose.config import load_dose_config
from radot.dose.mtp import compute_roi_dose_mtp

logger = logging.getLogger("radot")


def compute_roi_dose_file(
    image_path: Optional[Union[str, Path, List[Any]]] = None,
    labels_path: Optional[Union[str, Path]] = None,
    method: Optional[str] = None,
    isotope: Optional[str] = None,
    ct_path: Optional[Union[str, Path]] = None,
    total_activity_mbq: Optional[float] = None,
    rois: Optional[List[str]] = None,
    method_params: Optional[Dict[str, Any]] = None,
    output_path: Optional[Union[str, Path]] = None,
    resample_target: Optional[str] = None,
    verbose: bool = False,
    hu_threshold: Optional[float] = -800.0,
    dose_rate: bool = False,
    cascade: bool = False,
    config: Optional[Union[str, Path, Dict[str, Any]]] = None,
    timepoints: Optional[List[Any]] = None,
    dose_rate_map_path: Optional[Union[str, Path]] = None,
    dose_rate_map_activity_mbq: Optional[float] = None,
) -> Any:
    # 0. Check if image_path is actually a JSON config file
    if (
        image_path is not None
        and isinstance(image_path, (str, Path))
        and str(image_path).lower().endswith(".json")
        and config is None
        and method is None
        and isotope is None
    ):
        config = image_path
        image_path = None

    # Merge config file if provided
    if config is not None:
        cfg = load_dose_config(config)
        if timepoints is None:
            timepoints = cfg.get("timepoints")
        if image_path is None:
            image_path = cfg.get("image_path")
        if labels_path is None:
            labels_path = cfg.get("labels_path")
        if ct_path is None:
            ct_path = cfg.get("ct_path")
        if method is None:
            method = cfg.get("method")
        if isotope is None:
            isotope = cfg.get("isotope")
        if total_activity_mbq is None:
            total_activity_mbq = cfg.get("total_activity_mbq")
        if rois is None:
            rois = cfg.get("rois")
        if output_path is None:
            output_path = cfg.get("output_path")
        if resample_target is None:
            resample_target = cfg.get("resample_target")
        if not verbose and cfg.get("verbose"):
            verbose = True
        if hu_threshold == -800.0 and "hu_threshold" in cfg:
            hu_threshold = cfg["hu_threshold"]
        if not dose_rate and cfg.get("dose_rate"):
            dose_rate = True
        if not cascade and cfg.get("cascade"):
            cascade = True
        if dose_rate_map_path is None:
            dose_rate_map_path = cfg.get("dose_rate_map_path")
        if dose_rate_map_activity_mbq is None:
            dose_rate_map_activity_mbq = cfg.get("dose_rate_map_activity_mbq")

        cfg_mp = cfg.get("method_params", {})
        if method_params is not None:
            cfg_mp.update(method_params)
        method_params = cfg_mp

    # Multi-timepoint (MTP) dispatch
    if timepoints is not None or isinstance(image_path, (list, tuple)):
        if dose_rate_map_path is not None:
            raise ValueError("A dose-rate map cannot be used with multiple time points.")
        tp_list = timepoints if timepoints is not None else image_path
        extracted_paths = []
        extracted_labels = []
        extracted_cts = []
        extracted_times = []
        has_times = False
        for tp in tp_list:
            if isinstance(tp, (str, Path)):
                extracted_paths.append(Path(tp))
                extracted_labels.append(None)
                extracted_cts.append(None)
            elif isinstance(tp, dict):
                p = tp.get("image_path") or tp.get("image") or tp.get("path")
                if p:
                    extracted_paths.append(Path(p))
                lbl_tp = tp.get("labels_path") or tp.get("labels")
                extracted_labels.append(Path(lbl_tp) if lbl_tp else None)
                ct_tp = tp.get("ct_path") or tp.get("ct")
                extracted_cts.append(Path(ct_tp) if ct_tp else None)
                t = tp.get("time_h") if "time_h" in tp else tp.get("time")
                if t is not None:
                    has_times = True
                    extracted_times.append(float(t))
                else:
                    extracted_times.append(None)
            else:
                extracted_paths.append(tp)
                extracted_labels.append(None)
                extracted_cts.append(None)

        times_h = (
            extracted_times
            if (has_times and all(t is not None for t in extracted_times))
            else None
        )

        # Merge per-timepoint and global labels / CT:
        # If both are provided (global labels and one per tp, the one per tp is used)
        resolved_labels = [
            (lbl if lbl is not None else (Path(labels_path) if labels_path else None))
            for lbl in extracted_labels
        ]
        resolved_cts = [
            (c if c is not None else (Path(ct_path) if ct_path else None))
            for c in extracted_cts
        ]

        if not method:
            raise ValueError("Dosimetry method ('method') must be provided.")
        if not isotope:
            raise ValueError("Radionuclide ('isotope' / 'radionuclide') must be provided.")

        for idx, lbl in enumerate(resolved_labels):
            if lbl is None:
                raise ValueError(
                    f"Labels segmentation missing for timepoint {idx + 1}. "
                    "Provide a global 'labels' or specify 'labels' per timepoint."
                )

        ref_labels_path = resolved_labels[0]

        if output_path is not None:
            output_path = Path(output_path)

        if method_params is None:
            method_params = {}

        labels_map = check_dosimetry_inputs(
            method=method,
            isotope=isotope,
            labels_path=ref_labels_path,
            rois=rois,
            method_params=method_params,
        )

        for p in extracted_paths:
            if isinstance(p, Path) and not p.exists():
                raise FileNotFoundError(f"Activity image not found: {p}")
        for lbl in resolved_labels:
            if isinstance(lbl, Path) and not lbl.exists():
                raise FileNotFoundError(f"Labels image not found: {lbl}")
        for c in resolved_cts:
            if isinstance(c, Path) and not c.exists():
                raise FileNotFoundError(f"CT image not found: {c}")

        results = compute_roi_dose_mtp(
            images=extracted_paths,
            labels=resolved_labels,
            method=method,
            isotope=isotope,
            times_h=times_h,
            ct_image=resolved_cts if any(c is not None for c in resolved_cts) else None,
            rois=rois,
            labels_map=labels_map,
            labels_path=resolved_labels,
            ct_path=resolved_cts if any(c is not None for c in resolved_cts) else None,
            method_params=method_params,
            verbose=verbose,
        )

        if output_path:
            logger.info(f"Saving dose results to: {output_path}")
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, "w") as f:
                json.dump(results, f, indent=4)
        else:
            pretty_print_dose(results)

        return results

    if not image_path:
        raise ValueError("Activity image path ('image') must be provided.")
    if not method:
        raise ValueError("Dosimetry method ('method') must be provided.")
    if not isotope:
        raise ValueError("Radionuclide ('isotope' / 'radionuclide') must be provided.")

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
    if dose_rate_map_path is not None:
        dose_rate_map_path = Path(dose_rate_map_path)
        if not dose_rate_map_path.exists():
            raise FileNotFoundError(f"Dose-rate map not found: {dose_rate_map_path}")
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
        # Per-voxel quantities (Bq, Bq.s) keep their total when resampled
        preserve_total = unit is not None and unit.lower() in ("bq", "bqs")

        if resample_target == "like_activity":
            if labels is not None:
                labels = resample_image(
                    labels, image_type="labels", reference_image=image
                )
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
                raise ValueError(
                    "Cannot resample like labels because labels image is not provided."
                )
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
                raise ValueError(
                    "Cannot resample like CT because CT image is not provided."
                )
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
        dose_rate=dose_rate,
        cascade=cascade,
        dose_rate_map_path=dose_rate_map_path,
        dose_rate_map_activity_mbq=dose_rate_map_activity_mbq,
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
                        unit="gy/h" if dose_rate else "gy",
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
            if not output_path or output_path.name.endswith(
                (".nii", ".nii.gz", ".nrrd", ".mhd")
            ):
                pretty_print_dose(results_dict)
            return results
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
                unit="gy/h" if dose_rate else "gy",
                weight=orig_meta.get("weight"),
                activity=orig_meta.get("activity"),
                injection_date=orig_meta.get("injection_date"),
                acquisition_date=orig_meta.get("acquisition_date"),
            )
            return results

    # Print results to stdout
    if output_path:
        logger.info(f"Saving dose results to: {output_path}")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(results, f, indent=4)
    else:
        # Pretty print results as a table
        pretty_print_dose(results)

    return results


def pretty_print_dose(results: dict) -> None:
    """Pretty print dosimetry results as a table."""
    import sys

    is_rate = any("mean_dose_rate_gyh" in data for data in results.values())
    has_tia = any("tia_mbq_h" in data for data in results.values())
    if is_rate:
        headers = ["ROI", "Mean Dose Rate (Gy/h)", "Mass (g)", "Volume (mL)"]
    elif has_tia:
        headers = [
            "ROI",
            "Mean Dose (Gy)",
            "Mass (g)",
            "Volume (mL)",
            "Teff (h)",
            "TIA (MBq.h)",
        ]
    else:
        headers = ["ROI", "Mean Dose (Gy)", "Mass (g)", "Volume (mL)", "Teff (h)"]

    rows = []
    for roi, data in results.items():
        if is_rate:
            mean_val = data.get("mean_dose_rate_gyh", 0.0)
            rows.append(
                [
                    roi,
                    f"{mean_val:.3f}",
                    f"{data['roi_mass_g']:.3f}",
                    f"{data['roi_volume_ml']:.3f}",
                ]
            )
        elif has_tia:
            rows.append(
                [
                    roi,
                    f"{data['mean_dose_gy']:.3f}",
                    f"{data['roi_mass_g']:.3f}",
                    f"{data['roi_volume_ml']:.3f}",
                    f"{data['roi_teff_h']:.1f}" if "roi_teff_h" in data else "-",
                    f"{data.get('tia_mbq_h', 0.0):.2f}",
                ]
            )
        else:
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
