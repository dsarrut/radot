"""JSON dose configuration files."""

import json
import logging
from typing import Optional, Dict, Any, Tuple, Union
from pathlib import Path

logger = logging.getLogger("radot")

# Values of a dose computation given either explicitly (CLI option, Python argument) or by a dose config
# (keys of load_dose_config); method_params (per-ROI Teff, phantom) is merged per entry
DOSE_CONFIG_KEYS = (
    "timepoints",
    "image_path",
    "labels_path",
    "ct_path",
    "method",
    "isotope",
    "total_activity_mbq",
    "rois",
    "output_path",
    "resample_target",
    "verbose",
    "hu_threshold",
    "dose_rate",
    "cascade",
    "dose_rate_map_path",
    "dose_rate_map_activity_mbq",
)


def merge_dose_config(
    given: Dict[str, Any], cfg: Optional[Dict[str, Any]], given_source: str = "argument"
) -> Tuple[Dict[str, Any], Dict[str, str]]:
    """
    The one merge of the explicit values and a dose config (load_dose_config output): an explicit value
    (not None) always wins over the config. method_params entries (a ROI Teff, 'phantom') are merged per
    entry, the explicit entry winning. Returns (merged, sources): merged has every DOSE_CONFIG_KEYS key
    (None when given nowhere) and 'method_params'; sources[key] is given_source ('option' from the
    command line, 'argument' from the Python API) or 'config', for the keys with a value
    ('method_params.<entry>' for the method_params entries).
    """
    cfg = cfg or {}
    unknown = set(given) - set(DOSE_CONFIG_KEYS) - {"method_params"}
    if unknown:
        raise ValueError(f"merge_dose_config: unknown keys {sorted(unknown)}")
    merged: Dict[str, Any] = {}
    sources: Dict[str, str] = {}
    for key in DOSE_CONFIG_KEYS:
        if given.get(key) is not None:
            merged[key], sources[key] = given[key], given_source
        elif cfg.get(key) is not None:
            merged[key], sources[key] = cfg[key], "config"
        else:
            merged[key] = None
    method_params: Dict[str, Any] = {}
    for origin, entries in (("config", cfg.get("method_params")), (given_source, given.get("method_params"))):
        for name, value in (entries or {}).items():
            method_params[name] = value
            sources[f"method_params.{name}"] = origin
    merged["method_params"] = method_params
    logger.debug("Dose parameters: " + ", ".join(f"{k} ({v})" for k, v in sources.items()))
    return merged, sources


def load_dose_config(
    config: Union[str, Path, Dict[str, Any]],
) -> Dict[str, Any]:
    """
    Load and normalize a dosimetry configuration file or dictionary.
    Resolves relative file paths relative to the configuration file location.
    """
    if isinstance(config, (str, Path)):
        config_path = Path(config).resolve()
        if not config_path.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_path}")
        with open(config_path, "r") as f:
            raw_data = json.load(f)
        base_dir = config_path.parent
    elif isinstance(config, dict):
        raw_data = dict(config)
        base_dir = Path.cwd()
    else:
        raise TypeError(f"Expected path or dict for config, got {type(config).__name__}")

    def resolve_path(p: Optional[Union[str, Path]]) -> Optional[Path]:
        if p is None:
            return None
        p_obj = Path(p)
        if not p_obj.is_absolute():
            return (base_dir / p_obj).resolve()
        return p_obj.resolve()

    res: Dict[str, Any] = {}

    # Activity Image path
    img = raw_data.get("image_path") or raw_data.get("image")
    if img is not None:
        if isinstance(img, (list, tuple)):
            res["timepoints"] = [
                resolve_path(p) if isinstance(p, (str, Path)) else p for p in img
            ]
        else:
            res["image_path"] = resolve_path(img)

    # Timepoints (MTP)
    timepoints = raw_data.get("timepoints") or raw_data.get("images")
    if timepoints is not None and isinstance(timepoints, (list, tuple)):
        resolved_tps = []
        for tp in timepoints:
            if isinstance(tp, (str, Path)):
                resolved_tps.append(resolve_path(tp))
            elif isinstance(tp, dict):
                tp_dict = dict(tp)
                for k in ["image", "image_path", "path"]:
                    if k in tp_dict and tp_dict[k] is not None:
                        tp_dict[k] = resolve_path(tp_dict[k])
                for k in ["labels", "labels_path"]:
                    if k in tp_dict and tp_dict[k] is not None:
                        tp_dict[k] = resolve_path(tp_dict[k])
                for k in ["ct", "ct_path"]:
                    if k in tp_dict and tp_dict[k] is not None:
                        tp_dict[k] = resolve_path(tp_dict[k])
                resolved_tps.append(tp_dict)
            else:
                resolved_tps.append(tp)
        res["timepoints"] = resolved_tps

    # Labels path
    lbl = raw_data.get("labels_path") or raw_data.get("labels")
    if lbl is not None:
        res["labels_path"] = resolve_path(lbl)

    # CT path
    ct = raw_data.get("ct_path") or raw_data.get("ct")
    if ct is not None:
        res["ct_path"] = resolve_path(ct)

    # Dose-rate map (external Monte Carlo) and its activity
    drm = raw_data.get("dose_rate_map_path") or raw_data.get("dose_rate_map")
    if drm is not None:
        res["dose_rate_map_path"] = resolve_path(drm)
    if raw_data.get("dose_rate_map_activity_mbq") is not None:
        res["dose_rate_map_activity_mbq"] = float(raw_data["dose_rate_map_activity_mbq"])

    # Output path
    out = raw_data.get("output_path") or raw_data.get("output")
    if out is not None:
        res["output_path"] = resolve_path(out)

    # Method / Model
    method = raw_data.get("method") or raw_data.get("model")
    if method is not None:
        res["method"] = str(method)

    # Radionuclide / Isotope
    isotope = raw_data.get("isotope") or raw_data.get("radionuclide")
    if isotope is not None:
        res["isotope"] = str(isotope)

    # Total activity
    act = (
        raw_data.get("total_activity_mbq")
        if "total_activity_mbq" in raw_data
        else raw_data.get("total_activity")
    )
    if act is not None:
        res["total_activity_mbq"] = float(act)

    # Method params & Teff
    method_params: Dict[str, Any] = {}
    if "method_params" in raw_data and isinstance(raw_data["method_params"], dict):
        method_params = dict(raw_data["method_params"])

    # Phantom
    phantom = raw_data.get("phantom")
    if phantom is not None:
        res["phantom"] = str(phantom)
        if "phantom" not in method_params:
            method_params["phantom"] = str(phantom)

    # Resample target
    resample = raw_data.get("resample_target") or raw_data.get("resample")
    if resample is not None:
        res["resample_target"] = str(resample)

    # Cascade
    if "cascade" in raw_data:
        res["cascade"] = bool(raw_data["cascade"])

    # Dose rate
    if "dose_rate" in raw_data:
        res["dose_rate"] = bool(raw_data["dose_rate"])

    # HU threshold
    if "hu_threshold" in raw_data:
        res["hu_threshold"] = float(raw_data["hu_threshold"])

    # Verbose
    if "verbose" in raw_data:
        res["verbose"] = bool(raw_data["verbose"])

    # Quiet
    if "quiet" in raw_data:
        res["quiet"] = bool(raw_data["quiet"])

    # ROIs
    rois = raw_data.get("rois")
    if rois is not None:
        if isinstance(rois, str):
            rois_list = [r.strip() for r in rois.split(",") if r.strip()]
        elif isinstance(rois, dict):
            rois_list = []
            for roi_name, roi_val in rois.items():
                rois_list.append(str(roi_name))
                if isinstance(roi_val, dict):
                    if roi_name not in method_params:
                        method_params[roi_name] = {}
                    method_params[roi_name].update(roi_val)
                elif roi_val is not None:
                    method_params.setdefault(roi_name, {})["teff"] = _teff_hours(
                        roi_val, roi_name, "rois"
                    )
        elif isinstance(rois, (list, tuple)):
            rois_list = []
            for r in rois:
                if isinstance(r, str):
                    for item in r.split(","):
                        item_stripped = item.strip()
                        if item_stripped:
                            rois_list.append(item_stripped)
                elif isinstance(r, dict):
                    name = r.get("name") or r.get("roi")
                    if name:
                        rois_list.append(str(name))
                        params = {k: v for k, v in r.items() if k not in ("name", "roi")}
                        if params:
                            if name not in method_params:
                                method_params[name] = {}
                            method_params[name].update(params)
                else:
                    rois_list.append(str(r))
        else:
            rois_list = list(rois)
        res["rois"] = rois_list if rois_list else None

    teff = raw_data.get("teff")
    if teff is not None:
        if isinstance(teff, dict):
            for roi_name, val in teff.items():
                if isinstance(val, dict):
                    method_params.setdefault(roi_name, {}).update(val)
                else:
                    method_params.setdefault(roi_name, {})["teff"] = _teff_hours(
                        val, roi_name, "teff"
                    )
        else:
            items = [teff] if isinstance(teff, str) else teff
            if not isinstance(items, (list, tuple)):
                raise ValueError(
                    f"Config 'teff' must be a dict {{ROI: hours}}, or a list of 'ROI:hours', got {teff!r}."
                )
            for item in items:
                if not isinstance(item, str):
                    raise ValueError(f"Config 'teff' item {item!r} must be a 'ROI:hours' string.")
                for t in item.split(","):
                    t_strip = t.strip()
                    if not t_strip:
                        continue
                    parts = t_strip.split(":")
                    if len(parts) != 2 or not parts[0].strip():
                        raise ValueError(
                            f"Config 'teff' item '{t_strip}' must be in ROI:hours format (e.g. Liver:45.2)."
                        )
                    r_name = parts[0].strip()
                    method_params.setdefault(r_name, {})["teff"] = _teff_hours(
                        parts[1].strip(), r_name, "teff"
                    )

    if method_params:
        res["method_params"] = method_params

    return res


def _teff_hours(val, roi_name, key) -> float:
    """Teff in hours from a config value: a number, never silently dropped."""
    if isinstance(val, bool):
        val = None
    try:
        return float(val)
    except (ValueError, TypeError):
        raise ValueError(
            f"Config '{key}': the Teff of ROI '{roi_name}' must be a number of hours, got {val!r}."
        ) from None
