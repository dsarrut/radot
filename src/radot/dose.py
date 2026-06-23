from typing import List, Optional, Dict, Any
import SimpleITK as sitk
from pathlib import Path

def compute_roi_dose(
    image: sitk.Image,
    labels: sitk.Image,
    method: str,
    isotope: str,
    rois: Optional[List[str]] = None,
    ct_image: Optional[sitk.Image] = None,
    total_activity_mbq: Optional[float] = None,
    method_params: Optional[Dict[str, Any]] = None
) -> dict:
    """
    Computes average/voxel absorbed doses for specified regions of interest.

    Parameters:
    -----------
    image : sitk.Image
        The functional activity map (SPECT/PET).
    labels : sitk.Image
        The anatomical segmentation label map.
    method : str
        The mathematical model to apply (e.g., 'madsen2018', 'hanscheid2017').
    isotope : str
        The therapeutic radionuclide (e.g., 'Lu177', 'Y90'). Used to load OpenDose 
        S-values and physical decay constraints.
    rois : list[str], optional
        Target labels to calculate. If None, evaluates all labels in the image.
    ct_image : sitk.Image, optional
        Companion CT scan used for HU-to-density extraction and mass scaling factors.
    total_activity_mbq : float, optional
        Global calibration override. If provided, normalizes total image activity 
        to this baseline before computing dose.
    method_params : dict, optional
        Method-specific variables parsed per-ROI (e.g., {'Liver': {'teff': 45.2}}).

    Returns:
    --------
    dict
        Structured dosimetry summary containing calculated metrics per ROI.
    """
    # Placeholder implementation
    results = {}
    
    # We will just construct a dummy response showing the inputs are mapped correctly
    target_rois = rois if rois is not None else ["Full Image"]
    for roi in target_rois:
        results[roi] = {
            "mean_dose_gy": 0.0,
            "min_dose_gy": 0.0,
            "max_dose_gy": 0.0,
            "roi_mass_g": 0.0,
            "roi_activity_mbq": 0.0,
        }
    return results


def compute_roi_dose_file(
    image_path: Path,
    labels_path: Path,
    method: str,
    isotope: str,
    ct_path: Optional[Path] = None,
    total_activity_mbq: Optional[float] = None,
    rois: Optional[List[str]] = None,
    teff_params: Optional[Dict[str, float]] = None,
    output_path: Optional[Path] = None,
) -> None:
    """
    File interface to compute ROI dose and print/save the results.
    """
    if not image_path.exists():
        raise FileNotFoundError(f"Activity image not found: {image_path}")
    if not labels_path.exists():
        raise FileNotFoundError(f"Labels image not found: {labels_path}")

    print(f"Reading activity image: {image_path}")
    image = sitk.ReadImage(str(image_path))

    print(f"Reading labels image: {labels_path}")
    labels = sitk.ReadImage(str(labels_path))

    ct_image = None
    if ct_path:
        if not ct_path.exists():
            raise FileNotFoundError(f"CT image not found: {ct_path}")
        print(f"Reading companion CT image: {ct_path}")
        ct_image = sitk.ReadImage(str(ct_path))

    # Convert teff_params into method_params structure
    method_params = {}
    if teff_params:
        for roi, teff in teff_params.items():
            method_params[roi] = {"teff": teff}

    results = compute_roi_dose(
        image=image,
        labels=labels,
        method=method,
        isotope=isotope,
        rois=rois,
        ct_image=ct_image,
        total_activity_mbq=total_activity_mbq,
        method_params=method_params
    )

    # Print results to stdout
    import json
    if output_path:
        print(f"Saving dose results to: {output_path}")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(results, f, indent=4)
    else:
        print("\n--- Dose Computation Results ---")
        print(f"Method: {method} | Isotope: {isotope}")
        print(json.dumps(results, indent=4))
        print()
