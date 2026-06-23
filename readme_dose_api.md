def compute_roi_dose(
    image: sitk.Image,
    labels: sitk.Image,
    method: str,
    isotope: str,
    rois: list[str] | None = None,
    ct_image: sitk.Image | None = None,
    total_activity_mbq: float | None = None,
    method_params: dict | None = None
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
    pass


radot dose spect.nii.gz \
  -l labels.nii.gz \
  --ct ct.nii.gz \
  -m madsen2018 \
  --isotope Lu177 \
  --total-activity 250.0 \
  --rois Liver,Kidney_L \
  --teff Liver:45.2,Kidney_L:32.1


Test 1

radot dose pet.nii.gz -l labels_all.nii.gz --ct ct.nii.gz -m hanscheid2017 --isotope Lu177 --total-activity 2000  --rois Liver

