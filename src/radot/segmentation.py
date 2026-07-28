import os
from pathlib import Path
import logging

logger = logging.getLogger("radot")
import SimpleITK as sitk
from totalsegmentator.python_api import totalsegmentator

from radot.nomenclature import TS_TO_TG263, FALLBACK_UNMAPPED_ID_RANGE

def get_tg263_name(filename):
    p = Path(filename)
    while p.suffix:
        p = p.with_suffix("")

    name = p.name
    if name.endswith("_crop"):
        name = name[:-5]

    # Safely retrieve TG-263 or fallback to TS name
    res = TS_TO_TG263.get(name)
    if res is not None:
        return res[1]
    return name


def get_tg263_id(filename):
    p = Path(filename)
    while p.suffix:
        p = p.with_suffix("")

    name = p.name
    if name.endswith("_crop"):
        name = name[:-5]

    res = TS_TO_TG263.get(name)
    if res is not None:
        return res[0]
    return None


def segment_ct(ct_path: Path, output_path: Path, mode: str = "total", fast: bool = False) -> None:
    """
    Segment the given CT image using TotalSegmentator and group/name labels according to AAPM TG-263.
    Uses the multi-label (ml=True) output option to avoid writing multiple individual mask files.
    """
    import tempfile
    import json
    import numpy as np
    from totalsegmentator.map_to_binary import class_map

    if not ct_path.exists():
        raise FileNotFoundError(f"Input CT file not found: {ct_path}")

    # Read native CT template
    logger.info(f"Reading CT image from: {ct_path}")
    ct_image = sitk.ReadImage(str(ct_path))

    # Dictionary to map TotalSegmentator class IDs to names for this task
    ts_id_to_name = {}
    if mode in class_map:
        ts_id_to_name = class_map[mode]

    # Construct equivalent CLI command

    ts_cli = f"totalsegmentator -i {ct_path} -o {output_path} --ml"
    if mode != "total":
        ts_cli += f" --task {mode}"
    if fast:
        ts_cli += " --fast"
    logger.info(f"Equivalent TotalSegmentator command: {ts_cli}")

    logger.info(f"Running TotalSegmentator (mode: {mode}, fast: {fast}) in multi-label mode...")

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_output = Path(temp_dir) / "temp_multilabel.nii.gz"
        
        # Run totalsegmentator with ml=True
        totalsegmentator(
            input=str(ct_path),
            output=str(temp_output),
            task=mode,
            fast=fast,
            ml=True,
            quiet=False
        )

        if not temp_output.exists():
            raise RuntimeError("TotalSegmentator did not produce the multi-label file.")



        # Read the generated multi-label image
        seg_image = sitk.ReadImage(str(temp_output))
        seg_arr = sitk.GetArrayFromImage(seg_image)

    # Find the unique labels predicted (excluding background)
    unique_ids = np.unique(seg_arr)
    unique_ids = [int(i) for i in unique_ids if i > 0]

    # Remap the labels to the AAPM TG-263 system
    remapped_arr = np.zeros_like(seg_arr, dtype=np.uint16)
    json_map = {}
    fallback_lo, fallback_hi = FALLBACK_UNMAPPED_ID_RANGE
    fallback_id_counter = fallback_lo
    unmapped_class_names = []

    for old_id in unique_ids:
        class_name = ts_id_to_name.get(old_id)
        if class_name:
            tg263_id = get_tg263_id(class_name)
            tg263_name = get_tg263_name(class_name)
            if tg263_id is None:
                # Not yet present in TS_TO_TG263: assign an ID from the reserved
                # fallback range instead of colliding with DISCRETE_TISSUE_IDS.
                if fallback_id_counter > fallback_hi:
                    raise RuntimeError(
                        f"Exhausted FALLBACK_UNMAPPED_ID_RANGE {FALLBACK_UNMAPPED_ID_RANGE} "
                        f"while assigning an ID for unmapped TS class '{class_name}'."
                    )
                tg263_id = fallback_id_counter
                tg263_name = class_name
                fallback_id_counter += 1
                unmapped_class_names.append(class_name)
        else:
            tg263_id = old_id
            tg263_name = f"unknown_{old_id}"

        remapped_arr[seg_arr == old_id] = tg263_id
        json_map[int(tg263_id)] = tg263_name

    if unmapped_class_names:
        logger.warning(
            "TS_TO_TG263 coverage gap: the following TotalSegmentator classes "
            "are not yet mapped and were assigned fallback IDs from "
            f"{FALLBACK_UNMAPPED_ID_RANGE}: {unmapped_class_names}"
        )

    # Convert back to SimpleITK Image
    merged_image = sitk.GetImageFromArray(remapped_arr)
    merged_image.CopyInformation(seg_image)

    # Verify size alignment with the original CT template image
    if merged_image.GetSize() != ct_image.GetSize():
        logger.info("Resampling label map to match native CT coordinates...")
        resampler = sitk.ResampleImageFilter()
        resampler.SetReferenceImage(ct_image)
        resampler.SetInterpolator(sitk.sitkNearestNeighbor)
        merged_image = resampler.Execute(merged_image)

    # Write merged label map
    logger.info(f"Saving merged label map to: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(merged_image, str(output_path))

    # Write JSON mapping file matching the output path name
    json_path = output_path.with_suffix("").with_suffix(".json")
    logger.info(f"Saving label mapping to: {json_path}")
    # Sort JSON by ID key
    sorted_json_map = {str(k): json_map[k] for k in sorted(json_map.keys())}
    with open(json_path, "w") as f:
        json.dump(sorted_json_map, f, indent=4)

    logger.info("Segmentation and TG-263 mapping completed successfully.")

