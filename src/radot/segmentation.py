import os
from pathlib import Path
import SimpleITK as sitk
from totalsegmentator.python_api import totalsegmentator

# AAPM TG-263 Official Mapping for TotalSegmentator v2
# Constructed immediately when the module is imported
# Master Mapping: "TotalSegmentator_Name": (Label_ID, "AAPM_TG263_Name")
TS_TO_TG263 = {
    # --- Glands & Organs ---
    "adrenal_gland_left": (1, "Adrenal_L"),
    "adrenal_gland_right": (2, "Adrenal_R"),
    "brain": (3, "Brain"),
    "gallbladder": (4, "Gallbladder"),
    "kidney_left": (5, "Kidney_L"),
    "kidney_right": (6, "Kidney_R"),
    "kidney_cyst_left": (7, "Kidney_L_Cyst"),
    "kidney_cyst_right": (8, "Kidney_R_Cyst"),
    "liver": (9, "Liver"),
    "pancreas": (10, "Pancreas"),
    "prostate": (11, "Prostate"),
    "spleen": (12, "Spleen"),
    "thyroid_gland": (13, "Thyroid"),
    "urinary_bladder": (14, "Bladder"),
    # --- Gastrointestinal ---
    "colon": (15, "Bowel_Lrg"),
    "small_bowel": (16, "Bowel_Sm"),
    "duodenum": (17, "Duodenum"),
    "esophagus": (18, "Esophagus"),
    "stomach": (19, "Stomach"),
    # --- Cardiovascular ---
    "heart": (20, "Heart"),
    "atrial_appendage_left": (21, "Atrium_L_Append"),
    "aorta": (22, "Aorta"),
    "inferior_vena_cava": (23, "V_Cava_I"),
    "superior_vena_cava": (24, "V_Cava_S"),
    "pulmonary_vein": (25, "V_Pulmonary"),
    "brachiocephalic_trunk": (26, "A_Brachioceph"),
    "brachiocephalic_vein_left": (27, "V_Brachioceph_L"),
    "brachiocephalic_vein_right": (28, "V_Brachioceph_R"),
    "common_carotid_artery_left": (29, "A_Carotid_C_L"),
    "common_carotid_artery_right": (30, "A_Carotid_C_R"),
    "subclavian_artery_left": (31, "A_Subclavian_L"),
    "subclavian_artery_right": (32, "A_Subclavian_R"),
    "iliac_artery_left": (33, "A_Iliac_L"),
    "iliac_artery_right": (34, "A_Iliac_R"),
    "iliac_vena_left": (35, "V_Iliac_L"),
    "iliac_vena_right": (36, "V_Iliac_R"),
    "portal_vein_and_splenic_vein": (37, "V_Portal"),
    # --- Respiratory ---
    "trachea": (38, "Trachea"),
    "lung_lower_lobe_left": (39, "Lung_L_Lower"),
    "lung_lower_lobe_right": (40, "Lung_R_Lower"),
    "lung_middle_lobe_right": (41, "Lung_R_Middle"),
    "lung_upper_lobe_left": (42, "Lung_L_Upper"),
    "lung_upper_lobe_right": (43, "Lung_R_Upper"),
    # --- Central Nervous System ---
    "spinal_cord": (44, "SpinalCord"),
    # --- Body / External ---
    "body": (45, "External"),
    "body_trunc": (46, "External"),
    "skin": (47, "Skin"),
    # --- Muscles ---
    "autochthon_left": (48, "Muscle_Back_L"),
    "autochthon_right": (49, "Muscle_Back_R"),
    "gluteus_maximus_left": (50, "Muscle_GlutMax_L"),
    "gluteus_maximus_right": (51, "Muscle_GlutMax_R"),
    "gluteus_medius_left": (52, "Muscle_GlutMed_L"),
    "gluteus_medius_right": (53, "Muscle_GlutMed_R"),
    "gluteus_minimus_left": (54, "Muscle_GlutMin_L"),
    "gluteus_minimus_right": (55, "Muscle_GlutMin_R"),
    "iliopsoas_left": (56, "Muscle_Iliopsoas_L"),
    "iliopsoas_right": (57, "Muscle_Iliopsoas_R"),
    # --- Bones ---
    "clavicula_left": (58, "Bone_Clavicle_L"),
    "clavicula_right": (59, "Bone_Clavicle_R"),
    "femur_left": (60, "Bone_Femur_L"),
    "femur_right": (61, "Bone_Femur_R"),
    "humerus_left": (62, "Bone_Humerus_L"),
    "humerus_right": (63, "Bone_Humerus_R"),
    "scapula_left": (64, "Bone_Scapula_L"),
    "scapula_right": (65, "Bone_Scapula_R"),
    "sacrum": (66, "Bone_Sacrum"),
    "skull": (67, "Bone_Skull"),
    "sternum": (68, "Bone_Sternum"),
    "costal_cartilages": (69, "Bone_Cartilage_Costal"),
    "hip_left": (138, "Bone_Ilium_L"),
    "hip_right": (139, "Bone_Ilium_R"),
    # --- Head & Neck / Glands & Cavities ---
    "eye_left": (119, "Eye_L"),
    "eye_right": (120, "Eye_R"),
    "eye_lens_left": (121, "Lens_L"),
    "eye_lens_right": (122, "Lens_R"),
    "optic_nerve_left": (123, "OpticNrv_L"),
    "optic_nerve_right": (124, "OpticNrv_R"),
    "parotid_gland_left": (125, "Parotid_L"),
    "parotid_gland_right": (126, "Parotid_R"),
    "submandibular_gland_left": (127, "Glnd_Submand_L"),
    "submandibular_gland_right": (128, "Glnd_Submand_R"),
    "nasopharynx": (129, "Pharynx_Naso"),
    "oropharynx": (130, "Pharynx_Oro"),
    "hypopharynx": (131, "Pharynx_Hypo"),
    "nasal_cavity_left": (132, "Cavity_Nasal_L"),
    "nasal_cavity_right": (133, "Cavity_Nasal_R"),
    "auditory_canal_left": (134, "Canal_Auditory_L"),
    "auditory_canal_right": (135, "Canal_Auditory_R"),
    "soft_palate": (136, "Palate_Soft"),
    "hard_palate": (137, "Palate_Hard"),
}

# --- Dynamic Ribs (IDs 70-93) ---
curr_id = 70
for i in range(1, 13):
    TS_TO_TG263[f"rib_left_{i}"] = (curr_id, f"Bone_Rib_L_{i:02d}")
    curr_id += 1
    TS_TO_TG263[f"rib_right_{i}"] = (curr_id, f"Bone_Rib_R_{i:02d}")
    curr_id += 1

# --- Dynamic Vertebrae (IDs 94-118) ---
# curr_id starts at 94
for prefix, count in [("C", 7), ("T", 12), ("L", 5), ("S", 1)]:
    for i in range(1, count + 1):
        val = (curr_id, f"Bone_V_{prefix}{i:02d}")
        TS_TO_TG263[f"vertebra__{prefix}{i}"] = val  # Support double underscore if any
        TS_TO_TG263[f"vertebra_{prefix}{i}"] = val
        TS_TO_TG263[f"vertebrae_{prefix}{i}"] = val
        curr_id += 1


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


def segment_ct(ct_path: Path, output_path: Path, mode: str = "total", fast: bool = False, test_mode: bool = False) -> None:
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
    print(f"Reading CT image from: {ct_path}")
    ct_image = sitk.ReadImage(str(ct_path))

    # Dictionary to map TotalSegmentator class IDs to names for this task
    ts_id_to_name = {}
    if mode in class_map:
        ts_id_to_name = {v: k for k, v in class_map[mode].items()}

    print(f"Running TotalSegmentator (mode: {mode}, fast: {fast}, test_mode: {test_mode}) in multi-label mode...")
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_output = Path(temp_dir) / "temp_multilabel.nii.gz"
        
        if test_mode:
            print("Test mode enabled: resampling input CT to 5.0mm spacing...")
            original_spacing = ct_image.GetSpacing()
            original_size = ct_image.GetSize()
            target_spacing = (5.0, 5.0, 5.0)
            
            new_size = [
                max(1, int(round(original_size[i] * original_spacing[i] / target_spacing[i])))
                for i in range(3)
            ]
            
            resampler = sitk.ResampleImageFilter()
            resampler.SetSize(new_size)
            resampler.SetOutputSpacing(target_spacing)
            resampler.SetOutputOrigin(ct_image.GetOrigin())
            resampler.SetOutputDirection(ct_image.GetDirection())
            resampler.SetInterpolator(sitk.sitkLinear)
            resampler.SetDefaultPixelValue(-1000)
            
            resampled_ct = resampler.Execute(ct_image)
            totalseg_input_path = Path(temp_dir) / "ct_5mm.nii.gz"
            sitk.WriteImage(resampled_ct, str(totalseg_input_path))
            totalseg_input = str(totalseg_input_path)
        else:
            totalseg_input = str(ct_path)

        # Run totalsegmentator with ml=True
        totalsegmentator(
            input=totalseg_input,
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
    fallback_id_counter = 1000

    for old_id in unique_ids:
        class_name = ts_id_to_name.get(old_id)
        if class_name:
            tg263_id = get_tg263_id(class_name)
            tg263_name = get_tg263_name(class_name)
            if tg263_id is None:
                tg263_id = fallback_id_counter
                fallback_id_counter += 1
        else:
            tg263_id = old_id
            tg263_name = f"unknown_{old_id}"

        remapped_arr[seg_arr == old_id] = tg263_id
        json_map[int(tg263_id)] = tg263_name

    # Convert back to SimpleITK Image
    merged_image = sitk.GetImageFromArray(remapped_arr)
    merged_image.CopyInformation(seg_image)

    # Verify size alignment with the original CT template image
    if merged_image.GetSize() != ct_image.GetSize():
        print("Resampling label map to match native CT coordinates...")
        resampler = sitk.ResampleImageFilter()
        resampler.SetReferenceImage(ct_image)
        resampler.SetInterpolator(sitk.sitkNearestNeighbor)
        merged_image = resampler.Execute(merged_image)

    # Write merged label map
    print(f"Saving merged label map to: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(merged_image, str(output_path))

    # Write JSON mapping file matching the output path name
    json_path = output_path.with_suffix("").with_suffix(".json")
    print(f"Saving label mapping to: {json_path}")
    # Sort JSON by ID key
    sorted_json_map = {str(k): json_map[k] for k in sorted(json_map.keys())}
    with open(json_path, "w") as f:
        json.dump(sorted_json_map, f, indent=4)

    print("Segmentation and TG-263 mapping completed successfully.")

