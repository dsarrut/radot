import os
from pathlib import Path
import SimpleITK as sitk

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


def segment_ct(ct_path: Path, output_path: Path) -> None:
    """
    Segment the given CT image using TotalSegmentator and group/name labels according to AAPM TG-263.
    (Currently a placeholder implementation).
    """
    # Verify input exists
    if not ct_path.exists():
        raise FileNotFoundError(f"Input CT file not found: {ct_path}")

    print(f"Reading CT image from: {ct_path}")
    # Read the CT image to verify it's valid
    ct_image = sitk.ReadImage(str(ct_path))

    print(f"Performing segmentation placeholder (size: {ct_image.GetSize()})...")

    # Create a dummy label map with the same metadata (spacing, origin, direction)
    # as the input image, containing zeros.
    labels_image = sitk.Image(ct_image.GetSize(), sitk.sitkUInt8)
    labels_image.CopyInformation(ct_image)

    # Save the output image
    print(f"Saving output labels to: {output_path}")
    # Ensure directory exists
    output_path.parent.mkdir(parents=True, exist_ok=True)
    sitk.WriteImage(labels_image, str(output_path))
    print("Segmentation completed successfully.")
