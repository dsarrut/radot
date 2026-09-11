# RADOT ROI Nomenclature & Mapping Policy

This document defines the 3-stage cascading policy used by `radot` to bridge segmentations, clinical naming standards, and reference phantoms for internal dosimetry calculations.

## 1. Core Principles

* Master Anchor: AAPM TG-263 is the primary internal hub name standard.
* Lossless Fallback: If an organ or structure is missing from the TG-263 guidelines, `radot` retains the initial source name (TotalSegmentator or OpenDose) to prevent data loss and preserve anatomical context.
* Deterministic Mapping: Explicit translation dictionaries are preferred over automatic fuzzy string matching to guarantee total predictability across dose calculations.
* Data Decoupling: All nomenclature translation dictionaries must be stored as flat JSON configuration files inside the `radot/data/` folder and loaded from there dynamically.
* Architecture: ROI nomenclature resolution logic is isolated within a specific file (`radot/nomenclature.py`) and imported uniformly by other core modules.
* Interactive Guidance: Automated fuzzy mapping is banned during execution. Instead, fuzzy string distance matching is strictly reserved for user-facing interactive CLI commands (`match` and `verify`) to help users understand errors and explicitly propose corrections.

---

## 2. The 3-Stage Resolution Cascade

When the dosimetry engine looks up an S-value matrix ($S(\text{Target} \leftarrow \text{Source})$), it maps the user's segmented mask name to the OpenDose reference phantom name using the following priority order:

### Stage 1: The Standard TG-263 Route
The segmented name is translated to its official TG-263 alphanumeric code via `TS_TO_TG263`. This code is then mapped to the OpenDose string using the loaded `TG263_TO_OD` lookup.

*Example:* `kidney_left` → `Kidney_L` → `Left kidney`

### Stage 2: The Direct TS-to-OD Bypass Bridge
If a structure exists in both TotalSegmentator and OpenDose phantoms but lacks a standard or specific designation in TG-263 guidelines, the mapping completely bypasses the TG-263 layer using `TS_TO_OD_FALLBACK`.

*Example:* `portal_vein_and_splenic_vein` → `Blood vessels, trunk`

### Stage 3: Absolute Raw Retain Fallback & Safeguard
If no explicit mapping exists anywhere in the cross-walk dictionaries, the initial name is preserved exactly as it appears in the user sidecar (e.g., custom manual segmentations like `Tumor_Lesion_1`). 

Crucial Safety Rule (Strict Dosimetry Halt vs Loop Warning):
If a segmented label lacks S-value mappings in the active phantom (such as sex-specific organs like Prostate when using a female model, or unmapped bone targets):
- Generic Label Loop: When looping over all labels present in the image (no explicit --rois option specified), `radot` prints a warning to stderr and skips calculation for that specific ROI, allowing the other organs to compute successfully.
- Explicit ROI Request: When the user explicitly requests a specific ROI via the --rois parameter, `radot` halts execution immediately by raising a `click.ClickException`. To proceed, the user must supply a custom S-value override file (e.g., --custom-s-values override.json) or adjust their input sidecar mappings.

---

## 3. Reference Handling of Exceptional Regions

### Deep Musculoskeletal Elements (TS Mask Only)
Fine-grained muscle strings that are omitted by OpenDose are mapped directly to broad tissue pool categories:
* Deep back/pelvic muscles → `Muscle`

Skeletal Aggregation Prohibition:
Automatic grouping of individual bone structures (such as `rib_*` or `vertebra_*`) into a single `"Skeleton"` compartment is strictly prohibited to enforce dosimetry accuracy. OpenDose details 38 separate cortical/spongiosa sub-components and lacks a unified `"Skeleton"` S-value. If any unmapped bone fragment is encountered, the calculation is halted immediately and the user is requested to supply custom S-value overrides.

### Abstract MIRD Targets (OpenDose Only)
Structures available in OpenDose but not currently segmented by automatic segmentation tools (e.g., `Ureters`, `Thymus`) live silently inside the downloaded `.json.gz` matrix files. They act as passive targets for cross-fire deposition from active source organs, requiring no code mapping inside the segmentation layer.

---

## 4. User Tools for Nomenclature Querying & Verification

### List All Mappings
Displays the core explicit anatomical correspondences.

```bash
radot roi list --source all
```

You can narrow the overview by tracking origins and compliance checklists explicitly:

```bash
radot roi list --source ts   # Lists the TotalSegmentator Mapping Inventory checklist
radot roi list --source od   # Lists all 171 OpenDose phantom target states (Active, Passive, Inactive)
radot roi list --source tg   # Lists Supported AAPM TG-263 standard compliance labels
```

### Interactive Match (Fuzzy Suggestion Tool)
Queries a string to check for exact registration. If missing, it computes string distances to explicitly suggest correct nomenclature alternatives to the user.

```bash
radot roi match "kiney_l"
Output:

❌ "kiney_l" is not explicitly registered.
💡 Did you mean one of these?
   - Kidney_L  (TG-263)
   - kidney_left (TotalSegmentator)
```

### Sidecar Validation Pass
Scans a patient's multi-label sidecar mapping file before execution to capture and report unmapped dosimetry flags.

```bash
radot roi verify patient_001_seg.json
```
