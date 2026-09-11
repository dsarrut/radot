# SUL Conversion (Lean Body Mass Scaling)

This document describes how SUL (Standardized Uptake Value normalized by Lean Body Mass) is programmatically calculated within `radot` using stoichiometric integration of segmented CT images, rendering the calibration independent of patient weight (from dicom or a database).

---

## 1. Mathematical Derivation

Standard SUV (Standardized Uptake Value) scales voxel avidity by a patient's total body weight:

$$\text{SUV}_{\text{voxel}} = \frac{\text{Activity}_{\text{voxel}} \times \text{Weight}_{\text{patient}}}{\text{Injected\_Activity}}$$

SUL scales avidity by **Lean Body Mass (LBM)**:

$$\text{SUL}_{\text{voxel}} = \frac{\text{Activity}_{\text{voxel}} \times \text{LBM}}{\text{Injected\_Activity}}$$

Within `radot`, instead of relying on clinical LBM equations (like James or Janmahasatian, which require manual height and weight entries that are often missing or incorrect), LBM is calculated directly from the CT volume using a segmented density integration.

---

## 2. Weight Cancellation Property

By definition:

$$\text{SUL}_{\text{voxel}} = \text{SUV}_{\text{voxel}} \times \text{Ratio}_{\text{LBM}}$$

Where $\text{Ratio}_{\text{LBM}}$ is the lean fraction of the scanned patient body:

$$\text{Ratio}_{\text{LBM}} = \frac{\text{LBM}_{\text{scanned}}}{\text{Weight}_{\text{scanned}}}$$

Expanding $\text{SUV}_{\text{voxel}}$ in the equation:

$$\text{SUL}_{\text{voxel}} = \left( \frac{\text{Activity}_{\text{voxel}} \times \text{Weight}_{\text{patient}}}{\text{Injected\_Activity}} \right) \times \frac{\text{LBM}_{\text{scanned}}}{\text{Weight}_{\text{scanned}}}$$

Assuming the scanned weight of the patient ($\text{Weight}_{\text{scanned}}$) matches their actual physical weight ($\text{Weight}_{\text{patient}}$), the weight parameters cancel out completely:

$$\text{SUL}_{\text{voxel}} = \frac{\text{Activity}_{\text{voxel}} \times \text{LBM}_{\text{scanned}}}{\text{Injected\_Activity}}$$

### Clinical Advantage:
Because weight cancels out, **SUL is mathematically independent of the manually recorded patient weight value**. Any errors, typos, or missing weight parameters in `db.json` or DICOM headers will not affect SUL calculations.

---

## 3. LBM Estimation via CT Density Integration

`radot` integrates mass voxel-by-voxel within the segmented body boundary.

### 3.1 Stoichiometric Density Mapping
The Hounsfield Unit (HU) of each voxel is converted to physical density ($\rho$ in $\text{g/cm}^3$) using a multi-linear calibration function:

$$\rho(\text{HU}) = \begin{cases} 
1.000 + 0.00100 \times \text{HU} & \text{for } -1000 \le \text{HU} \le 0 \quad \text{(Air to Water)} \\
1.000 + 0.00092 \times \text{HU} & \text{for } 0 < \text{HU} \le 100 \quad \text{(Soft Tissue / Fat)} \\
1.018 + 0.00059 \times \text{HU} & \text{for } \text{HU} > 100 \quad \text{(Skeletal / Bone)} 
\end{cases}$$

### 3.2 Mass Integrals
Using the TotalSegmentator label map:
* **Fat Mass ($M_{\text{fat}}$)**: Summed over voxels belonging to fat classes (e.g., `subcutaneous_fat`, `torso_fat`).
* **Lean Mass ($\text{LBM}_{\text{scanned}}$)**: Summed over all remaining non-fat anatomical classes inside the body outline.

$$\text{LBM}_{\text{scanned}} = V_{\text{voxel}} \times \sum_{i \in \text{Lean}} \rho(\text{HU}_i)$$
$$\text{Weight}_{\text{scanned}} = \text{LBM}_{\text{scanned}} + M_{\text{fat}}$$

The final ratio $\text{Ratio}_{\text{LBM}} = \frac{\text{LBM}_{\text{scanned}}}{\text{Weight}_{\text{scanned}}}$ is applied to the input SUV image to produce the SUL volume.

---

## 4. Bibliography

1. **Stoichiometric Conversion Calibration**:
   * Schneider W, Bortfeld T, Schlegel W. *Correlation between CT numbers and tissue parameters for Monte Carlo dose calculations in radiotherapy.* Phys Med Biol. 2000;45(2):459-478. [doi:10.1088/0031-9155/45/2/312](https://doi.org/10.1088/0031-9155/45/2/312)

2. **LBM/SUL Clinical Standard**:
   * Wahl RL, Jacene A, Kasamon Y, Lodge MA. *From RECIST to PERCIST: Evolving Considerations for PET Response Criteria in Solid Tumors.* J Nucl Med. 2009;50 Suppl 1(Suppl 1):122S-50S. [doi:10.2967/jnumed.108.057307](https://doi.org/10.2967/jnumed.108.057307)
   * Janmahasatian S, Duffull SB, Ash S, Ward LC, Byrne NM, Green B. *Quantification of lean bodyweight.* Clin Pharmacokinet. 2005;44(10):1051-1065. [doi:10.2165/00003088-200544100-00004](https://doi.org/10.2165/00003088-200544100-00004)
