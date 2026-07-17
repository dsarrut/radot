# Design Proposal: SUL (Lean Body Mass SUV) Conversion

This document describes the CLI interface, equations, and integration plan for converting Standardized Uptake Values (SUV) to Lean Body Mass normalized values (SUL) in `radot`.

---

## 1. Mathematical Model & LBM Alternatives

SUL is computed by normalizing the activity concentration by the patient's Lean Body Mass (LBM) instead of their total Body Weight (BW):

$$
\text{SUL} = \text{SUV}_{\text{bw}} \cdot \frac{\text{LBM}}{\text{Body Weight}}
$$

`radot` supports three alternative strategies to estimate/calculate LBM:

### Strategy A: James (1976) Formula
The standard clinical formula used in PET protocols (e.g., PERCIST):

$$
\text{LBM}_{\text{male}} = 1.10 \cdot W - 128 \cdot \left(\frac{W}{H}\right)^2
$$

$$
\text{LBM}_{\text{female}} = 1.07 \cdot W - 148 \cdot \left(\frac{W}{H}\right)^2
$$

Where:
*   $W$ is the patient weight in kg.
*   $H$ is the patient height in cm.

### Strategy B: Janmahasatian (2005) Formula
A newer population standard derived from a larger patient cohort that behaves better in obese patients ($BMI > 30\text{ kg/m}^2$, where the James formula starts to under-predict LBM):

$$
\text{LBM}_{\text{male}} = \frac{9.27 \cdot 10^3 \cdot W}{6.68 \cdot 10^3 + 216 \cdot \text{BMI}}
$$

$$
\text{LBM}_{\text{female}} = \frac{9.27 \cdot 10^3 \cdot W}{8.78 \cdot 10^3 + 244 \cdot \text{BMI}}
$$

Where $\text{BMI} = \frac{W}{(H/100)^2}$ with $H$ in cm and $W$ in kg.

### Strategy C: Segmentation-Based LBM (Fat Subtraction)
If a patient CT segmentation map is available and contains fat/adipose tissue classes (e.g. from TotalSegmentator's sub-segmentation of `subcutaneous_fat` and `visceral_fat`), LBM can be determined directly by subtracting the segmented Fat Mass from the patient's total Body Weight:

$$
\text{LBM}_{\text{segmentation}} = W - M_{\text{fat}}
$$

Where Fat Mass ($M_{\text{fat}}$ in kg) is estimated from the sum of the volumes of segmented fat/adipose tissues ($V_{\text{fat}}$ in mL) multiplied by the physical density of adipose tissue ($\approx 0.92\text{ g/mL}$):

$$
M_{\text{fat}} = V_{\text{fat}} \cdot 0.92 \cdot 10^{-3}
$$

---


## 2. Metadata Schema Updates

To support SUL calculation, companion JSON sidecars will support two new metadata keys:
*   `sex`: Options are `male` or `female`.
*   `height`: Floating point height in cm.
*   `unit`: Extended to support `sul` as a first-class unit choice.

Metadata can be saved to the sidecar:
```bash
radot meta spect.nii.gz --sex male --height 175.5 --weight 75.0
```

---

## 3. CLI Interface: `convert-units`

We extend the `convert-units` command to accept `sul` as a valid source/target unit. The estimation method can be selected via `--lbm-method`. Height and sex parameters can be provided explicitly as CLI arguments or resolved automatically from the companion JSON sidecar.

### CLI Option Reference
```text
Usage: radot convert-units [OPTIONS] IMAGE_PATH

  Convert units of a radioactive image (PET/SPECT).

Options:
  -f, --from-unit [bq|bqml|suv|sul]       Source unit of the image. (Resolved from sidecar if omitted)
  -t, --to-unit [bq|bqml|suv|sul]         Target unit of the image. (Required)
  -o, --output PATH                       Output path.
  --weight FLOAT                          Patient weight in kg (required for SUV/SUL).
  --activity FLOAT                        Injected activity in MBq (required for SUV/SUL).
  --height FLOAT                          Patient height in cm (required for SUL analytical methods).
  --sex [male|female]                     Patient sex (required for SUL analytical methods).
  --lbm-method [james|janmahasatian|segmentation]
                                          Method to compute Lean Body Mass (LBM).
                                          Defaults to james if height/sex are available,
                                          or segmentation if label map contains fat classes.
  -l, --labels PATH                       Optional label map path (required for segmentation-based LBM).
  -h, --help                              Show this message and exit.
```

### Examples

#### A. Converting from SUV (Body Weight) to SUL (Lean) using James Formula
If `spect.json` already contains `weight`, `height`, and `sex`:
```bash
radot convert-units spect_suv.nii.gz --to-unit sul
```

To provide or override details explicitly:
```bash
radot convert-units spect_suv.nii.gz --to-unit sul --weight 75.0 --height 175.0 --sex male --lbm-method james
```

#### B. Direct Activity Concentration to SUL Conversion using Janmahasatian Formula
```bash
radot convert-units spect_bqml.nii.gz --to-unit sul --activity 250.0 --weight 75.0 --height 175.0 --sex male --lbm-method janmahasatian
```

#### C. Converting using Segmentation-Based LBM (Fat Subtraction)
If `labels_all.nii.gz` contains segmented fat/adipose tissue regions, LBM is calculated directly from fat volume subtraction:
```bash
radot convert-units spect_suv.nii.gz --to-unit sul --weight 75.0 --lbm-method segmentation --labels labels_all.nii.gz
```

