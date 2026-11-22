# Dosimetry: Time Integration

!!! warning "Partly implemented"
    Implemented: `radot integrate` at **ROI level** on **activity images** (fit first), with the fitted
    models `f1`–`f6a` (§3.1), the error models (§3.2), model selection (§4) and the ROI table (§5); and
    `radot dose` on the ROI table (§6). Planned: `trapezoid` (§3.3), single time point models (§3.4),
    voxel level, dose-rate images (fit after).

Time integration turns a series of images acquired at several times after injection into one
time-integrated image or one value per ROI:

| Input series | Unit in | Result | Unit out |
|--------------|---------|--------|----------|
| activity | `bq` | time-integrated activity (TIA) | `bqs` (Bq·s per voxel) |
| activity concentration | `bqml` | TIA concentration | `bqsml` (Bq·s/mL) |
| dose rate | `gy/h`, `gy/s` | absorbed dose | `gy` |

$$
\tilde{A} = \int_0^{\infty} A(t)\, dt
$$

The integral is computed from a model of the time-activity curve (TAC) $A(t)$: a fitted function, a
piecewise curve through the measured points, or, with a single time point, an assumed curve shape.

`radot integrate` does only this step. The absorbed dose is then computed by `radot dose` on the
result (see [overview](6a_dose.md)): a TIA image with `svalue` or `ldm`, or directly a dose when the
series is a dose rate. The order is the user's choice:

* **fit first**: integrate the activity images, then compute the dose once on the TIA;
* **fit after**: compute a dose rate at each time point (`radot dose -m ldm --dose-rate`, or a Monte
  Carlo map), then integrate the dose rates. This uses the anatomy (masses, CT) of each time point.

With the same anatomy at every time point, the same fit level and a model linear in the data
(`trapezoid`), both orders give the same dose: the dose engines are linear in activity. They differ
through the non-linearity of the fit, changes of the anatomy between time points, and registration.

---

## 1. Inputs

* **Images**: two or more images of the same quantity (activity or dose rate), each with its sidecar.
  The time from injection of each is read from its sidecar (`injection_date` and `acquisition_date`, or
  `time_from_injection`), as everywhere in radot. All images must have the same unit and radionuclide.
* **Decay correction** of each image: sidecar key `decay_correction` (section 2). Required.
* **Fit level**:
    * **ROI** (`-l labels.nii.gz`): the ROI total (Bq, or mass-weighted mean dose rate) at each time point
      is fitted, one curve per ROI. One label image for registered time points, or one per time point
      (JSON config). Output: a table.
    * **Voxel** (no labels): one curve per voxel. The images must be registered (same grid, aligned
      anatomy): radot does not register them. Output: an image.
* **Radionuclide**: its physical decay constant $\lambda_{\text{phys}}$ is known and built into the
  models (from the sidecar `radionuclide`, or `--isotope`).

---

## 2. Decay correction

"Decay correction" means two different things, both sources of large errors:

**1. In the images.** The reconstruction may correct the activity for the physical decay to a
reference time. DICOM tag Decay Correction (0054,1102), written by `radot dicom` in the sidecar key
`decay_correction`:

| `decay_correction` | The image holds | radot uses it as |
|--------------------|-----------------|------------------|
| `start` | the activity at the start of the acquisition | is |
| `none` | the activity during the acquisition, not corrected | is (the decay during the acquisition is neglected) |
| `admin` | the activity decay-corrected to the injection time: $A(t)\, e^{\lambda_{\text{phys}} t}$ | multiplied by $e^{-\lambda_{\text{phys}} t}$, to get back $A(t)$ |

A missing or unknown `decay_correction` is an error wherever an activity at a time $t$ is used
(`radot integrate`, `radot dose` on an activity image, SUV/SUL conversion): set it once with
`radot meta IMAGE --decay-correction start|none|admin`. It is not needed for TIA, dose, CT or label
images. `radot dicom` writes it from the DICOM tag when present (PET); for SPECT the tag is usually
absent and it must be set by hand. In the SNMMI ¹⁷⁷Lu challenge, the images of one patient were corrected to the
injection time and the other's to the scan start; participants who missed it overestimated the TIA by
43 % to 204 % (Kurkowska 2025).

**2. Before the fit.** Some methods decay-correct the measured points to $t = 0$, fit the biological
curve, then multiply back by the physical decay. radot does not: the models are fitted to the
activities at the acquisition times, with $\lambda_{\text{phys}}$ inside the model (Kletting 2013). Not
decay-correcting before the fit is recommended (Ivashchenko 2024): the physical decay is a known
parameter, not a fitted one. Exception: `jackson2013` is defined on decay-corrected data and keeps it.

**Checks** (`radot check`; `radot integrate` repeats the first one):

* `meta.decay_correction`: the key is present and valid for each activity image (missing: an error in
  `radot integrate`); a series that mixes conventions is a WARN (a WARNING in `radot integrate`:
  allowed, but a common source of error).
* `activity.mtp_trend` (exists): after the conversion above, the total activity must decrease at least
  as fast as the physical decay. A slower decrease suggests images decay-corrected without saying so.
* In the fit: a terminal clearance at the bound $\lambda_k = 0$ (pure physical decay) is reported.

---

## 3. Models

Notation of Kletting et al. 2013 (NUKFIT): $t$ the time from injection, $\lambda_{\text{phys}} = \ln 2
/ T_{\text{phys}}$ known, $A_k$ the fitted amplitudes and $\lambda_k$ the fitted biological clearance
rates, $K$ the number of fitted parameters. (Sarrut et al. 2017 use the same models with $C_k$ for the
amplitudes.)

### 3.1 Fitted models (sums of exponentials)

| Model | $A(t)$ | $K$ | Shape |
|-------|--------|-----|-------|
| `f1` | $A_1 e^{-\lambda_{\text{phys}} t}$ | 1 | physical decay only |
| `f2` | $A_1 e^{-(\lambda_1+\lambda_{\text{phys}}) t}$ | 2 | mono-exponential |
| `f3` | $A_1 e^{-(\lambda_1+\lambda_{\text{phys}}) t} + A_2 e^{-\lambda_{\text{phys}} t}$ | 3 | clearance + retained fraction |
| `f4` | $A_1 e^{-(\lambda_1+\lambda_{\text{phys}}) t} + A_2 e^{-(\lambda_2+\lambda_{\text{phys}}) t}$ | 4 | bi-exponential |
| `f4a` | $A_1 \left[e^{-(\lambda_1+\lambda_{\text{phys}}) t} - e^{-(\lambda_2+\lambda_{\text{phys}}) t}\right]$ | 3 | uptake then washout, $A(0) = 0$ |
| `f5` | $A_1 e^{-(\lambda_1+\lambda_{\text{phys}}) t} + A_2 e^{-(\lambda_2+\lambda_{\text{phys}}) t} + A_3 e^{-\lambda_{\text{phys}} t}$ | 5 | bi-exponential + retained fraction |
| `f6` | $\sum_{k=1}^{3} A_k e^{-(\lambda_k+\lambda_{\text{phys}}) t}$ | 6 | tri-exponential |
| `f6a` | `f6` with $A_3 = -(A_1 + A_2)$ | 5 | tri-exponential, $A(0) = 0$ |
| `jackson2013` | tri-exponential on decay-corrected data (Jackson et al. 2013) | — | exactly 3 time points |

The TIA is analytic: $\tilde{A} = \sum_k A_k / (\lambda_k + \lambda_{\text{phys}})$, from 0 to ∞ (with
$\lambda_k = 0$ for a physical-decay term).

Kletting 2013 also defines `f4b`, `f6b` (amplitudes summing to 100 % of the injected activity) and
`f4c` (thyroid model with an input function): to add if needed.

* **Constraint**: $\lambda_k \ge 0$ (the effective decay is at least the physical decay). A fit at the
  bound is reported.
* **Number of points** (Kletting 2013, Table I), with $J = K$ for absolute weighting and $J = K + 1$
  for relative weighting (the variance scale is also estimated, §3.2): the fit needs $J \le N$ (else an
  error), the parameter standard errors $J + 1 \le N$ (else not given, and reported), the AICc
  $J + 2 \le N$ (a candidate of `--model auto` that does not satisfy it is an error).
* **Fit failure** (no convergence; a non-physical curve, negative after injection, or a TIA $\le 0$):
  for a ROI requested with `--rois` an error, else the ROI is skipped with a WARNING (as a ROI absent at
  a time point). Planned for voxels: the failed voxels counted, reported and written in a mask image.
* **Fit**: starting points from a logarithmic grid of the rates $\lambda_k \ge 0$ (half-lives from
  about 4 min to 7000 h, and 0), the amplitudes being linear for given rates (weighted least squares,
  closed form). Exchangeable terms (e.g. the two exponentials of `f4`) are kept in one order; `f4a` has
  $\lambda_2 > \lambda_1$ (uptake faster than washout). Then all the parameters are refined together
  from the 3 best starting points (`scipy.optimize.least_squares`, trust region with the bound
  $\lambda_k \ge 0$, analytic Jacobian). A rate below $10^{-6}\lambda_{\text{phys}}$ is set to 0 (at the
  bound).

### 3.2 Error model

The variance of each measured point (Kletting 2013):

$$
\vartheta_i^2 = \nu \left(a + b\, f(t_i)^{c}\right)
$$

$a$ constant part, $b$ multiplicative part, $c$ the weighting, $\nu$ a scaling factor. radot uses
data-based weighting ($f(t_i)$ replaced by the measured $y_i$, Kletting 2013 Eq. 2), so the fit is a
weighted least squares with fixed weights. `--error-model` (no default: it changes the weights of the
points and the model selection):

| `--error-model` | $\vartheta_i$ | $\nu$ | Weighting (Kletting 2013) | $J$ |
|-----------------|---------------|-------|---------------------------|-----|
| `constant` | 1 | estimated | relative | $K + 1$ |
| `proportional` | $y_i$ | estimated | relative | $K + 1$ |
| `fsd:X` | $X\, y_i$ (e.g. `fsd:0.1`: 10 % SD) | 1 | absolute | $K$ |
| `sd:X` | $X$ (MBq) | 1 | absolute | $K$ |

With $\nu$ estimated, the likelihood uses $\nu = \text{WRSS}/N$ (Kletting 2013 Eq. 3) and the
parameter covariance $\nu (J_p^T W J_p)^{-1}$ with $\nu = \text{WRSS}/(N-K)$ (Ivashchenko 2024 Eq. 8);
with absolute weighting $\nu = 1$. The TIA standard error is propagated with the covariances
(Gaussian propagation). `proportional` and `fsd:X` need values $> 0$ at every time point.

### 3.3 Piecewise model (`trapezoid`)

No fitted function between the scans: the curve goes through the measured points. Three intervals, as in
the SNMMI challenge analysis (Kurkowska 2025):

| Part | Interval | Options |
|------|----------|---------|
| head | $0 \to t_1$ | `zero` (no activity before $t_1$), `ramp` (linear from $A(0)=0$: $A_1 t_1/2$), `const` ($A_1 t_1$), `extrap` (mono-exponential through the first two points, back to $t=0$) |
| body | $t_1 \to t_N$ | trapezoids between the points |
| tail | $t_N \to \infty$ | `phys` ($A_N/\lambda_{\text{phys}}$), `fit` (mono-exponential `f2` on the last `--tail-points` points, to give) |

`--head` and `--tail` have no default: they change the result and must be given. The physical-decay
tail is conservative but often inaccurate (it overestimates the TIA when the activity clears;
Ivashchenko 2024). A fitted tail slower than the physical decay is an error: use `--tail phys`.

### 3.4 Single time point (one image)

The curve shape is assumed, not measured:

| Model | $\tilde{A}$ | Needs |
|-------|-------------|-------|
| `physical` | $A(t)\, e^{\lambda_{\text{phys}} t} / \lambda_{\text{phys}}$ | $t$ |
| `hanscheid2018` | $A(t) \cdot 2t / \ln 2$ | $t$ (assumes $T_{\text{eff}} = 2t$ when not given) |
| `madsen2018` | $A(t)\, e^{\lambda_{\text{eff}} t} / \lambda_{\text{eff}}$, $\lambda_{\text{eff}} = \ln 2 / T_{\text{eff}}$ | $t$, $T_{\text{eff}}$ per ROI |

These are the time factors used today inside `radot dose` (see [organ level](6b_dose_organ.md)). With
`integrate` they produce a TIA, so that `radot dose -m hanscheid2018` is the same as `radot integrate
--model hanscheid2018` followed by `radot dose -m svalue`. Hänscheid 2017 is not here: it is an
empirical dose formula, not a time model. Planned: single-time-point with population-based model
selection and nonlinear mixed-effects modelling (Hardiansyah 2021, 2024).

---

## 4. Model selection (`--model auto`)

With several candidate models (`--model auto --candidates f2,f3,f4a,f4`), each is fitted and compared
with the corrected Akaike information criterion (Glatting 2007, Kletting 2013). With $N$ points, $J$
estimated parameters (§3.2) and the objective function $-2\ln P$ at the minimum (Kletting 2013 Eq. 1,
5–6):

$$
-2\ln P = \sum_{i=1}^{N} \left[\ln(2\pi\vartheta_i^2) + \ln\nu + \frac{(y_i - f(t_i))^2}{\nu\,\vartheta_i^2}\right],
\qquad
\text{AICc} = -2\ln P + 2J + \frac{2J(J+1)}{N-J-1}
$$

Every candidate must satisfy $J + 2 \le N$, else it is an error (remove it from `--candidates`); radot
does not fall back to the AIC.

The Akaike weight of each model, its probability of being the best of the candidates:

$$
w_i = \frac{e^{-\Delta_i/2}}{\sum_j e^{-\Delta_j/2}}, \qquad \Delta_i = \text{AICc}_i - \min_j \text{AICc}_j
$$

Two ways to use them (`--select`, no default):

* `best`: the model with the lowest AICc.
* `average`: the TIA averaged over the models with $w_i > 0.01$, weighted by $w_i$ (renormalised), as
  NUKFIT does (Kletting 2013 Eq. 10); standard error
  $\sum_i w_i \sqrt{\text{SE}_i^2 + (\tilde{A}_i - \tilde{A})^2}$ (Eq. 11).

Only fits that pass the quality criteria enter the selection (Kletting 2013, Table I): every parameter
with a coefficient of variation below 50 % (below 25 %: precise). The others are listed in `excluded`
with the reason; none left is an error. A parameter correlation $\ge 0.8$ and a rate at the bound 0
are reported in `issues`, not excluded. With a single `--model`, the criteria are reported, never
applied.

At voxel level, the selection is done per voxel (VoMM, Sarrut 2017): each voxel gets its own model.

---

## 5. Outputs

**Image** (voxel level): the TIA image (`bqs` or `bqsml`) or the dose (`gy`), with its sidecar
(radionuclide and patient data of the first time point, unit, model, time points, and the provenance).
With `--model auto`, also a label image of the selected model per voxel and a mask of the failed voxels.

**Table** (ROI level, JSON, `-o`): the settings, then one entry per ROI (`radot.integrate.read_integrate`
reads it; an unknown format version is an error):

```json
{
  "radot_integrate": 1, "radot_version": "0.1.0", "quantity": "activity",
  "radionuclide": "Lu177", "lambda_phys_h": 0.004345,
  "model": "auto", "candidates": ["f2", "f3", "f4a"], "select": "average",
  "error_model": "proportional", "times_h": [4.0, 24.0, 96.0, 168.0, 240.0, 336.0],
  "images": ["tp1.nii.gz", "..."], "decay_correction": ["start", "..."],
  "rois": {
    "Kidney_L": {
      "tia_mbq_h": 50.8, "tia_mbq_s": 1.83e5, "tia_se_mbq_h": 2.1, "model": "average",
      "times_h": [4.0, 24.0, 96.0, 168.0, 240.0, 336.0], "activities_mbq": [3.9, 3.1, 1.2, 0.5, 0.3, 0.1],
      "volumes_ml": [151.2, "..."],
      "fits": {
        "f4a": {"params": {"A1": 4.1, "lambda1": 0.012, "lambda2": 0.31},
                "params_se": {"...": 0}, "params_cv": {"A1": 0.04, "lambda1": 0.09, "lambda2": 0.21},
                "correlation_max": 0.62, "tia_mbq_h": 51.0, "tia_se_mbq_h": 2.0, "aicc": -12.4,
                "weight": 0.81, "r2": 0.998, "nu": 0.0021, "n_points": 6, "n_estimated": 4,
                "at_bound": [], "teff_terminal_h": 52.3, "weighted_residuals": [0.4, "..."], "issues": []},
        "f2": {"...": "..."}
      },
      "excluded": {"f3": "parameter CV >= 50%: A2 71%"}
    }
  }
}
```

`activities_mbq` are the activities at the acquisition times (after undoing an `admin` decay
correction). Planned: the time-integrated activity coefficient (TIAC, $\tilde{a} = \tilde{A} / A_0$,
in h) when the injected activity is in the sidecar.

---

## 6. Command line

```bash
# ROI level, model averaging; then S-value doses on the ROI TIA table (implemented)
radot integrate tp1.nii.gz tp2.nii.gz tp3.nii.gz tp4.nii.gz tp5.nii.gz tp6.nii.gz -l labels.nii.gz \
    --model auto --candidates f2,f3,f4a --select average --error-model proportional -o tia_rois.json
radot dose tia_rois.json -m svalue -l labels.nii.gz --ct ct.nii.gz -o dose.json

# ROI level, one model, known 10 % standard deviation of each point (implemented)
radot integrate tp1.nii.gz tp2.nii.gz tp3.nii.gz -l labels.nii.gz --model f2 --error-model fsd:0.1 \
    -o tia_rois.json

# voxel level (registered images), then LDM on the TIA image (the second line works today)
radot integrate tp1.nii.gz tp2.nii.gz tp3.nii.gz tp4.nii.gz --model f2 --error-model proportional -o tia.nii.gz
radot dose tia.nii.gz -l labels.nii.gz --ct ct.nii.gz -m ldm

# piecewise: trapezoids, linear ramp before the first scan, physical decay after the last one
radot integrate tp1.nii.gz tp2.nii.gz tp3.nii.gz -l labels.nii.gz \
    --model trapezoid --head ramp --tail phys -o tia_rois.json

# fit after: dose rate per time point (works today), then integration of the dose rates
radot dose tp1.nii.gz --ct ct1.nii.gz -m ldm --dose-rate -o dr1.nii.gz     # same for tp2, tp3
radot integrate dr1.nii.gz dr2.nii.gz dr3.nii.gz -l labels.nii.gz --model auto --candidates f2,f4a --select best \
    --error-model proportional -o dose_rois.json

# one time point, Hänscheid 2018 assumption
radot integrate tp1.nii.gz -l labels.nii.gz --model hanscheid2018 -o tia_rois.json

# images decay-corrected to the injection time: say it once in the sidecar
radot meta tp1.nii.gz --decay-correction admin

# list the models and their options
radot integrate --list-models
```

**Dose on the ROI TIA table** (`radot dose tia_rois.json -m svalue -l labels.nii.gz [--ct ct.nii.gz]`):
the OpenDose self S-value of each ROI, scaled to the ROI mass, as `svalue` on a TIA image
([organ level](6b_dose_organ.md)):

* the method is `svalue` (the time integration is done); `--teff`, `--total-activity`, `--dose-rate`,
  `--cascade`, `--hu-threshold`, a dose-rate map are errors;
* ROI mass from `-l` (required) and `--ct` (CT density, voxel-wise); without a CT, the label volume with
  water density 1 g/mL (logged). `--resample like_labels` resamples the CT to the labels grid;
* radionuclide from the table (`--isotope` that contradicts it is an error); phantom from `--phantom`,
  else the patient sex recorded in the table (from the image sidecars), else an error;
* ROIs: all those of the table, or `--rois` (a requested ROI missing from the table or the labels, or
  without OpenDose data, is an error; without `--rois`, skipped with a WARNING);
* results: as any `radot dose` ROI table, plus `mean_dose_fit_se_gy` (the dose uncertainty from the
  TIA fit only: S-value and mass not included), `tia_se_mbq_h` and `tia_model`.

Time points with one label image each: repeat `-l`, one per image, in the same order as the images
(they are then sorted by time from injection). `--resample like_labels` resamples each image to the grid
of its labels; without it a grid mismatch is an error. `--rois` selects ROIs: a requested ROI absent at
a time point, or whose fit fails, is an error; without `--rois` it is skipped with a WARNING.

---

## 7. Choosing a model

* More time points allow more parameters. With relative weighting ($\nu$ estimated, $J = K + 1$),
  the AICc of `f2` needs 5 points and `f3`, `f4a` 6; with absolute weighting ($J = K$), `f2` needs 4
  and `f3`, `f4a` 5. With 3 points, a single `f2` is fitted (no AICc; its standard error needs absolute
  weighting). Kurkowska 2025: do not fit a 4-parameter model to 4
  points.
* Organs with an almost instantaneous uptake are often well described by `f2`; lesions often show
  uptake then washout (`f4a`): a mono-exponential on a lesion can be strongly biased (Kurkowska 2025).
* The tail matters most when the last scan is early: a last measurement after several effective
  half-lives reduces the error (Ivashchenko 2024).
* In the SNMMI ¹⁷⁷Lu challenge (Kurkowska 2025, same data as the [SNMMI validation](8_validation_snmmi.md)),
  the variability of the absorbed doses due to the fitting was small for the kidneys (below 1 % for one
  patient, below 10 % for the other) and up to 25 % for the lesions. The large outliers came from errors
  (decay correction, lesions inside the liver ROI), not from the choice of the fit.
* Always look at the fits: the ROI table gives the parameters, their coefficients of variation and the
  measured values.

---

## References

<!-- Draft: Sarrut 2017 pages and Gear 2018 from memory (no full text in docs/bib): check, then add all
to dosimetry.bib. -->

* Kletting P, Schimmel S, Kestler HA, Hänscheid H, Luster M, Fernández M, et al. Molecular radiotherapy:
  the NUKFIT software for calculating the time-integrated activity coefficient. *Med Phys.*
  2013;40(10):102504. doi:10.1118/1.4820367
* Glatting G, Kletting P, Reske SN, Hohl K, Ring C. Choosing the optimal fit function: comparison of the
  Akaike information criterion and the F-test. *Med Phys.* 2007;34(11):4285–4292.
* Sarrut D, Halty A, Badel JN, Ferrer L, Bardiès M. Voxel-based multimodel fitting method for modeling
  time activity curves in SPECT images. *Med Phys.* 2017;44:6280–6288.
* Ivashchenko OV, O'Doherty J, Hardiansyah D, Cremonesi M, Tran-Gia J, Hippeläinen E, et al.
  Time-activity data fitting in molecular radiotherapy: methodology and pitfalls. *Phys Med.*
  2024;117:103192.
* Hardiansyah D, Riana A, Kletting P, Zaid NRR, Eiber M, Pawiro SA, et al. A population-based method to
  determine the time-integrated activity in molecular radiotherapy. *EJNMMI Phys.* 2021;8:82.
  doi:10.1186/s40658-021-00427-x
* Hardiansyah D, Yousefzadeh-Nowshahr E, Kind F, Beer AJ, Ruf J, Glatting G, Mix M. Single-time-point
  renal dosimetry using nonlinear mixed-effects modeling and population-based model selection in
  [177Lu]Lu-PSMA-617 therapy. *J Nucl Med.* 2024;65:566–572.
* Kurkowska S, Brosch-Lenz J, Dewaraja YK, Frey E, Sunderland J, Uribe C. An international study of
  factors affecting variability of dosimetry calculations, part 4: impact of fitting functions in
  estimated absorbed doses. *J Nucl Med.* 2025;66:441–448.
* Sjögreen Gleisner K, Chouin N, Minguez Gabina P, Cicone F, Gnesin S, Stokke C, et al. EANM dosimetry
  committee recommendations for dosimetry of 177Lu-labelled somatostatin-receptor- and PSMA-targeting
  ligands. *Eur J Nucl Med Mol Imaging.* 2022;49:1778–1809. doi:10.1007/s00259-022-05727-7
* Gear JI, Cox MG, Gustafsson J, Sjögreen Gleisner K, Murray I, Glatting G, et al. EANM practical
  guidance on uncertainty analysis for molecular radiotherapy absorbed dose calculations. *Eur J Nucl
  Med Mol Imaging.* 2018;45:2456–2474.
* Jackson PA, Beauregard JM, Hofman MS, Kron T, Hogg A, Hicks RJ. An automated voxelized dosimetry tool
  for radionuclide therapy based on serial quantitative SPECT/CT imaging. *Med Phys.* 2013;40:112503.
