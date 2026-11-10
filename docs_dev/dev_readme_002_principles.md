# 002 — Core principles of radot

The rules that follow from `dev_readme_001_values.md`: what radot is, what it is not, and how its outputs must be trusted. Status per item: **[done]** in the code today, **[partial]**, **[todo]**.

---

## 1. Scope: a small trusted core

### 1.1 The inclusion criterion
> A feature belongs in radot if a silent error in it changes the dose, and an independent rewrite (by an agent or a person) would plausibly do it wrong or differently.

Test: *"if every lab had an agent rewrite this, would their doses differ without anyone noticing?"* Yes → core, with validation or invariant tests. No (or the error is loud: crash, visibly wrong image) → glue, outside radot.

Second bound: **the core does not grow faster than what can be validated and maintained.** Each public function is either compared with a reference (Monte Carlo, challenge, published values) or covered by an invariant test (activity conserved, units consistent).

### 1.2 Three layers
1. **Core**: physics data, dose models and time integration, unit conversions and resampling invariants, nomenclature, QC, provenance, schemas (sidecar, results). Small, strict, validated, versioned, citable.
2. **Adapters** (optional extras): DICOM pixel I/O, TotalSegmentator, OpenGATE dose-rate maps, MCP server. Thin, replaceable without touching the core.
3. **Glue** (outside radot): study pipelines, orchestration, batch, reporting. Written by users or agents on top of the core.

### 1.3 In the core
| Area | Note |
|---|---|
| Dose models, time integration, strategies, MTP fit, uncertainty (`dose/`) | the dose itself |
| Reference data (`db`): decay chains, OpenDose masses and S-values | curated, licensed, checksummed |
| Nomenclature TS → TG-263 → OpenDose | wrong mapping = wrong S-value, silently |
| Units, SUV/SUL, LBM (`convert`) | formula and unit traps |
| Activity-preserving resampling | the invariants are the value, SimpleITK does the rest |
| Mass from CT density, SUVpeak, ROI stats; DVH **[todo]** (mass-weighted; today only in vvv, per voxel) | definitions must be identical everywhere |
| DICOM **metadata** extraction | injection time and activity, decay-correction flags; vendor SPECT scaling **[todo]** |
| Sidecar and results schemas | the contract between tools |
| `radot check` | QC |
| Provenance | see section 6 |
| Validation suite (repo, not package) | justifies the rest |

### 1.4 Not in radot
- Generic image processing (registration, filtering, cropping): elastix, ANTs, SimpleITK. radot only *checks* the result.
- Segmentation models, lesion segmentation: adapters and label contract only.
- Monte Carlo: OpenGATE; radot reads the dose-rate map, never runs GATE.
- Viewers, GUI: vvv.
- Pipelines specific to a study or centre, workflow engines, batch runners, PACS/database connections, PDF reports.
- Methods from the literature that cannot be validated or are not needed.
- Automatic repair of inputs (guessing a missing date, fuzzy ROI matching during a run).
- Heavy runtime dependencies (only optional extras).

---

## 2. The sidecar is the contract **[done, schema todo]**

- Every image has a companion JSON sidecar (`image.nii.gz` → `image.json`) holding the therapy metadata: modality, unit, weight, height, sex, injected activity (MBq), injection and acquisition dates, radionuclide. Image headers are not used for these.
- **Read the sidecar first.** Any parameter needed by a tool is looked up in the sidecar; a CLI option or API argument overrides it explicitly.
- **No static defaults** in the CLI for these parameters (`None`, resolved dynamically), and not `required=True`.
- **Explicit failure**: missing in both option and sidecar → clear `ValueError` naming the attribute.
- **Inheritance**: every tool that creates an image (resample, convert-units, dose) writes a matching sidecar, updated (unit, modality) and linked to its parents (section 6).
- **[todo]** A versioned, formal schema (`schema_version` field, JSON Schema file shipped in `data/`), so that sidecars stay readable across radot versions and other tools (vvv, agents) can validate them. Consider alignment with BIDS (`GeneratedBy`, sidecar conventions).

---

## 3. Explicit over silent **[done, to keep]**

- **Fail loudly** on anything that would make the dose wrong: unknown radionuclide, missing unit or dates, unmapped explicitly requested ROI, unmapped bone (no "Skeleton" grouping).
- **Warn and skip** only in generic loops (all labels of an image), and record the skip (provenance).
- **No fuzzy matching during a run.** Fuzzy suggestions only in interactive help (`match`, `verify`).
- **Errors state the fix**: which attribute, which command (`radot meta ... --activity`), which file.
- **Never guess an input.** radot reports; the user (or agent) fixes explicitly.

---

## 4. Physics and data invariants **[done]**

- Total activity preserved when resampling absolute-activity images (scale by voxel-volume ratio).
- Linear interpolation for continuous images, nearest neighbour for labels; anti-aliasing on downsampling of continuous images only.
- ICRP 107: emissions and $\Delta_{\text{local}}$. IAEA LiveChart: decay-chain topology and cumulative yields.
- AAPM TG-263 is the naming hub; 3-stage cascade (TG-263 route, direct TS→OD bypass, raw retain).
- Reference data in flat JSON in `src/radot/data/`, never hardcoded. Downloaded data checksummed, cached outside the package (`$RADOT_DATA_DIR` or `~/.cache/radot`), never written into the install.
- Licences respected and stated: OpenDose CC BY-SA (citation), SNMMI data CC BY-NC (never committed or redistributed).

---

## 5. Validation as the guarantee **[done, to extend]**

- Any change to a dose model, S-values, masses or decay data keeps `pytest validation` passing, and the generated pages (`docs/7_*`, `docs/8_*`) are regenerated.
- Single source of truth for benchmarks and tolerances (`benchmarks.py`, `cases.py`), never duplicated in `tests/`.
- A new method enters the core with a validation or an invariant test, or does not enter.
- Validation coverage is published with the method (docs state what is validated, against what, within what).
- **[todo]** Extend: lesions (note 085), multi-time-point time integration, more SNMMI tasks.

---

## 6. Provenance **[partial]**

Every result carries the record of how it was produced, so any dose can be traced, reproduced and explained without asking who ran it.

Today: `convert-units` writes the output `sha256` and the decay-correction keys (`decay_correction_factor`, `decay_correction_reference`) in the sidecar; downloaded OpenDose data are checksummed. No `provenance` block yet; these keys move into it.

### 6.1 Content
A `provenance` block in every output sidecar and results JSON:
- `generated_by`: radot version, git commit, command line or API call, timestamp, Python and key dependency versions.
- `data`: versions and SHA-256 of the reference data used (OpenDose record/version/phantom/file, decay chains, nomenclature tables).
- `inputs`: file basename, SHA-256, unit for each input image.
- `parameters`: each value with its **source** (`sidecar`, `cli_option`, `method_params`, `config`, `default`).
- `decisions`: every choice made on the user's behalf: resampling (grid, interpolation), ROI mapping stage, skipped ROIs and reason, mass source (CT density or phantom), mass scaling, cascade, fallbacks.
- `qc`: `radot check` findings at compute time.
- `parents`: sidecars of the inputs, giving the lineage DICOM → convert → resample → dose.

### 6.2 Rules
- Collected by a `Provenance` object created per top-level call; lower code adds `param(...)` and `decision(...)`.
- Decisions recorded from the same code paths that log warnings, so log and record never diverge.
- **Privacy**: basenames or hashes only by default (paths often contain patient names); never copy identifying DICOM fields.
- Hashing optional if it becomes slow (`--no-hash`).
- Later: `radot replay result.json` reruns from the provenance and compares.

---

## 7. Quality control **[done, to extend]**

- `radot check` runs before dosimetry; it only reads and reports, never modifies, never raises on a suspicious value.
- Levels FAIL / WARN / INFO / OK; thresholds in `data/check_thresholds.json`, user-overridable; JSON report; exit code for scripts and agents.
- New checks go through the registry (`register_check`, code `category.name`, declared inputs).
- **[todo]** Run automatically (or recorded) at dose time; image-based checks (misregistration SPECT/CT, labels outside the body, activity outside the ROIs).

---

## 8. Agent-legible by design **[partial]**

- Structured outputs: every command can write JSON **[partial]**: today only through `-o *.json` (`dose`, `stats`, `check`); `meta`, `radionuclide`, `opendose`, `roi` print text only; no `--json` option yet.
- Stable names and schemas (sidecar, results, check report), versioned.
- Errors that state the fix (section 3).
- A **user-facing** agent guide (skill / `SKILL.md`): pipeline, method choice, traps, `radot check` first **[todo]**. `AGENTS.md` stays the developer guide.
- Optional MCP server (`radot[mcp]`): thin, generated from the Python API, file paths not images, stop rules kept as hard errors **[todo]**.

---

## 9. Engineering

- PyPA src-layout; public dose API re-exported by `radot/dose/__init__.py`.
- Light core dependencies (SimpleITK, pydicom, click, numpy, tqdm, icrp107-database); everything heavy is an extra.
- `logging`, no `print()` in the library.
- Tests run without patient data and without TotalSegmentator.
- Releases tagged and archived on Zenodo (DOI per version), so that a result can cite the exact radot version.
- AI-assisted development is allowed and disclosed; the guarantee is validation and tests, not authorship.
- Research and education use only, not clinical management of patients (as OpenDose).
