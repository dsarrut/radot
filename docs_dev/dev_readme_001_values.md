# 001 — The values of radot

Why radot is worth developing and using now that AI coding agents can write dosimetry code in minutes. This note sets the direction; the concrete rules are in `dev_readme_002_principles.md`.

---

## 1. What changed

An agent can write an LDM, a Madsen 2018 formula, a DICOM converter or a Bq→SUV conversion in a few minutes. **Code is no longer the scarce part.** A toolkit whose only value is "you don't have to write it" is worth much less than before.

What stays scarce is **trust**: knowing that a dose value is right, that it means the same thing in another centre, and being able to explain how it was obtained. An agent produces code, not trust.

So radot's value is not its code. It is everything that makes a dose value trustworthy, comparable and explainable.

---

## 2. The values

### 2.1 Validation evidence (the main asset)
- A script written by an agent for one study is not validated. radot is compared with OpenGATE Monte Carlo (p142: Lu177, Y90, I131, Ac225) and with the SNMMI 177Lu Dosimetry Challenge participants, and `pytest validation` keeps these comparisons as regression tests.
- A statement such as "liver self-dose within 0.4 % of OpenGATE Monte Carlo for Lu177 with LDM; Ac225 with the decay cascade within 2.7 %" is worth more than the code that produces it.
- Consequence: **radot is worth as much as its validation covers.** Where a method is not validated, radot is barely better than an agent's script.

### 2.2 Domain decisions, written down and reviewed
- Total activity kept when resampling Bq images; nearest neighbour for labels; anti-aliasing only for continuous images.
- No automatic "Skeleton" grouping; stop on an explicitly requested unmapped ROI.
- ICRP 107 for emissions, IAEA LiveChart for decay chains; organ-mass convention; mass scaling of S-values.
- These are the mistakes an agent makes without noticing. The code is cheap to rewrite; the reviewed decisions are not.

### 2.3 Curated, licensed reference data
- Decay chains (`decay_chains.json`), OpenDose masses and S-values (Zenodo hosting, checksums, CC BY-SA), TotalSegmentator → TG-263 → OpenDose mappings.
- Slow and error-prone to build, and must be checked by a person. Asked for "the S-value liver ← kidney", an agent often invents a number.

### 2.4 Input quality control
- Wrong doses come more often from the data than from the models (units, GBq vs MBq, injection time, weight, misregistration).
- `radot check` encodes expert knowledge that an agent does not think to apply.

### 2.5 Comparable and reproducible results between centres
- The SNMMI challenge showed a large spread between participants. If every group has an agent write its own pipeline, the spread grows.
- A shared, versioned, citable tool ("radot x.y, madsen2018, OpenDose Zenodo record 23841128 v1.0") makes results comparable. This matters **more** now that code is cheap, not less.

### 2.6 Provenance (to build)
- Every output records how it was obtained: versions, data checksums, parameters and where each came from, fallbacks and decisions, QC findings.
- Turns "trust radot" into "check what radot did". One-off agent scripts almost never do this. See `dev_readme_002_principles.md`, section 6.

### 2.7 A trusted core for agents
- The likely main users of radot are agents writing patient- or study-specific pipelines.
- The model: **radot provides a small, strict, validated core; agents write the glue.** The glue no longer needs to be in radot.
- radot must therefore be easy for agents to use: structured (JSON) outputs, errors that state the fix, stable schemas, a user-facing agent guide, possibly an MCP server.

---

## 3. What has lost value

- Convenience wrappers and CLI comfort as an end in themselves.
- Breadth of features: adding methods without validation.
- Study-specific pipelines, batch runners, reporting: agents now write these cheaply.

---

## 4. Consequences for development

- **Invest in:** validation (lesions, multi-time-point, more SNMMI tasks), `radot check`, curated data, provenance, schemas (sidecar, results JSON), agent interfaces.
- **Invest less in:** new unvalidated methods, generic image processing, workflow features.
- **Refuse:** anything that belongs to the glue (see `dev_readme_002_principles.md`, section 1).

### Agent interfaces, in order
1. A **user-facing** agent guide (skill / `SKILL.md`): typical pipeline, which method for which data, traps, always run `radot check` first. (`AGENTS.md` is for agents *developing* radot, not *using* it.)
2. JSON output everywhere (`--json`) and provenance.
3. A thin MCP server (`radot[mcp]`, local stdio, file paths not images): first the read-only lookups (radionuclide data, decay chains, S-values, masses, ROI name resolution), then the processing tools (check, meta, dose, stats). Generated from the Python API, no duplicated logic; the stop rules stay hard errors.

The lookups are useful even to agents writing their own code: they ground the physics constants instead of letting the agent invent them.

---

## 5. Publication angle

A research paper on radot (EJNMMI Physics, original research) where this reasoning explains the design, backed by experiments:
1. **Agents with and without radot:** several models × several runs compute organ doses from the SNMMI TIA maps from scratch, then with radot. Measure the spread and the typical errors.
2. **Error injection:** realistic input errors (GBq/MBq, swapped dates, wrong radionuclide, missing decay correction, shifted labels); fraction caught by `radot check`, and dose impact if missed.
3. **Provenance case study:** the SNMMI differences (organ-mass convention, mapping stage) read directly from the provenance.
4. Existing validation (OpenGATE, SNMMI).

Prerequisites: tagged release + Zenodo DOI, provenance implemented, AI use in radot development disclosed (consistent with the thesis: validation, not authorship, is the guarantee). A short JOSS paper can run in parallel.

