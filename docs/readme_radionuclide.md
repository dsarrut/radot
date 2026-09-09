# Radionuclide Database and Nuclear Decay Information

The `radot` library integrates nuclear decay physics to support targeted radionuclide dosimetry (including Local Deposition Model / LDM and MIRD formalisms). Nuclear decay properties and radiation emission spectra are resolved using two authoritative physical sources: **ICRP Publication 107** and the **IAEA LiveChart (ENSDF)**.

---

## Physical Data Sources & Prioritization

Dosimetry calculations for radionuclides—especially alpha-emitting therapy radionuclides like $^{225}\text{Ac}$, $^{223}\text{Ra}$, or $^{212}\text{Pb}$—require two distinct types of nuclear data:
1. Complete transition radiation emission spectra (energies and yields of alphas, betas, conversion electrons, Auger electrons, photons, etc.).
2. Decay chain topology and multi-generation branching ratios to account for dose delivered by short-lived progeny in secular equilibrium.

To achieve maximum clinical accuracy and consistency, `radot` adopts a clear separation of roles and prioritization policy:

| Physical Quantity | Primary Source | Secondary / Validation Source | Rationale |
| :--- | :--- | :--- | :--- |
| **Emission Spectra & $\Delta_{\text{local}}$** | **ICRP Publication 107** (`icrp107_database`) | ENSDF / IAEA | ICRP 107 is the international standard reference for medical internal dosimetry (used by OpenDose and MIRD). Local non-penetrating energy per disintegration $\Delta_{\text{local}}$ is calculated directly from its categorized transitions. |
| **Decay Chain Topology & Branching $Y_i$** | **IAEA LiveChart** (`scripts/update_decay_chains.py` $\to$ `radot/data/decay_chains.json`) | ICRP 107 | IAEA LiveChart provides full decay network trees with multi-branch percentages for complex cascades down to stable isotopes. |
| **Parent Physical Half-Life ($T_{1/2}$)** | **ICRP Publication 107** (`radot.db.get_physical_half_life_h`) | IAEA LiveChart | Ensures exact kinetic consistency with OpenDose S-values and clinical SPECT/CT acquisition scheduling. |
| **Progeny Half-Lives** | **IAEA LiveChart** | ICRP 107 | Used in the cascade table to verify the validity of the secular equilibrium condition ($T_{1/2}^{\text{daughter}} \ll T_{1/2}^{\text{parent}}$). |

### Combined Cascade Energy Formula

For radionuclides with decay cascades under secular equilibrium ($T_{1/2}^{\text{progeny}} \ll T_{1/2}^{\text{parent}}$), the total non-penetrating energy deposited per parent decay is computed as:

$$\Delta_{\text{local, cascade}} = \sum_{i \in \text{chain}} Y_i \cdot \Delta_{\text{local, ICRP107}}(i)$$

where:
- $Y_i$ is the cumulative yield (branching weight) of progeny nuclide $i$ computed from IAEA decay branching paths ($Y_{\text{parent}} = 1.0$).
- $\Delta_{\text{local, ICRP107}}(i)$ is the non-penetrating energy constant for nuclide $i$ evaluated from ICRP 107 transitions (alpha, beta-, beta+, Auger electrons, internal conversion electrons, and recoil nuclei).
- The **Cascade Factor** is defined as:
  $$\text{Cascade Factor} = \frac{\Delta_{\text{local, cascade}}}{\Delta_{\text{local, parent}}}$$

For example, $^{225}\text{Ac}$ emits a single alpha particle ($\Delta_{\text{parent}} \approx 5.92\text{ MeV}$), but its progeny cascade emits 3 additional alpha particles plus several beta decays, yielding $\Delta_{\text{cascade}} \approx 28.70\text{ MeV}$ (a cascade factor of $4.85\times$).

---

## CLI Usage

### Listing Available Radionuclides
To display all supported isotopes present in the ICRP 107 database:
```bash
radot radionuclide list
```

### Retrieving Decay Properties and Emissions Summary
To query half-life and radiation emissions breakdown for a specific radionuclide:
```bash
radot radionuclide info Lu-177
```

Example output (single-isotope decay):
```text
Radionuclide: Lu-177
Physical Half-life: 6.647 d (6.647 days / 159.53 hours)

Emissions Summary:
Type         | Num   | Total Yield | Total Energy (MeV) | Avg Energy (MeV)
─────────────────────────────────────────────────────────────────────────
beta-        | 4     | 1.0000      | 0.1342             | 0.1342          
gamma        | 6     | 0.1793      | 0.0258             | 0.1439          
X            | 60    | 1.4505      | 0.0121             | 0.0083          
auger        | 15    | 0.9529      | 0.0019             | 0.0020          
IE           | 36    | 0.1444      | 0.0102             | 0.0706          
```

### Querying Cascading Decay Chains (e.g., Ac-225)
When a radionuclide belongs to a curated decay chain, `radot radionuclide info` displays both the parent's ICRP 107 emissions table and the IAEA LiveChart secular equilibrium cascade breakdown:

```bash
radot radionuclide info Ac-225
```

Example output:
```text
Radionuclide: Ac-225
Physical Half-life: 10.0 d (10.000 days / 240.00 hours)

Emissions Summary:
Type         | Num   | Total Yield | Total Energy (MeV) | Avg Energy (MeV)
─────────────────────────────────────────────────────────────────────────
alpha        | 44    | 1.0000      | 5.7872             | 5.7872          
gamma        | 93    | 0.0723      | 0.0104             | 0.1444          
X            | 71    | 7.6045      | 0.0067             | 0.0009          
auger        | 15    | 6.6385      | 0.0050             | 0.0007          
IE           | 511   | 0.6742      | 0.0198             | 0.0294          
alpha recoil | 44    | 1.0000      | 0.1048             | 0.1048          

Decay Chain Cascade (Secular Equilibrium - Source: IAEA LiveChart):
Progeny    | T1/2         | Yield (Y_i)  | Delta_local (MeV)  | Contrib (MeV) 
──────────────────────────────────────────────────────────────────────────
Ac225      | 9.920 d      | 1.0000       | 5.9168             | 5.9168        
Fr221      | 4.9 m        | 1.0000       | 6.4288             | 6.4288        
At217      | 32.6 ms      | 1.0000       | 7.2009             | 7.2009        
Ra221      | 28 s         | 0.0010       | 6.8605             | 0.0069        
Bi213      | 45.59 m      | 0.9999       | 0.5686             | 0.5685        
Rn217      | 0.54 ms      | 0.0010       | 7.8856             | 0.0079        
Po213      | 3.706 us     | 0.9795       | 8.5370             | 8.3622        
Tl209      | 2.162 m      | 0.0214       | 0.6875             | 0.0147        
Pb209      | 3.234 h      | 1.0009       | 0.1974             | 0.1976        
Bi209      | 2.01E19 Y    | 1.0009       | stable             | 0.0000        
──────────────────────────────────────────────────────────────────────────
Total Cascade Non-Penetrating Energy: 28.7043 MeV (Parent: 5.9168 MeV, Cascade Factor: 4.85x)
```

---

## Application to Dosimetry (LDM)

In Local Deposition Model (LDM) calculations, the decay cascade is accounted for automatically:
- **CLI**: `radot dose ldm ... --cascade` (default is `--cascade`; use `--no-cascade` to isolate parent-only dose).
- **Python**:
  ```python
  from radot.dose import compute_dose_ldm

  # Computes dose including full daughter cascade (default cascade=True)
  dose_img = compute_dose_ldm(activity_img, "Ac-225", cascade=True)
  ```

---

## Programmatic Lookup

Radionuclide properties can be queried directly via `radot.db`:

```python
from radot.db import (
    get_physical_half_life_h,
    get_local_energy_constant_mev,
    get_decay_chain,
    get_cascade_energy_constant_mev,
)

# 1. Retrieve physical half-life in hours (from ICRP 107)
half_life_hours = get_physical_half_life_h("Lu-177")
print(f"Lu-177 Half-life: {half_life_hours:.2f} hours")

# 2. Retrieve local energy constant of parent only (MeV/Bq-s)
delta_parent = get_local_energy_constant_mev("Ac-225")
print(f"Ac-225 Parent Delta_local: {delta_parent:.4f} MeV")

# 3. Retrieve decay chain structure (from IAEA LiveChart)
chain = get_decay_chain("Ac-225")
print(f"Progeny count: {len(chain['progeny'])}")

# 4. Retrieve total cascade energy constant (summed over all daughters with ICRP 107 Delta_local)
delta_cascade = get_cascade_energy_constant_mev("Ac-225", cascade=True)
print(f"Ac-225 Total Cascade Delta_local: {delta_cascade:.4f} MeV")
```

---

## Updating Decay Chains (Maintainer Script)

Decay chain trees and branching ratios are maintained in `src/radot/data/decay_chains.json`. This curated JSON file can be refreshed from the official IAEA LiveChart REST API using the maintainer utility script:

```bash
python scripts/update_decay_chains.py
```

The script:
1. Queries the IAEA LiveChart API endpoint (`https://nds.iaea.org/relnsd/v0/data?fields=decay_rads&nuclide=...`).
2. Recursively parses decay modes, branching percentages, and daughter ground states.
3. Computes cumulative secular equilibrium yields ($Y_i$) along all decay branches down to stable isotopes.
4. Updates `src/radot/data/decay_chains.json`.

---

## Emission Types Documented

The emissions summary groups nuclear transformations into standard categories defined by ICRP Publication 107:
- `alpha`: Alpha particles (local deposition)
- `beta-`: Negatrons / electrons (local deposition)
- `beta+`: Positrons (local deposition)
- `gamma`: Gamma-ray photons (penetrating)
- `X`: Characteristic X-rays (penetrating)
- `auger`: Auger electrons (local deposition)
- `IE`: Internal conversion electrons (local deposition)
- `neutron`: Neutrons
- `alpha recoil`: Recoil daughter nuclei (local deposition)
- `annihilation`: Annihilation photons (penetrating)
- `fission`: Spontaneous fission products
- `betaD`: Delayed beta decay particles

