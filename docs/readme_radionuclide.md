# Radionuclide Database and Nuclear Decay Information

The radot library integrates nuclear decay metadata sourced from the ICRP Publication 107 reference database. This module supports dosimetry calculations by resolving physical decay properties and emission spectra.

## CLI Usage

### Listing Available Radionuclides
To display all supported isotopes present in the database:
```bash
radot radionuclide list
```

### Retrieving Decay properties and Emissions Summary
To query half-life and radiation emissions breakdown for a specific radionuclide:
```bash
radot radionuclide info Lu-177
```

Example output:
```text
Radionuclide: Lu-177
Physical Half-life: 6.647 d (6.647 days / 159.53 hours)

Emissions Summary:
Type         | Num   | Total Yield | Total Energy (MeV) | Avg Energy (MeV)
─────────────┼───────┼─────────────┼────────────────────┼─────────────────
beta-        | 4     | 1.0000      | 0.1342             | 0.1342          
gamma        | 6     | 0.1793      | 0.0258             | 0.1439          
X            | 60    | 1.4505      | 0.0121             | 0.0083          
auger        | 15    | 0.9529      | 0.0019             | 0.0020          
IE           | 36    | 0.1444      | 0.0102             | 0.0706          
```

## Programmatic Lookup

Radionuclide physical half-life lookup is exposed via the database layer:

```python
from radot.db import get_physical_half_life_h

# Retrieve physical half-life in hours
half_life_hours = get_physical_half_life_h("Lu-177")
print(f"Lu-177 Half-life: {half_life_hours} hours")
```

## Emission Types Documented

The emissions summary groups nuclear transformations into standard categories:
- alpha: Alpha particles
- beta-: Negatrons (electrons)
- beta+: Positrons
- gamma: Gamma-ray photons
- X: Characteristic X-rays
- auger: Auger electrons
- IE: Internal conversion electrons
- neutron: Neutrons
- alpha recoil: Recoil nuclei
- annihilation: Annihilation radiation
- fission: Spontaneous fission products
- betaD: Delayed beta decay particles
