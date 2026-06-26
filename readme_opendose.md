# OpenDose Database Module Documentation

The radot opendose module maintains and queries the local cache of reference organ masses and S-value matrices from the OpenDose project. It handles lazy on-demand downloads from a central Creatis server and manages the ICRP-107 radionuclide decay database.

---

## 1. CLI Commands

### Downloading and Managing: radot opendose

Downloads and caches reference masses and S-values for a specified phantom and radionuclide.

```bash
radot opendose -m <phantom_model> -r <radionuclide> [options]
```

#### Options:
* -m, --model <AM|AF|male|female>: Shortcut for the target phantom. normalized to icrp_110_male or icrp_110_female.
* -r, --radionuclide <isotope_name>: Name of the radionuclide (e.g., Lu177, Y90, I131).
* --force: Forcefully overwrites any existing local JSON files.
* --upload: Uploads the local database contents to a remote server.
* --upload-dest <scp_path>: Target destination path for scp (e.g. user@server:public_html/radot/opendose).

---

### Database Status and Inventory: radot opendose status

Displays the locally cached database status, phantoms, radionuclides, and file sizes.

```bash
radot opendose status [options]
```

#### Options:
* -m, --model <phantom_model>: If specified, prints the list of organs and their cached reference masses.
* -r, --radionuclide <isotope_name>: If specified along with --model, prints the entire S-value matrix (Source Organ, Target Organ, S-value in mGy/MBq.s).

#### Example Inventory Output:
```text
Local OpenDose Database Inventory:
Path: /usr/local/lib/python3.14/site-packages/radot/data/opendose

├── icrp_110_female/
│   ├── [OK] masses.json (41 organs)
│   └── [OK] lu177.json.gz (41 source organs, 1681 S-values)
└── icrp_110_male/
    ├── [OK] masses.json (41 organs)
    └── [OK] lu177.json.gz (41 source organs, 1681 S-values)

Total space used: 184.2 KB (2 Radionuclides, 2 Phantoms)
```

---

## 2. Python API

You can programmatically query phantoms, half-lives, and retrieve S-value files.

### get_physical_half_life_h

Retrieves the physical half-life of any radionuclide in hours from the ICRP-107 database.

```python
from radot.db import get_physical_half_life_h

# Look up Lu177 decay half-life in hours
tp_hours = get_physical_half_life_h("Lu177")
print(f"Lu177 Half-Life: {tp_hours:.2f} hours")  # Prints: 159.36 hours
```

---

### resolve_organ_name

Resolves custom or AAPM TG-263 organ name strings to standard OpenDose names.

```python
from radot.db import resolve_organ_name

# Map Kidney_L to standard OpenDose organ name
name = resolve_organ_name("Kidney_L")
print(name)  # Prints: Left kidney

# Map a custom spelling
name_2 = resolve_organ_name("stomach_wall")
print(name_2)  # Prints: Stomach wall
```

---

### download_opendose_data

Triggers an explicit download of the reference masses and S-values.

```python
from radot.db import download_opendose_data

download_opendose_data(
    phantom="male",
    isotope="Lu-177",
    force=False
)
```

---

## 3. Storage and Caching Behavior

1. Lazy Loading: The dosimetry calculation module automatically calls the lazy-downloading mechanism on-the-fly. If an S-value file or masses file is missing locally during a calculation, radot retrieves it from the central repository (https://www.creatis.insa-lyon.fr/~dsarrut/radot/opendose/) and saves it in the local package data path.
2. In-Memory Caching: Downloaded masses and S-values are cached in-memory (using _masses_cache and _s_values_cache) so that subsequent calculations or lookups do not read from disk repeatedly.
3. Direct OpenDose Fallback: If downloading from the central Creatis registry fails, radot attempts to query and rebuild the required data by scraping the official OpenDose registry. If both sources are unavailable, a clear error is raised.

DEFAULT_DATABASE_URL = "https://www.creatis.insa-lyon.fr/~dsarrut/radot/opendose/"
