# OpenDose Database

radot uses the reference organ masses and S-values of the [OpenDose](https://opendose.org) collaboration for its organ-level dosimetry (S-value methods: `madsen2018`, `svalue`, and the cross-doses). The ICRP 110 adult male (AM) and female (AF) phantoms are available.

## Data, licence and disclaimer

- **Source**: OpenDose collaboration. Please cite: Chauvin M. et al., "OpenDose: open-access resource for nuclear medicine dosimetry", J Nucl Med 2020, [doi:10.2967/jnumed.119.240366](https://doi.org/10.2967/jnumed.119.240366).
- **Licence**: the OpenDose data are under [Creative Commons Attribution-ShareAlike 4.0](https://creativecommons.org/licenses/by-sa/4.0/) (CC BY-SA 4.0). The copy used by radot is published on Zenodo under the same licence: [doi:10.5281/zenodo.23841128](https://doi.org/10.5281/zenodo.23841128). The radot code is under LGPL-3.0.
- **Disclaimer (OpenDose)**: the data are intended for the nuclear medicine research community, for educational and research purposes only. They are not intended to be used in the clinical management of patients.

## Where the data come from

1. **Shipped with radot**: the organ masses of both phantoms and the Lu-177 S-values.
2. **Downloaded on first use**: the other radionuclides. The whole archive (12 MB) is downloaded once from Zenodo, its SHA-256 is checked, and it is extracted into the cache folder.
   The archive contains 36 radionuclides: Ac-225, Ag-111, At-211, Bi-212, Bi-213, Cu-64, Cu-67, Er-169, F-18, Ga-68, Ho-166, I-123, I-124, I-125, I-131, In-111, Lu-177, P-32, Pb-203, Pb-212, Ra-223, Re-186, Re-188, Sc-44, Sc-47, Sm-153, Sr-89, Tb-149, Tb-152, Tb-155, Tb-161, Tc-99m, Th-227, Y-86, Y-90, Zr-89. The S-values are those of the radionuclide alone, without its progeny.
3. **OpenDose website** (fallback): if the archive cannot be downloaded, or for a radionuclide it does not contain, radot reads the S-values from the OpenDose website (one radionuclide at a time, slower).

The cache folder is `$RADOT_DATA_DIR/opendose/` if the environment variable is set, else the user cache: `~/.cache/radot/opendose/` (Linux), `~/Library/Caches/radot/opendose/` (macOS), `%LOCALAPPDATA%\radot\opendose\` (Windows). radot never writes into its installation folder. Without network access, download the archive by hand from Zenodo and unzip the content of its top folder into the cache folder.

S-values are in mGy/(MBq·s), masses in g.

## CLI commands

### Download the archive: `radot opendose download`

```bash
radot opendose download            # all phantoms and radionuclides of the archive, into the cache
radot opendose download --force    # overwrite the files already in the cache
```

### Download one radionuclide from the OpenDose website: `radot opendose -m -r`

For a radionuclide not in the archive:

```bash
radot opendose -m <AM|AF|male|female> -r <radionuclide> [--force]
```

### Inventory: `radot opendose status` (or `list`)

```bash
radot opendose status                  # files shipped with radot, then the cache
radot opendose status -m AM            # organs and reference masses of the phantom
radot opendose status -m AM -r lu177   # S-value matrix (source organ, target organ, mGy/MBq.s)
```

## Python API

```python
from radot.db import get_s_value, get_reference_mass, resolve_organ_name, get_physical_half_life_h

get_s_value("male", "lu177", "Liver", "Liver")     # mGy/(MBq.s), self-dose S-value
get_reference_mass("female", "Kidney_L")           # g
resolve_organ_name("Kidney_L")                     # 'Left kidney' (TG-263 -> OpenDose name)
get_physical_half_life_h("Lu177")                  # h, from ICRP 107
```

Other functions: `download_opendose_archive(force=False)` (the Zenodo archive), `download_opendose_data(phantom, isotope, force=False)` (one radionuclide from the OpenDose website), `data_dir()` and `opendose_cache_dir()` (cache location).

Masses and S-values are loaded on demand and kept in memory for the following lookups.
