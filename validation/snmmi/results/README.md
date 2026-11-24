# radot results on the SNMMI challenge data

`t5_radot.json` (task 5, written by `t5.py`) and `t4_radot.json` (task 4, written by `t4.py`): radot
values per patient and VOI (activities, TIA, masses, doses, differences with the participants). No image,
mask or voxel-level output.

They are committed so that the page `docs/8_validation_snmmi.md` is generated and the tests run without
the data; with the data, the tests check that a rerun reproduces them.

Source: derived from the SNMMI 177Lu Dosimetry Challenge data, University of Michigan Deep Blue
collection [hm50ts030](https://deepblue.lib.umich.edu/data/collections/hm50ts030), license
[CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/) (research only, no commercial use).
Participants' statistics: Brosch-Lenz et al., J Nucl Med 2023;64:1109 (`../challenge_stats.json`).

> We acknowledge Yuni Dewaraja at the University of Michigan for providing access to Lu-177 patient
> imaging data and segmentations obtained with support from R01EB022075 awarded by NIBIB and R01CA240706
> awarded by NCI, NIH. Data were shared via the University of Michigan Deep Blue Data sharing repository.
