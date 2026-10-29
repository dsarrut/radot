"""
SNMMI 177Lu Dosimetry Challenge (2021): patients, VOIs and data layout.

Single source of truth for the SNMMI evaluation (see dev_docs/dev_readme_082_snmmi_challenge.md).
The data (University of Michigan Deep Blue, CC BY-NC 4.0) are not redistributed: download the
zips by hand into DATA_DIR (default <repo>/data/snmmi_challenge, or $RADOT_SNMMI_DATA), then run
prepare_data.py.

Acknowledgement requested by the data owners for any publication or presentation:
"We acknowledge Yuni Dewaraja at the University of Michigan for providing access to Lu-177
patient imaging data and segmentations obtained with support from R01EB022075 awarded by
NIBIB and R01CA240706 awarded by NCI, NIH. Data were shared via the University of Michigan
Deep Blue Data sharing repository."
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List

REPO_DIR = Path(__file__).resolve().parents[2]
DATA_DIR = Path(os.environ.get("RADOT_SNMMI_DATA", REPO_DIR / "data" / "snmmi_challenge"))
EXTRACTED_DIR = DATA_DIR / "extracted"
NIFTI_DIR = DATA_DIR / "nifti"
RESULTS_DIR = Path(__file__).resolve().parent / "results"

# Deep Blue works (zip name = work id)
ZIPS = {
    "spect_ct": "dn39x173n",  # Multi-time point SPECT/CT, doi:10.7302/0n8e-rz46
    "voi": "4t64gn408",  # Lesion and organ VOIs, doi:10.7302/vhrh-qg23
    "tia": "gb19f607x",  # SPECT TIA maps (Bq/mL.s, registered to scan 1), doi:10.7302/9nct-bk44
    "planar": "1j92g763g",  # Planar whole body (not used)
    "projections": "xp68kg44x",  # SPECT projections + mu maps (not used)
    "diagnostic": "js956g027",  # Pre-therapy diagnostic images (not used)
}

# VOI mask folder name -> radot ROI name (TG-263); lesions keep their name
ORGAN_NAMES = {
    "liver": "Liver",
    "ltkidney": "Kidney_L",
    "rtkidney": "Kidney_R",
    "spleen": "Spleen",
}


@dataclass
class Patient:
    deep_blue_id: str  # folder name in the Deep Blue data
    challenge_id: str  # name in the challenge papers
    sex: str  # OpenDose phantom
    injected_mbq: float
    days: List[int]  # imaging days after injection
    organs: List[str]
    lesions: List[str] = field(default_factory=list)


# Patient 4 = A and patient 6 = B, from the injected activities (DICOM header), the
# imaging days and the spleen (Brosch-Lenz 2023: B is splenectomised). The VOI readme
# sentence "Patient 4 does not have a spleen" is wrong: patient 4 has a spleen mask.
PATIENTS: Dict[str, Patient] = {
    "patient_4": Patient(
        "patient_4", "A", "male", 7210.4, [0, 1, 4, 5],
        ["liver", "ltkidney", "rtkidney", "spleen"], ["lesion1", "lesion2"],
    ),
    "patient_6": Patient(
        "patient_6", "B", "female", 7305.0, [0, 1, 4, 8],
        ["liver", "ltkidney", "rtkidney"], ["lesion1", "lesion2", "lesion3", "lesion4"],
    ),
}


# radot VOI name (t5.py results) -> VOI name in the challenge statistics
CHALLENGE_VOI = {
    "Liver": "Liver",
    "Liver (no lesions)": "Liver",
    "Spleen": "Spleen",
    "Kidney_R": "R Kidney",
    "Kidney_L": "L Kidney",
    "Kidneys (total)": "Total Kidney",
    "Lesion1": "Lesion 1",
    "Lesion2": "Lesion 2",
    "Lesion3": "Lesion 3",
    "Lesion4": "Lesion 4",
}

CHALLENGE_STATS_FILE = Path(__file__).resolve().parent / "challenge_stats.json"


def challenge_stats(task: str, patient: str) -> Dict[str, dict]:
    """Participants' statistics (Gy) for a task ('T1'..'T5') and a Deep Blue patient name."""
    import json

    stats = json.loads(CHALLENGE_STATS_FILE.read_text())
    return stats["tasks"][task][PATIENTS[patient].challenge_id]


def nifti_dir(patient: str) -> Path:
    return NIFTI_DIR / patient


def has_nifti(patient: str) -> bool:
    d = nifti_dir(patient)
    return (d / "tia.nii.gz").exists() and (d / "ct_scan1.nii.gz").exists()
