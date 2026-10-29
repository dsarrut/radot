#!/usr/bin/env python3
"""
Extract the SNMMI challenge zips (downloaded by hand from Deep Blue) and convert the
images needed for T5 to NIfTI with radot:

  <DATA_DIR>/nifti/<patient>/tia.nii.gz           TIA map (Bq/mL.s), sidecar unit bqsml
  <DATA_DIR>/nifti/<patient>/spect_scan1.nii.gz   SPECT of time point 1 (Bq/mL)
  <DATA_DIR>/nifti/<patient>/ct_scan1.nii.gz      CT of time point 1
  <DATA_DIR>/nifti/<patient>/masks/<voi>.nii.gz   organ and lesion masks of time point 1

Each zip is extracted into its own new folder; existing files are not overwritten.
"""

import json
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cases import DATA_DIR, EXTRACTED_DIR, PATIENTS, ZIPS, nifti_dir  # noqa: E402

from radot.dicom import convert_dicom  # noqa: E402
from radot.log import configure_logging  # noqa: E402


def extract(zip_path: Path, out_dir: Path) -> None:
    """Extract a zip and the per-patient zips inside it (skip existing folders)."""
    if not out_dir.exists():
        print(f"Extracting {zip_path.name} -> {out_dir}")
        with zipfile.ZipFile(zip_path) as z:
            z.extractall(out_dir)
    for inner in sorted(out_dir.glob("patient_*.zip")):
        d = inner.with_suffix("")
        if not d.exists():
            with zipfile.ZipFile(inner) as z:
                z.extractall(d)


def convert(dicom_dir: Path, out: Path) -> None:
    if out.exists():
        return
    out.parent.mkdir(parents=True, exist_ok=True)
    convert_dicom(dicom_dir, out)


def set_unit(image_path: Path, unit: str, modality: str = "spect") -> None:
    """MIM exports the SPECT and TIA as PT/BQML: store the real modality and unit."""
    sidecar = image_path.parent / (image_path.name[: -len(".nii.gz")] + ".json")
    meta = json.loads(sidecar.read_text()) if sidecar.exists() else {}
    meta.update({"modality": modality, "unit": unit})
    sidecar.write_text(json.dumps(meta, indent=4))


def main() -> None:
    configure_logging(quiet=True)
    for key in ("tia", "voi", "spect_ct"):
        zip_path = DATA_DIR / f"{ZIPS[key]}.zip"
        if not zip_path.exists():
            sys.exit(f"Missing {zip_path}: download it from Deep Blue (see cases.py).")
        extract(zip_path, EXTRACTED_DIR / ZIPS[key])

    for name, patient in PATIENTS.items():
        print(f"Converting {name} (challenge patient {patient.challenge_id})")
        out = nifti_dir(name)
        src = lambda key, *parts: EXTRACTED_DIR.joinpath(ZIPS[key], name, name, *parts)  # noqa: E731
        convert(src("tia", "TIA"), out / "tia.nii.gz")
        set_unit(out / "tia.nii.gz", "bqsml")
        convert(src("spect_ct", "SPECT_Cts", "scan1", "spect"), out / "spect_scan1.nii.gz")
        set_unit(out / "spect_scan1.nii.gz", "bqml")
        convert(src("spect_ct", "SPECT_Cts", "scan1", "ct"), out / "ct_scan1.nii.gz")
        for kind, vois in (("Organs", patient.organs), ("Lesions", patient.lesions)):
            for voi in vois:
                convert(src("voi", "VOI", "Scan1", "Masks", kind, voi), out / "masks" / f"{voi}.nii.gz")


if __name__ == "__main__":
    main()
