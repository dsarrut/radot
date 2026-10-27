#!/usr/bin/env python3
"""Download the p142 images and OpenGATE outputs from Zenodo into p142/data/."""

import sys

from benchmarks import DATA_DIR, DATA_URL, download_data, has_data

if __name__ == "__main__":
    force = "--force" in sys.argv
    if has_data() and not force:
        print(f"Data already in {DATA_DIR} (use --force to download again)")
    else:
        print(f"Downloading {DATA_URL} ...")
        download_data(force=force)
        print(f"Data in {DATA_DIR}")
