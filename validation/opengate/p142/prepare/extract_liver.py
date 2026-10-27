#!/usr/bin/env python3
import json
from pathlib import Path
import click
import SimpleITK as sitk


@click.command()
@click.option(
    "-l",
    "--labels",
    default="labels_all_4mm.nii.gz",
    show_default=True,
    help="Input labels image",
)
@click.option(
    "-o",
    "--output",
    default="liver_4mm.nii.gz",
    show_default=True,
    help="Output activity image",
)
@click.option(
    "--label-id",
    default=None,
    type=int,
    help="Liver label ID (auto-detected if companion JSON exists)",
)
def main(labels, output, label_id):
    """Extract liver from labels image and create an activity image (0 everywhere, 1 in liver)."""
    labels_file = Path(labels)
    if label_id is None:
        json_file = labels_file.with_name(labels_file.name.split(".")[0] + ".json")
        if json_file.exists():
            with open(json_file) as f:
                mapping = json.load(f)
            label_id = next(
                (int(k) for k, v in mapping.items() if v.lower() == "liver"), 9
            )
        else:
            label_id = 9

    click.echo(f"Reading {labels}...")
    img = sitk.ReadImage(str(labels_file))

    click.echo(f"Extracting liver with label ID {label_id}...")
    activity = sitk.Cast(img == label_id, sitk.sitkFloat32)

    sitk.WriteImage(activity, output)
    click.echo(f"Saved activity image: {output}")


if __name__ == "__main__":
    main()
