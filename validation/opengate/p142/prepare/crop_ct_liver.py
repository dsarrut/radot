#!/usr/bin/env python3
import click
import numpy as np
import SimpleITK as sitk


@click.command()
@click.option("--ct", default="ct_4mm.nii.gz", show_default=True, help="Input CT image")
@click.option(
    "--liver",
    default="liver_4mm.nii.gz",
    show_default=True,
    help="Input liver activity image",
)
@click.option(
    "--output-ct",
    default="ct_4mm_cropped.nii.gz",
    show_default=True,
    help="Output cropped CT",
)
@click.option(
    "--output-liver",
    default="liver_4mm_cropped.nii.gz",
    show_default=True,
    help="Output cropped liver",
)
@click.option(
    "--factor",
    default=2.0,
    show_default=True,
    help="Height multiplier relative to liver height",
)
def main(ct, liver, output_ct, output_liver, factor):
    """Crop CT and liver activity images along Z to a multiple of the liver height."""
    click.echo(f"Reading liver activity image: {liver}...")
    liver_img = sitk.ReadImage(liver)
    liver_arr = sitk.GetArrayViewFromImage(liver_img)

    z_indices = np.where(liver_arr > 0)[0]
    if len(z_indices) == 0:
        raise click.ClickException("Liver mask is empty (no voxels > 0).")

    z_min, z_max = int(z_indices.min()), int(z_indices.max())
    liver_height = z_max - z_min + 1
    target_height = int(round(liver_height * factor))
    z_size = liver_img.GetSize()[2]

    # Center crop around liver along Z
    z_center = (z_min + z_max) / 2.0
    z_start = int(round(z_center - target_height / 2.0))
    z_start = max(0, min(z_start, z_size - target_height))
    z_end = min(z_size, z_start + target_height)

    click.echo(f"Liver Z range: [{z_min}, {z_max}] (height: {liver_height} slices)")
    click.echo(
        f"Cropping Z range: [{z_start}, {z_end}] ({z_end - z_start} slices, factor: {factor}x)"
    )

    liver_cropped = liver_img[:, :, z_start:z_end]
    sitk.WriteImage(liver_cropped, output_liver)
    click.echo(f"Saved cropped liver image: {output_liver}")

    click.echo(f"Reading and cropping CT image: {ct}...")
    ct_img = sitk.ReadImage(ct)
    ct_cropped = ct_img[:, :, z_start:z_end]
    sitk.WriteImage(ct_cropped, output_ct)
    click.echo(f"Saved cropped CT image: {output_ct}")


if __name__ == "__main__":
    main()
