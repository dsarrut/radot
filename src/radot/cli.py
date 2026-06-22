import click
from pathlib import Path

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.version_option(version="0.1.0")
def main():
    """radot: Radiation Dosimetry for Therapy command-line tools."""
    pass


#### HELLO WORLD ####
@main.command()
@click.option("--name", default="World", help="Name to greet.")
def hello(name):
    """Simple hello world command to verify setup."""
    click.echo(f"Hello, {name}! Welcome to radot.")


#### SEGMENTATION ####
@main.command()
@click.argument("ct_path", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-o",
    "--output",
    "output_path",
    type=click.Path(path_type=Path),
    help="Output label image path. Defaults to labels_all.nii.gz or labels_head.nii.gz depending on mode.",
)
@click.option(
    "--mode",
    default="total",
    help="TotalSegmentator mode (e.g. total, head_glands_cavities, body, etc.).",
)
@click.option(
    "--fast", is_flag=True, help="Run TotalSegmentator in fast mode (lower resolution)."
)
def segment(ct_path, output_path, mode, fast):
    """Segment CT and group/name labels according to AAPM TG-263."""
    from radot.segmentation import segment_ct

    if output_path is None:
        if mode == "total":
            output_path = Path("labels_all.nii.gz")
        elif mode == "head_glands_cavities":
            output_path = Path("labels_head.nii.gz")
        else:
            output_path = Path(f"labels_{mode}.nii.gz")

    try:
        segment_ct(ct_path, output_path, mode, fast)
    except Exception as e:
        raise click.ClickException(str(e))


#### ROIS STATISTICS ####
@main.command()
@click.argument("labels_path", required=False, type=click.Path(exists=True, path_type=Path))
@click.option(
    "-i",
    "--image",
    "image_path",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Input image (CT/PET/SPECT) to compute stats on.",
)
@click.option(
    "-t",
    "--type",
    "image_type",
    type=click.Choice(["ct", "pet", "spect"], case_sensitive=False),
    help="Image modality type.",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    type=click.Path(path_type=Path),
    help="Optional path to output stats as JSON. Defaults to printing to console.",
)
@click.option(
    "--labels", help="Comma-separated list of label IDs to filter (e.g. 5,12)."
)
@click.option(
    "--rois", help="Comma-separated list of ROI names to filter (e.g. Kidney_L,Spleen)."
)
@click.option(
    "--resample-to",
    default="labels",
    type=click.Choice(["labels", "image"], case_sensitive=False),
    help="Grid alignment target (labels or image).",
)
@click.option(
    "--sigma",
    type=float,
    help="Optional standard deviation (in mm) of the Gaussian smoothing kernel used for peak computation.",
)
@click.option(
    "--unit",
    type=click.Choice(["hu", "suv", "bqml", "bq"], case_sensitive=False),
    help="Current unit of the image. Required to determine if activity should be preserved (e.g. bq) during resampling.",
)
def stats(labels_path, image_path, image_type, output_path, labels, rois, resample_to, sigma, unit):
    """Compute volume and intensity statistics for labels."""
    from radot.stats import compute_stats_file

    # Parse labels list
    labels_list = None
    if labels:
        try:
            labels_list = [int(x.strip()) for x in labels.split(",") if x.strip()]
        except ValueError:
            raise click.BadParameter(
                "Labels must be a comma-separated list of integers."
            )

    # Parse rois list
    rois_list = None
    if rois:
        rois_list = [x.strip() for x in rois.split(",") if x.strip()]

    try:
        compute_stats_file(
            labels_path=labels_path,
            image_path=image_path,
            image_type=image_type.lower() if image_type else None,
            output_path=output_path,
            labels=labels_list,
            rois=rois_list,
            resample_to=resample_to.lower(),
            sigma=sigma,
            unit=unit,
        )
    except Exception as e:
        raise click.ClickException(str(e))


#### RESAMPlING  ####
@main.command()
@click.argument("image_path", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-s",
    "--spacing",
    required=True,
    type=float,
    nargs=3,
    help="Target spacing in mm (e.g. -s 4 4 4).",
)
@click.option(
    "-t",
    "--type",
    "image_type",
    type=click.Choice(["ct", "labels", "pet", "spect"], case_sensitive=False),
    help="Image modality type.",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    type=click.Path(path_type=Path),
    help="Optional path to output the resampled image. Defaults to <name>_<spacing>mm.<ext>.",
)
@click.option(
    "--unit",
    type=click.Choice(["hu", "suv", "bqml", "bq"], case_sensitive=False),
    help="Units for images. Defaults to HU for CT, and SUV for PET/SPECT.",
)
@click.option(
    "--preserve-total",
    is_flag=True,
    help="Preserve total activity/value instead of concentration (scales by voxel volume ratio).",
)
@click.option(
    "--gaussian/--no-gaussian",
    default=True,
    help="Apply a Gaussian anti-aliasing filter when downsampling (continuous images only). Enabled by default.",
)
def resample(
    image_path, spacing, image_type, output_path, unit, preserve_total, gaussian
):
    """Resample an image to a target spacing."""

    image_type_lower = image_type.lower() if image_type else None

    if output_path is None:
        suffix = "".join(image_path.suffixes)
        base = image_path.name[: -len(suffix)] if suffix else image_path.name

        # Format spacing to remove trailing decimals if integers (e.g. 4.0 -> 4)
        if len(set(spacing)) == 1:
            sp_str = f"{spacing[0]:g}"
        else:
            sp_str = "x".join(f"{s:g}" for s in spacing)

        output_path = Path(f"{base}_{sp_str}mm{suffix}")

    try:
        from radot.resampling import resample_image_file

        resample_image_file(
            image_path=image_path,
            spacing=spacing,
            image_type=image_type_lower,
            output_path=output_path,
            unit=unit,
            preserve_total=preserve_total,
            gaussian=gaussian,
        )
    except Exception as e:
        raise click.ClickException(str(e))


#### CONVERT UNITS ####
@main.command(name="convert-units")
@click.argument("image_path", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-f",
    "--from-unit",
    type=click.Choice(["bq", "bqml", "suv"], case_sensitive=False),
    help="Current unit of the image.",
)
@click.option(
    "-t",
    "--to-unit",
    required=True,
    type=click.Choice(["bq", "bqml", "suv"], case_sensitive=False),
    help="Target unit of the image.",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    type=click.Path(path_type=Path),
    help="Output path. Defaults to a modified filename like <base>_<to-unit>.<ext>.",
)
@click.option(
    "--weight",
    type=float,
    help="Patient weight in kg (required for SUV conversions).",
)
@click.option(
    "--activity",
    type=float,
    help="Injected activity in MBq (required for SUV conversions).",
)
def convert_units(image_path, from_unit, to_unit, output_path, weight, activity):
    """Convert units of a radioactive image (PET/SPECT)."""
    from radot.convert import convert_image_units_file

    try:
        convert_image_units_file(
            image_path=image_path,
            from_unit=from_unit,
            to_unit=to_unit,
            output_path=output_path,
            weight=weight,
            activity=activity,
        )
    except Exception as e:
        raise click.ClickException(str(e))


#### IMAGE METADATA ####
def common_meta_options(f):
    import click
    f = click.argument("image_path", type=click.Path(exists=True, path_type=Path))(f)
    f = click.option(
        "--modality",
        type=click.Choice(["ct", "pet", "spect"], case_sensitive=False),
        help="Modality type to store (ct, pet, spect).",
    )(f)
    f = click.option(
        "--unit",
        type=click.Choice(["hu", "suv", "bqml", "bq"], case_sensitive=False),
        help="Physical unit of voxel values (hu, suv, bqml, bq).",
    )(f)
    f = click.option(
        "--weight",
        type=float,
        help="Patient weight in kg.",
    )(f)
    f = click.option(
        "--activity",
        type=float,
        help="Injected activity in MBq.",
    )(f)
    f = click.option(
        "--injection-date",
        "--injection-time",
        "injection_date",
        help="Injection date and time (YYYY-MM-DD HH:MM:SS).",
    )(f)
    f = click.option(
        "--acquisition-date",
        "--acquisition-time",
        "acquisition_date",
        help="Acquisition date and time (YYYY-MM-DD HH:MM:SS).",
    )(f)
    f = click.option(
        "--time-from-injection",
        "time_from_injection",
        help="Time elapsed from injection to acquisition (day:hours:min).",
    )(f)
    return f


@main.command(name="meta")
@common_meta_options
def meta(
    image_path,
    modality,
    unit,
    weight,
    activity,
    injection_date,
    acquisition_date,
    time_from_injection,
):
    """Read or write image metadata fields."""
    _run_meta_command(
        image_path,
        modality,
        unit,
        weight,
        activity,
        injection_date,
        acquisition_date,
        time_from_injection,
    )


@main.command(name="info")
@common_meta_options
def info(
    image_path,
    modality,
    unit,
    weight,
    activity,
    injection_date,
    acquisition_date,
    time_from_injection,
):
    """Alias for 'meta' command to read or write image metadata fields."""
    _run_meta_command(
        image_path,
        modality,
        unit,
        weight,
        activity,
        injection_date,
        acquisition_date,
        time_from_injection,
    )


def _run_meta_command(
    image_path,
    modality,
    unit,
    weight,
    activity,
    injection_date,
    acquisition_date,
    time_from_injection,
):
    from radot.meta import manage_metadata

    try:
        manage_metadata(
            image_path=image_path,
            modality=modality,
            unit=unit,
            weight=weight,
            activity=activity,
            injection_date=injection_date,
            acquisition_date=acquisition_date,
            time_from_injection=time_from_injection,
        )
    except Exception as e:
        raise click.ClickException(str(e))


#### DICOM CONVERSION ####
@main.command()
@click.argument("dicom_path", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-o",
    "--output",
    "output_path",
    type=click.Path(path_type=Path),
    help="Output image path (e.g. ct.nii.gz).",
)
@click.option(
    "--all-tags",
    "all_tags",
    is_flag=True,
    help="Print all DICOM tags present in the image header.",
)
def dicom(dicom_path, output_path, all_tags):
    """Convert DICOM directory or single file to another format (NIfTI, NRRD)."""
    from radot.dicom import convert_dicom

    try:
        convert_dicom(dicom_path, output_path, all_tags)
    except Exception as e:
        raise click.ClickException(str(e))


@main.command(name="dicom-list")
@click.argument("folder_path", type=click.Path(exists=True, path_type=Path))
def dicom_list(folder_path):
    """List all DICOM series in a folder and subfolders recursively, sorted by date."""
    from radot.dicom import list_dicom_series

    try:
        list_dicom_series(folder_path)
    except Exception as e:
        raise click.ClickException(str(e))


if __name__ == "__main__":
    main()
