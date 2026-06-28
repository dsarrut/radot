import click
from pathlib import Path

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])


@click.group(context_settings=CONTEXT_SETTINGS)
@click.option("-v", "--verbose", is_flag=True, help="Show verbose debugging logs and Python stack trace on errors.")
@click.option("-q", "--quiet", is_flag=True, help="Silence progress/status output, showing only errors.")
@click.version_option(version="0.1.0")
@click.pass_context
def main(ctx, verbose, quiet):
    """radot: Radiation Dosimetry for Therapy command-line tools."""
    ctx.ensure_object(dict)
    ctx.obj["verbose"] = verbose
    ctx.obj["quiet"] = quiet
    from radot.log import configure_logging
    configure_logging(verbose=verbose, quiet=quiet)


def handle_cli_exception(e: Exception) -> None:
    """Helper to handle exceptions in CLI commands, showing tracebacks in verbose mode."""
    ctx = click.get_current_context(silent=True)
    verbose = False
    if ctx and ctx.obj:
        verbose = ctx.obj.get("verbose", False)
    if ctx and ctx.params and ctx.params.get("verbose"):
        verbose = True
    if verbose:
        raise e
    else:
        raise click.ClickException(str(e))


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
        import os
        import sys
        if "pytest" not in sys.modules and "_pytest" not in sys.modules:
            os._exit(0)
    except Exception as e:
        handle_cli_exception(e)


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
        handle_cli_exception(e)


#### RESAMPlING  ####
@main.command()
@click.argument("image_path", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-s",
    "--spacing",
    type=float,
    nargs=3,
    help="Target spacing in mm (e.g. -s 4 4 4).",
)
@click.option(
    "-l",
    "--like",
    "like_path",
    type=click.Path(exists=True, path_type=Path),
    help="Reference image path to match grid layout (spacing, size, origin, direction).",
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
    help="Optional path to output the resampled image. Defaults to <name>_<spacing>mm.<ext> or <name>_like_<ref>.<ext>.",
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
    image_path, spacing, like_path, image_type, output_path, unit, preserve_total, gaussian
):
    """Resample an image to a target spacing or to match another image grid."""

    if not spacing and not like_path:
        raise click.UsageError("Either --spacing (-s) or --like (-l) must be specified.")
    if spacing and like_path:
        raise click.UsageError("Cannot specify both --spacing (-s) and --like (-l).")

    image_type_lower = image_type.lower() if image_type else None

    if output_path is None:
        suffix = "".join(image_path.suffixes)
        base = image_path.name[: -len(suffix)] if suffix else image_path.name

        if like_path:
            like_suffix = "".join(like_path.suffixes)
            like_base = like_path.name[: -len(like_suffix)] if like_suffix else like_path.name
            output_path = Path(f"{base}_like_{like_base}{suffix}")
        else:
            # Format spacing to remove trailing decimals if integers (e.g. 4.0 -> 4)
            if len(set(spacing)) == 1:
                sp_str = f"{spacing[0]:g}"
            else:
                sp_str = "x".join(f"{s:g}" for s in spacing)

            output_path = Path(f"{base}_{sp_str}mm{suffix}")

    try:
        if like_path:
            from radot.resampling import resample_image_file_like

            resample_image_file_like(
                image_path=image_path,
                like_image_path=like_path,
                output_path=output_path,
                type=image_type_lower,
                unit=unit,
                preserve_total=preserve_total,
                gaussian=gaussian,
            )
        else:
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
        handle_cli_exception(e)


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
        handle_cli_exception(e)


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
        handle_cli_exception(e)


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
        handle_cli_exception(e)


@main.command(name="dicom-list")
@click.argument("folder_path", type=click.Path(exists=True, path_type=Path))
def dicom_list(folder_path):
    """List all DICOM series in a folder and subfolders recursively, sorted by date."""
    from radot.dicom import list_dicom_series, print_dicom_series

    click.echo(f"Scanning folder recursively for DICOM series: {folder_path}\n")
    try:
        series_list = list_dicom_series(folder_path)
        print_dicom_series(series_list)
    except Exception as e:
        handle_cli_exception(e)


@main.command(name="dose")
@click.argument("image_path", type=click.Path(exists=True, path_type=Path))
@click.option(
    "-l",
    "--labels",
    "labels_path",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Labels segmentation image path.",
)
@click.option(
    "--ct",
    "ct_path",
    type=click.Path(exists=True, path_type=Path),
    help="Companion CT image path.",
)
@click.option(
    "-m",
    "--method",
    required=True,
    help="Dosimetry model / method (e.g. madsen2018).",
)
@click.option(
    "--isotope",
    required=True,
    help="Radionuclide (e.g. Lu177).",
)
@click.option(
    "--total-activity",
    "total_activity",
    type=float,
    help="Global calibration activity override (in MBq).",
)
@click.option(
    "--rois",
    multiple=True,
    help="Target ROI names. Can be specified multiple times or comma-separated.",
)
@click.option(
    "--teff",
    multiple=True,
    help="Effective half-life params per ROI (e.g. --teff Liver:45.2 --teff Kidney_L:32.1). Can be specified multiple times or comma-separated.",
)
@click.option(
    "--resample",
    type=click.Choice(["like_activity", "like_labels", "like_ct"]),
    help="Optionally resample mismatched images to a specific target grid layout.",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    type=click.Path(path_type=Path),
    help="Optional path to output dose JSON file.",
)
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help="Print detailed calculation log including all parameters and units.",
)
@click.option(
    "-q",
    "--quiet",
    is_flag=True,
    help="Silence progress/status output, showing only errors.",
)
@click.option(
    "--phantom",
    default="male",
    help="OpenDose phantom name (e.g. male, female).",
)
def dose(
    image_path,
    labels_path,
    ct_path,
    method,
    isotope,
    total_activity,
    rois,
    teff,
    output_path,
    resample,
    verbose,
    phantom,
    quiet,
):
    """Compute absorbed doses for specified regions of interest (ROI)."""
    from radot.dose import compute_roi_dose_file
    from radot.log import configure_logging

    ctx = click.get_current_context(silent=True)
    cli_quiet = quiet
    cli_verbose = verbose
    if ctx and ctx.obj:
        if not cli_quiet and ctx.obj.get("quiet"):
            cli_quiet = True
        if not cli_verbose and ctx.obj.get("verbose"):
            cli_verbose = True
    configure_logging(verbose=cli_verbose, quiet=cli_quiet)

    rois_list = []
    if rois:
        for r in rois:
            for item in r.split(","):
                item_stripped = item.strip()
                if item_stripped:
                    rois_list.append(item_stripped)
    if not rois_list:
        rois_list = None

    method_params = {"phantom": phantom}
    if teff:
        for t in teff:
            for item in t.split(","):
                item_stripped = item.strip()
                if not item_stripped:
                    continue
                if ":" not in item_stripped:
                    raise click.BadParameter(
                        f"Teff parameter '{item_stripped}' must be in ROI:value format (e.g. Liver:45.2)."
                    )
                parts = item_stripped.split(":")
                if len(parts) != 2:
                    raise click.BadParameter(
                        f"Teff parameter '{item_stripped}' must be in ROI:value format (e.g. Liver:45.2)."
                    )
                roi_name = parts[0].strip()
                try:
                    teff_val = float(parts[1].strip())
                    method_params[roi_name] = {"teff": teff_val}
                except ValueError:
                    raise click.BadParameter(
                        f"Invalid teff value: {parts[1]} for ROI: {roi_name}"
                    )

    try:
        compute_roi_dose_file(
            image_path=image_path,
            labels_path=labels_path,
            method=method,
            isotope=isotope,
            ct_path=ct_path,
            total_activity_mbq=total_activity,
            rois=rois_list,
            method_params=method_params,
            output_path=output_path,
            resample_target=resample,
            verbose=verbose,
        )
    except Exception as e:
        handle_cli_exception(e)


#### OPENDOSE DATA MAINTENANCE ####
@main.group(name="opendose", invoke_without_command=True)
@click.option(
    "-m",
    "--model",
    help="Phantom model shortcut (AM, AF, male, female, icrp_110_male, icrp_110_female).",
)
@click.option(
    "-r",
    "--radionuclide",
    help="Radionuclide name (e.g. Lu177, y-90).",
)
@click.option(
    "--force",
    is_flag=True,
    help="Force overwrite existing local S-value and mass JSON files.",
)
@click.option(
    "--upload",
    is_flag=True,
    help="Upload local database to a remote server using SCP.",
)
@click.option(
    "--upload-dest",
    help="SCP destination path (e.g. user@server:public_html/radot/opendose).",
)
@click.pass_context
def opendose(ctx, model, radionuclide, force, upload, upload_dest):
    """Download and maintain OpenDose reference masses and S-values."""
    if ctx.invoked_subcommand is not None:
        return

    if upload:
        if not upload_dest:
            raise click.UsageError("Destination path (--upload-dest) is required when using --upload.")
        
        # Normalize protocol separator if user wrote host://path instead of host:path
        if "://" in upload_dest:
            upload_dest = upload_dest.replace("://", ":")
            
        import subprocess
        import importlib.resources
        from pathlib import Path
        
        db_dir = importlib.resources.files("radot").joinpath("data").joinpath("opendose")
        db_path = Path(str(db_dir))
        
        # Get all subdirectories and files inside the database folder to upload them directly
        src_paths = [str(p) for p in db_path.iterdir() if not p.name.startswith(".")]
        if not src_paths:
            raise click.ClickException("No database files found to upload.")
            
        cmd = ["scp", "-r"] + src_paths + [upload_dest]
        click.echo(f"RADOT: Uploading local database contents using: {' '.join(cmd)}")
        try:
            subprocess.run(cmd, check=True)
            click.echo("RADOT: Database upload completed successfully!")
        except Exception as e:
            raise click.ClickException(f"Failed to upload database via SCP: {e}")
        return

    # Normal downloading validation
    if not model:
        raise click.UsageError("Missing option '-m' / '--model'.")
    if not radionuclide:
        raise click.UsageError("Missing option '-r' / '--radionuclide'.")

    from radot.db import download_opendose_data

    try:
        download_opendose_data(
            phantom=model,
            isotope=radionuclide,
            force=force,
        )
    except Exception as e:
        handle_cli_exception(e)


@opendose.command(name="status")
@click.option(
    "-m",
    "--model",
    help="Phantom model shortcut (AM, AF, male, female, icrp_110_male, icrp_110_female).",
)
@click.option(
    "-r",
    "--radionuclide",
    help="Radionuclide name (e.g. Lu177, y-90).",
)
def opendose_status(model, radionuclide):
    """List local OpenDose database inventory or print downloaded values."""
    import importlib.resources
    from pathlib import Path
    import json
    import gzip
    from radot.db import normalize_phantom, normalize_isotope

    db_dir = importlib.resources.files("radot").joinpath("data").joinpath("opendose")
    db_path = Path(str(db_dir))

    # If model is provided, check if we just want to print details of that model/radionuclide
    if model:
        try:
            norm_phantom = normalize_phantom(model)
        except ValueError as e:
            handle_cli_exception(e)

        phantom_path = db_path / norm_phantom
        if not phantom_path.exists():
            click.echo(f"Phantom '{norm_phantom}' is not downloaded/present locally.")
            return

        if radionuclide:
            norm_iso = normalize_isotope(radionuclide)
            file_gz = phantom_path / f"{norm_iso}.json.gz"
            file_json = phantom_path / f"{norm_iso}.json"
            
            target_file = None
            if file_gz.exists():
                target_file = file_gz
            elif file_json.exists():
                target_file = file_json
                
            if not target_file:
                click.echo(f"S-value matrix for isotope '{norm_iso}' under phantom '{norm_phantom}' is not downloaded/present locally.")
                return
                
            try:
                if target_file.suffix == ".gz":
                    with gzip.open(target_file, "rt", encoding="utf-8") as f:
                        s_matrix = json.load(f)
                else:
                    with open(target_file, "r", encoding="utf-8") as f:
                        s_matrix = json.load(f)
            except Exception as e:
                raise click.ClickException(f"Failed to read S-value matrix: {e}")
                
            click.echo(f"S-values for Phantom: {norm_phantom}, Radionuclide: {norm_iso} (Path: {target_file})")
            click.echo(f"{'Source Organ':<25} {'Target Organ':<25} {'S-value (mGy/MBq.s)':<20}")
            click.echo("-" * 72)
            for src, targets in sorted(s_matrix.items()):
                for tgt, val in sorted(targets.items()):
                    click.echo(f"{src:<25} {tgt:<25} {val:<20.4e}")
            return
        else:
            masses_file = phantom_path / "masses.json"
            if not masses_file.exists():
                click.echo(f"Reference masses for phantom '{norm_phantom}' are not downloaded/present locally.")
                return
                
            try:
                with open(masses_file, "r", encoding="utf-8") as f:
                    masses = json.load(f)
            except Exception as e:
                raise click.ClickException(f"Failed to read masses: {e}")
                
            click.echo(f"Reference Masses for Phantom: {norm_phantom} (Path: {masses_file})")
            click.echo(f"{'Organ':<30} {'Mass (g)':<15}")
            click.echo("-" * 45)
            for organ, mass in sorted(masses.items()):
                click.echo(f"{organ:<30} {mass:<15.3f}")
            return

    if radionuclide and not model:
        raise click.UsageError("Missing option '-m' / '--model' when using '-r' / '--radionuclide'.")

    # Print header
    click.echo("Local OpenDose Database Inventory:")
    click.echo(f"Path: {db_path}\n")

    if not db_path.exists():
        click.echo("Total space used: 0 B (0 Radionuclides, 0 Phantoms)")
        return

    # Find phantoms (subdirectories)
    phantoms = sorted([p for p in db_path.iterdir() if p.is_dir() and not p.name.startswith(".")])

    total_bytes = 0
    phantom_count = 0
    radionuclide_count = 0

    for i, phantom in enumerate(phantoms):
        is_last_phantom = (i == len(phantoms) - 1)
        phantom_prefix = "└── " if is_last_phantom else "├── "
        click.echo(f"{phantom_prefix}{phantom.name}/")
        
        # Find files in phantom
        files = sorted([f for f in phantom.iterdir() if f.is_file() and not f.name.startswith(".")])
        # Filter to only masses.json and S-value matrices
        files = [f for f in files if f.name == "masses.json" or f.name.endswith(".json.gz") or f.name.endswith(".json")]
        
        file_indent = "    " if is_last_phantom else "│   "
        
        phantom_has_files = False
        for j, file in enumerate(files):
            is_last_file = (j == len(files) - 1)
            file_prefix = "└── " if is_last_file else "├── "
            
            # Read details of the file
            details = ""
            try:
                if file.name == "masses.json":
                    import json
                    with open(file, "r", encoding="utf-8") as f:
                        data = json.load(f)
                    details = f" ({len(data)} organs)"
                elif file.name.endswith(".json.gz"):
                    import gzip
                    import json
                    with gzip.open(file, "rt", encoding="utf-8") as f:
                        data = json.load(f)
                    num_sources = len(data)
                    num_svalues = sum(len(v) for v in data.values())
                    details = f" ({num_sources} source organs, {num_svalues} S-values)"
                elif file.name.endswith(".json"):
                    import json
                    with open(file, "r", encoding="utf-8") as f:
                        data = json.load(f)
                    num_sources = len(data)
                    num_svalues = sum(len(v) for v in data.values())
                    details = f" ({num_sources} source organs, {num_svalues} S-values)"
            except Exception as e:
                details = f" (error: {e})"
                
            click.echo(f"{file_indent}{file_prefix}[OK] {file.name}{details}")
            
            # Track statistics
            total_bytes += file.stat().st_size
            if file.name != "masses.json":
                radionuclide_count += 1
            phantom_has_files = True
            
        if phantom_has_files:
            phantom_count += 1

    # Format space used
    if total_bytes >= 1024 * 1024:
        space_str = f"{total_bytes / (1024 * 1024):.1f} MB"
    elif total_bytes >= 1024:
        space_str = f"{total_bytes / 1024:.1f} KB"
    else:
        space_str = f"{total_bytes} B"

    click.echo(f"\nTotal space used: {space_str} ({radionuclide_count} Radionuclides, {phantom_count} Phantoms)")


@opendose.command(name="list")
@click.option(
    "-m",
    "--model",
    help="Phantom model shortcut (AM, AF, male, female, icrp_110_male, icrp_110_female).",
)
@click.option(
    "-r",
    "--radionuclide",
    help="Radionuclide name (e.g. Lu177, y-90).",
)
@click.pass_context
def opendose_list(ctx, model, radionuclide):
    """List local OpenDose database inventory or print downloaded values."""
    ctx.forward(opendose_status)

@main.group(name="roi")
def roi():
    """Nomenclature mapping and ROI verification tools."""
    pass


@roi.command(name="list")
@click.option(
    "--source",
    type=click.Choice(["all", "ts", "od", "tg"]),
    default="all",
    help="Filter mappings by origin ('all' for TG-263 Hub, 'ts' for TotalSegmentator checklist, 'od' for OpenDose phantom capability, 'tg' for AAPM compliance).",
)
def roi_list(source):
    """List all explicit anatomical correspondences in the nomenclature mapping."""
    import re
    from radot.nomenclature import (
        TS_TO_TG263,
        OD_TO_TG263,
        TS_TO_OD_FALLBACK,
        TG263_TO_OD,
        TS_GROUPS,
        OPENDOSE_PASSIVE_TARGETS,
    )
    
    if source == "all":
        click.echo("\nRegistered ROI Nomenclature Maps (TG-263 Hub)\n")
        click.echo(f"{'TotalSegmentator Name':<30} {'TG-263 Label':<20} {'OpenDose Target':<25}")
        click.echo("─" * 78)
        
        for ts_name, val in sorted(TS_TO_TG263.items()):
            tg263_label = val[1]
            od_target = TG263_TO_OD.get(tg263_label, "")
            click.echo(f"{ts_name:<30} {tg263_label:<20} {od_target:<25}")
        click.echo()
        
    elif source == "ts":
        click.echo("\nTotalSegmentator Mapping Inventory\n")
        click.echo(f"{'TotalSegmentator Name':<30} {'TG-263 Label':<20} {'Mapping Route':<20} {'OpenDose Target':<25}")
        click.echo("─" * 98)
        
        for ts_name, val in sorted(TS_TO_TG263.items()):
            tg263_label = val[1]
            od_target = TG263_TO_OD.get(tg263_label, "")
            
            if od_target:
                route = "[STAGE 1] Hub"
                resolves_to = od_target
            elif ts_name in TS_TO_OD_FALLBACK:
                route = "[STAGE 2] Bypass"
                resolves_to = TS_TO_OD_FALLBACK[ts_name]
            elif ts_name in TS_GROUPS:
                route = "[STAGE 2] Group"
                resolves_to = TS_GROUPS[ts_name]
            else:
                route = "[STAGE 3] UNMAPPED"
                resolves_to = "⛔ CRASH SAFEGUARD"
                
            click.echo(f"{ts_name:<30} {tg263_label:<20} {route:<20} {resolves_to:<25}")
        click.echo()
        
    elif source == "od":
        # Build mapping from OpenDose target -> list of TS Names
        od_to_ts_names = {}
        for ts_name, val in TS_TO_TG263.items():
            tg263_label = val[1]
            od_target = TG263_TO_OD.get(tg263_label, "")
            if od_target:
                od_to_ts_names.setdefault(od_target, []).append(ts_name)
                
        for ts_name, od_val in TS_TO_OD_FALLBACK.items():
            od_to_ts_names.setdefault(od_val, []).append(ts_name)
            
        for ts_name, od_val in TS_GROUPS.items():
            od_to_ts_names.setdefault(od_val, []).append(ts_name)
            
        from radot.db import _find_matching_key, _load_masses
        
        try:
            masses_male = _load_masses("icrp_110_male")
            masses_female = _load_masses("icrp_110_female")
            standard_od_organs = set(masses_male.keys()) | set(masses_female.keys())
        except Exception:
            standard_od_organs = set(OD_TO_TG263.keys()) | set(TS_TO_OD_FALLBACK.values()) | set(OPENDOSE_PASSIVE_TARGETS)
            
        full_od_to_ts_names = {}
        for od_organ in standard_od_organs:
            matched_base = _find_matching_key(od_to_ts_names, od_organ)
            if matched_base:
                full_od_to_ts_names[od_organ] = od_to_ts_names[matched_base]
            else:
                for base_target, ts_list in od_to_ts_names.items():
                    if od_organ.lower().startswith(base_target.lower()) or base_target.lower().startswith(od_organ.lower()):
                        full_od_to_ts_names[od_organ] = ts_list
                        break
        
        click.echo("\nOpenDose Reference Phantom Targets\n")
        click.echo(f"{'OpenDose Organ Name':<50} {'Status':<20} {'Mapped From (TS)':<25}")
        click.echo("─" * 98)
        
        for od_name in sorted(standard_od_organs):
            ts_list = full_od_to_ts_names.get(od_name, [])
            if ts_list:
                status = "[ACTIVE]"
                mapped_from = ", ".join(ts_list)
                if len(mapped_from) > 40:
                    mapped_from = mapped_from[:37] + "..."
            elif od_name in OPENDOSE_PASSIVE_TARGETS:
                status = "[PASSIVE TARGET]"
                mapped_from = "-- (Cross-fire only)"
            else:
                is_passive = False
                for p_target in OPENDOSE_PASSIVE_TARGETS:
                    if od_name.lower().startswith(p_target.lower()) or p_target.lower().startswith(od_name.lower()):
                        is_passive = True
                        break
                if is_passive:
                    status = "[PASSIVE TARGET]"
                    mapped_from = "-- (Cross-fire only)"
                else:
                    status = "[INACTIVE]"
                    mapped_from = "--"
                
            click.echo(f"{od_name:<50} {status:<20} {mapped_from:<25}")
        click.echo()
        
    elif source == "tg":
        tg263_labels = set(TG263_TO_OD.keys()) | {val[1] for val in TS_TO_TG263.values()}
        
        click.echo("\nSupported AAPM TG-263 Standard Labels\n")
        click.echo(f"{'TG-263 Label':<25} {'Description':<35} {'Status':<15}")
        click.echo("─" * 78)
        
        for label in sorted(tg263_labels):
            desc = label.replace("_", " ")
            od_name = TG263_TO_OD.get(label, "")
            if od_name:
                desc = f"TG-263 standard for {od_name.lower()}"
            else:
                desc = f"TG-263 standard for {desc.lower()}"
                
            click.echo(f"{label:<25} {desc:<35} {'Registered':<15}")
        click.echo()


@roi.command(name="match")
@click.argument("name")
def roi_match(name):
    """Query an ROI name to find exact matches or closest nomenclature suggestions."""
    from radot.nomenclature import resolve_roi_name_to_opendose, find_closest_matches, TS_TO_TG263, TG263_TO_OD, TS_TO_OD_FALLBACK
    
    lower_name = name.lower()
    is_registered = False
    
    if lower_name in [k.lower() for k in TS_TO_TG263.keys()]:
        is_registered = True
    elif lower_name in [k.lower() for k in TG263_TO_OD.keys()]:
        is_registered = True
    elif lower_name in [k.lower() for k in TG263_TO_OD.values()]:
        is_registered = True
    elif lower_name in [k.lower() for k in TS_TO_OD_FALLBACK.keys()]:
        is_registered = True
        
    if is_registered:
        resolved = resolve_roi_name_to_opendose(name)
        click.echo(f"\n✅ \"{name}\" is registered. Resolves to OpenDose: \"{resolved}\"\n")
    else:
        click.echo(f"\n❌ \"{name}\" is not explicitly registered.")
        closest = find_closest_matches(name)
        if closest:
            click.echo("💡 Did you mean one of these?")
            for match_name, source in closest:
                click.echo(f"   - {match_name}  ({source})")
        click.echo()


@roi.command(name="verify")
@click.argument("path", type=click.Path(exists=True, file_okay=True, dir_okay=False))
def roi_verify(path):
    """Scan a patient sidecar JSON mapping file before execution to report unmapped structures."""
    from radot.nomenclature import verify_sidecar_mappings
    import json
    
    try:
        with open(path, "r", encoding="utf-8") as f:
            raw_map = json.load(f)
            labels_map = {int(k): v for k, v in raw_map.items()}
    except Exception as e:
        raise click.ClickException(f"Failed to parse label mapping JSON file: {e}")
        
    report = verify_sidecar_mappings(labels_map)
    
    click.echo(f"\nScanning sidecar mappings from: {path}\n")
    
    if report["mapped"]:
        click.echo("Mapped structures:")
        for lid, (name, resolved) in sorted(report["mapped"].items()):
            click.echo(f"   - Label {lid}: \"{name}\" -> \"{resolved}\" (OpenDose)")
            
    if report["unmapped"]:
        click.echo("\n❌ Unmapped structures (requires override or manual S-value definition):")
        for lid, name in sorted(report["unmapped"].items()):
            click.echo(f"   - Label {lid}: \"{name}\"")
        click.echo()
        raise click.ClickException("Verification failed: sidecar contains unmapped structures.")
    else:
        click.echo("\n✅ Verification successful: all structures successfully mapped!")
        click.echo()


@main.group(name="radionuclide")
def radionuclide():
    """Radionuclide information tools."""
    pass


@radionuclide.command(name="list")
def radionuclide_list():
    """List all available radionuclides in the ICRP-107 database."""
    import importlib.resources
    import sys
    
    try:
        icrp_dir = importlib.resources.files("icrp107_database").joinpath("icrp107")
        available = []
        for p in icrp_dir.iterdir():
            if p.name.endswith(".json"):
                available.append(p.name[:-5])
        
        if not available:
            sys.stdout.write("No radionuclides found in the icrp107-database.\n")
            return
            
        available.sort()
        sys.stdout.write("\nAvailable Radionuclides:\n")
        # Format as 6 columns
        col_width = 12
        num_cols = 6
        for i in range(0, len(available), num_cols):
            chunk = available[i:i + num_cols]
            line = "".join(f"{name:<{col_width}}" for name in chunk) + "\n"
            sys.stdout.write(line)
        sys.stdout.write("\n")
    except Exception as e:
        handle_cli_exception(e)


@radionuclide.command(name="info")
@click.argument("isotope")
def radionuclide_info(isotope):
    """Print physical half-life and emissions summary for a radionuclide."""
    import importlib.resources
    import json
    import re
    import sys
    
    try:
        iso_clean = re.sub(r'[-_\s]', '', isotope).lower()
        icrp_dir = importlib.resources.files("icrp107_database").joinpath("icrp107")
        
        matched_path = None
        for p in icrp_dir.iterdir():
            if p.name.endswith(".json"):
                name_without_ext = p.name[:-5]
                if re.sub(r'[-_\s]', '', name_without_ext).lower() == iso_clean:
                    matched_path = p
                    break
                    
        if not matched_path:
            raise ValueError(f"Isotope '{isotope}' is not available in the icrp107-database package.")
            
        with matched_path.open("r", encoding="utf-8") as f:
            data = json.load(f)
            if isinstance(data, str):
                data = json.loads(data)
                
        name = data.get("name", matched_path.name[:-5])
        half_life = data.get("half_life")
        time_unit = data.get("time_unit", "s")
        
        # Display physical half-life with original units and converted units
        from radot.db import get_physical_half_life_h
        try:
            hl_hours = get_physical_half_life_h(name)
            hl_days = hl_hours / 24.0
            if hl_days >= 1.0:
                hl_str = f"{half_life} {time_unit} ({hl_days:.3f} days / {hl_hours:.2f} hours)"
            else:
                hl_str = f"{half_life} {time_unit} ({hl_hours:.4f} hours)"
        except Exception:
            hl_str = f"{half_life} {time_unit}"
            
        sys.stdout.write(f"\nRadionuclide: {name}\n")
        sys.stdout.write(f"Physical Half-life: {hl_str}\n\n")
        
        emissions = data.get("emissions", {})
        if not emissions:
            sys.stdout.write("No emissions data available for this radionuclide.\n")
            return
            
        sys.stdout.write("Emissions Summary:\n")
        sys.stdout.write(f"{'Type':<12} | {'Num':<5} | {'Total Yield':<11} | {'Total Energy (MeV)':<18} | {'Avg Energy (MeV)':<16}\n")
        sys.stdout.write("─" * 73 + "\n")
        
        # We order standard emission types
        emission_types = ["alpha", "beta-", "beta+", "gamma", "X", "auger", "IE", "neutron", "alpha recoil", "annihilation", "fission", "betaD"]
        for etype in emission_types:
            val = emissions.get(etype, [])
            if not val:
                continue
                
            num = len(val)
            total_yield = 0.0
            total_energy = 0.0
            for item in val:
                if len(item) >= 2:
                    energy, yield_val = item[0], item[1]
                    total_yield += yield_val
                    total_energy += energy * yield_val
            
            avg_energy = total_energy / total_yield if total_yield > 0 else 0.0
            sys.stdout.write(f"{etype:<12} | {num:<5} | {total_yield:<11.4f} | {total_energy:<18.4f} | {avg_energy:<16.4f}\n")
        sys.stdout.write("\n")
        
    except Exception as e:
        handle_cli_exception(e)


if __name__ == "__main__":
    main()

