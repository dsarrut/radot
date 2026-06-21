import click
from pathlib import Path

@click.group()
@click.version_option(version="0.1.0")
def main():
    """radot: Radiation Dosimetry for Therapy command-line tools."""
    pass

@main.command()
@click.option("--name", default="World", help="Name to greet.")
def hello(name):
    """Simple hello world command to verify setup."""
    click.echo(f"Hello, {name}! Welcome to radot.")

@main.command()
@click.argument("ct_path", type=click.Path(exists=True, path_type=Path))
@click.option("-o", "--output", "output_path", type=click.Path(path_type=Path), help="Output label image path. Defaults to labels_all.nii.gz or labels_head.nii.gz depending on mode.")
@click.option("--mode", default="total", help="TotalSegmentator mode (e.g. total, head_glands_cavities, body, etc.).")
@click.option("--fast", is_flag=True, help="Run TotalSegmentator in fast mode (lower resolution).")
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

@main.command()
@click.argument("labels_path", type=click.Path(exists=True, path_type=Path))
@click.option("-i", "--image", "image_path", required=True, type=click.Path(exists=True, path_type=Path), help="Input image (CT/PET/SPECT) to compute stats on.")
@click.option("-t", "--type", "image_type", required=True, type=click.Choice(["ct", "pet", "spect"], case_sensitive=False), help="Image modality type.")
@click.option("-o", "--output", "output_path", type=click.Path(path_type=Path), help="Optional path to output stats as JSON. Defaults to printing to console.")
@click.option("--labels", help="Comma-separated list of label IDs to filter (e.g. 5,12).")
@click.option("--rois", help="Comma-separated list of ROI names to filter (e.g. Kidney_L,Spleen).")
def stats(labels_path, image_path, image_type, output_path, labels, rois):
    """Compute volume and intensity statistics for labels."""
    from radot.stats import compute_stats
    
    # Parse labels list
    labels_list = None
    if labels:
        try:
            labels_list = [int(x.strip()) for x in labels.split(",") if x.strip()]
        except ValueError:
            raise click.BadParameter("Labels must be a comma-separated list of integers.")

    # Parse rois list
    rois_list = None
    if rois:
        rois_list = [x.strip() for x in rois.split(",") if x.strip()]

    try:
        compute_stats(labels_path, image_path, image_type.lower(), output_path, labels_list, rois_list)
    except Exception as e:
        raise click.ClickException(str(e))

if __name__ == "__main__":
    main()
