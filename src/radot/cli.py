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
@click.option("-o", "--output", "output_path", required=True, type=click.Path(path_type=Path), help="Output label image path.")
@click.option("--mode", default="total", help="TotalSegmentator mode (e.g. total, head_glands_cavities, body, etc.).")
@click.option("--fast", is_flag=True, help="Run TotalSegmentator in fast mode (lower resolution).")
def segment(ct_path, output_path, mode, fast):
    """Segment CT and group/name labels according to AAPM TG-263."""
    from radot.segmentation import segment_ct
    try:
        segment_ct(ct_path, output_path, mode, fast)
    except Exception as e:
        raise click.ClickException(str(e))

if __name__ == "__main__":
    main()
