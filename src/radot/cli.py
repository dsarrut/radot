import click

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

if __name__ == "__main__":
    main()
