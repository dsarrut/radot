import pytest
from click.testing import CliRunner
from radot.cli import main

def test_radionuclide_list():
    runner = CliRunner()
    result = runner.invoke(main, ["radionuclide", "list"])
    assert result.exit_code == 0
    assert "Lu-177" in result.output
    assert "Y-90" in result.output


def test_radionuclide_info():
    runner = CliRunner()
    result = runner.invoke(main, ["radionuclide", "info", "Lu-177"])
    assert result.exit_code == 0
    assert "Radionuclide: Lu-177" in result.output
    assert "Physical Half-life" in result.output
    assert "beta-" in result.output
    assert "gamma" in result.output


def test_radionuclide_info_invalid():
    runner = CliRunner()
    result = runner.invoke(main, ["radionuclide", "info", "invalid_iso"])
    assert result.exit_code != 0
    assert "is not available in the icrp107-database" in result.output
