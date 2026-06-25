import pytest
from radot.db import (
    get_physical_half_life_h,
    normalize_phantom,
    normalize_isotope,
    resolve_organ_name,
    get_reference_mass,
    get_s_value
)

def test_get_physical_half_life():
    # Test standard known isotopes
    assert get_physical_half_life_h("Lu177") == pytest.approx(159.528)
    assert get_physical_half_life_h("lu-177") == pytest.approx(159.528)
    assert get_physical_half_life_h("Y90") == pytest.approx(64.1)
    assert get_physical_half_life_h("y-90") == pytest.approx(64.1)
    assert get_physical_half_life_h("I131") == pytest.approx(192.4968)
    assert get_physical_half_life_h("i-131") == pytest.approx(192.4968)
    # Ac225 should work now that icrp107-database is required
    assert get_physical_half_life_h("Ac225") == pytest.approx(240.0)

def test_unknown_isotope_error():
    # Test that a truly unknown isotope raises ValueError
    with pytest.raises(ValueError, match="is not available in the icrp107-database package"):
        get_physical_half_life_h("Xyz999")

def test_normalize_phantom():
    assert normalize_phantom("AM") == "icrp_110_male"
    assert normalize_phantom("male") == "icrp_110_male"
    assert normalize_phantom("icrp_110_male") == "icrp_110_male"
    assert normalize_phantom("AF") == "icrp_110_female"
    assert normalize_phantom("female") == "icrp_110_female"
    assert normalize_phantom("icrp_110_female") == "icrp_110_female"
    
    with pytest.raises(ValueError, match="Unknown phantom"):
        normalize_phantom("unknown_phantom")

def test_normalize_isotope():
    assert normalize_isotope("Lu-177") == "lu177"
    assert normalize_isotope("lu177") == "lu177"
    assert normalize_isotope("Y_90") == "y90"

def test_resolve_organ_name():
    # Test TG-263 nomenclature mappings
    assert resolve_organ_name("Kidney_L") == "Left kidney"
    assert resolve_organ_name("kidney_r") == "Right kidney"
    assert resolve_organ_name("bladder") == "Urinary bladder wall"
    assert resolve_organ_name("liver") == "Liver"
    
    # Test fallback formatting
    assert resolve_organ_name("some_random_organ") == "Some random organ"

def test_get_reference_mass():
    # Check male masses from fake database
    assert get_reference_mass("male", "Kidney_L") == pytest.approx(152.983)
    assert get_reference_mass("AM", "liver") == pytest.approx(1799.776)
    
    # Check female masses from fake database
    assert get_reference_mass("female", "Kidney_R") == pytest.approx(125.530)
    assert get_reference_mass("AF", "spleen") == pytest.approx(129.996)

def test_get_s_value():
    # Check male S-values from fake matrix
    assert get_s_value("male", "lu177", "left_kidney", "right_kidney") == pytest.approx(9.0850e-8)
    assert get_s_value("AM", "Y90", "liver", "spleen") == pytest.approx(2.1234e-7)
    
    # Check female S-values from fake matrix
    assert get_s_value("female", "lu177", "left_kidney", "left_kidney") == pytest.approx(1.6034e-4)
    assert get_s_value("AF", "y90", "liver", "liver") == pytest.approx(6.9876e-5)

def test_db_errors():
    with pytest.raises(ValueError, match="not found in phantom"):
        get_reference_mass("male", "non_existent_organ")
        
    with pytest.raises(ValueError, match="not found in S-value matrix"):
        get_s_value("male", "lu177", "non_existent_organ", "left_kidney")
        
    with pytest.raises(ValueError, match="not found in S-value matrix"):
        get_s_value("male", "lu177", "left_kidney", "non_existent_organ")

def test_cli_opendose():
    from click.testing import CliRunner
    from radot.cli import main
    
    runner = CliRunner()
    
    # Test single isotope download (which runs either real or fallback offline downloader)
    result = runner.invoke(main, ["opendose", "-m", "AM", "-r", "Lu177", "--force"])
    assert result.exit_code == 0
    assert "icrp_110_male" in result.output
    
    # Test group core download
    result = runner.invoke(main, ["opendose", "-m", "AF", "--group", "core", "--force"])
    assert result.exit_code == 0
    assert "icrp_110_female" in result.output
    
    # Test error cases
    result = runner.invoke(main, ["opendose", "-m", "AM"])
    assert result.exit_code != 0
    assert "Either --radionuclide (-r) or --group <group> is required" in result.output


