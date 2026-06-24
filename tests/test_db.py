import pytest
from radot.db import get_physical_half_life_h

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

