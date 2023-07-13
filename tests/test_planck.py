from LTEpy import lte
import pytest

def test_set_temp():
    """Test setting the temperature.

    """

    planck = lte.Planck(1000)  # new instance of Planck

    # try changing to a positive temp
    new_temp = 5000
    planck.set_temp(new_temp)
    assert planck.temp == pytest.approx(new_temp, abs=1e-3)

    # try changing to a negative temp
    with pytest.raises(ValueError):
        new_temp = -1000
        planck.set_temp(new_temp)

    # try changing to zero temp
    with pytest.raises(ValueError):
        new_temp = 0
        planck.set_temp(new_temp)


def test_compute_B_nu():
    """Test computation of B_nu.

    """

    planck = lte.Planck(1000)  # new instance of Planck



if __name__ == '__main__':
    test_set_temp()