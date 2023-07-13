from LTEpy import lte
import pytest

def test_set_temp():
    """Test changing the temperature.

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
    """Test computing B_nu.

    """

    planck = lte.Planck(1000)  # new instance of Planck

    nu = 6e14  # frequency in Hz
    B_nu_expected = 9.94e-16  # from Wolfram Alpha
    B_nu = planck.compute_B_nu(nu)

    assert B_nu == pytest.approx(B_nu_expected, abs=0.1*B_nu_expected)


def test_compute_B_lambda():
    """Test computing B_lambda.

    """

    planck = lte.Planck(1000)  # new instance of Planck

    wl = 5e-5  # wavelength in cm
    B_lambda_expected = 0.001213  # from Wolfram Alpha
    B_lambda = planck.compute_B_lambda(wl)

    assert B_lambda == pytest.approx(B_lambda_expected, abs=0.1*B_lambda_expected)


if __name__ == '__main__':
    test_set_temp()
    test_compute_B_nu()
    test_compute_B_lambda()