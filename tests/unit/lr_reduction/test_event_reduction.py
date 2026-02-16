from unittest.mock import MagicMock

import mantid.simpleapi as mtd_api
import numpy as np
import pytest

from lr_reduction import event_reduction
from lr_reduction.utils import amend_config


@pytest.fixture
def make_mock_ws():
    def _make_mock_ws(properties=None, xi_reference=445, start_time="2025-07-14T10:00:00"):
        properties = properties or {}

        mock_instrument = MagicMock()
        mock_instrument.hasParameter.return_value = True
        mock_instrument.getNumberParameter.return_value = [xi_reference]

        def get_property_side_effect(name):
            val = properties.get(name)
            if val is None:
                # Mantid raises RuntimeError if property not found for workspace run
                raise RuntimeError(f"Property {name} not in mock properties")
            mock = MagicMock()
            mock.value = [val]  # keep consistent with .value[0]
            return mock

        mock_run = MagicMock()
        mock_run.getProperty.side_effect = get_property_side_effect

        mock_start_time = MagicMock()
        mock_start_time.value = start_time
        mock_run.__getitem__.side_effect = lambda key: mock_start_time if key == "start_time" else None

        ws = MagicMock()
        ws.getInstrument.return_value = mock_instrument
        ws.getRun.return_value = mock_run

        return ws

    return _make_mock_ws


def test_trapezoidal_cdf_analytic_basic():
    x = np.array([-1, 0, 1])
    L_, l_ = 1.0, 0.5
    cdf = event_reduction._trapezoidal_cdf_analytic(x, L_, l_)
    assert cdf.shape == x.shape
    assert np.all(cdf >= 0) and np.all(cdf <= 1)
    assert cdf[0] == 0
    assert cdf[-1] == 1


def test_find_sigma_68_basic():
    L_, l_ = 1.0, 0.5
    sigma = event_reduction._find_sigma_68(L_, l_)
    assert 0 < sigma < L_
    assert sigma == pytest.approx(0.5101, rel=1e-3)


def test_trapezoid_gaussian_function(make_mock_ws):
    """
    Test the trapezoidal distribution params function in event_reduction using mock instrument parameters.
    """
    mock_run_properties = {"ths": 0.369, "S1VHeight": 0.29, "SiVHeight": 0.145, "BL4B:Mot:xi.RBV": 300}

    ws = make_mock_ws(properties=mock_run_properties)

    L_bottom, l_top, sigma_equiv, dth_over_th = event_reduction.trapezoidal_distribution_params(ws)
    print(L_bottom, l_top, sigma_equiv, dth_over_th)
    # Calculate expected values independently here or reuse logic from event_reduction
    # For simplicity, just check types and ranges:
    assert isinstance(L_bottom, float)
    assert isinstance(l_top, float)
    assert isinstance(sigma_equiv, float)
    assert sigma_equiv == pytest.approx(0.00434, rel=1e-3)
    assert dth_over_th == pytest.approx(0.01176, rel=1e-3)
    assert L_bottom == pytest.approx(0.0093, rel=1e-3)
    assert l_top == pytest.approx(0.0031, rel=1e-3)


@pytest.mark.datarepo
@pytest.mark.parametrize("q_summing, theta_expected, dtheta_expected", [(False, 1.1825526, 0.011999676), (True, 1.1825526, 0.031941999)])
def test_compute_theta_resolution(q_summing, theta_expected, dtheta_expected, nexus_dir):
    """Test the compute_theta_resolution function."""
    with amend_config(data_dir=nexus_dir):
        ws = mtd_api.LoadEventNexus("REF_L_201287")

    theta, dtheta = event_reduction.compute_theta_resolution(ws, q_summing=q_summing)
    assert theta == pytest.approx(theta_expected, rel=1e-6)
    assert dtheta == pytest.approx(dtheta_expected, rel=1e-6)


def test_compute_wavelength_resolution():
    """Test the compute_wavelength_resolution function."""
    resolution_function_str = "name=UserFunction, Formula=L - A * exp(-k * x), L=0.07564423, A=0.13093263, k=0.34918918"

    wavelengths = np.array([11.4254409079602, 11.20141265486294, 10.981777112610725, 10.766448149618357, 10.555341323155254, 10.34837384623064, 10.145464555128079, 9.946533877576549, 9.751503801545635, 9.560297844652583, 9.3728410241692, 9.189059827616862, 9.0088821839381, 8.83223743523343, 8.659056309052383, 8.489270891227827, 8.322814599242967, 8.159622156120555, 7.999629564824074, 7.842774083160856, 7.688994199177309])

    # Calculate expected resolution using the default formula:
    # L - A * exp(-k * x), L=0.07564423, A=0.13093263, k=0.34918918
    A = 0.13093263
    L = 0.07564423
    k = 0.34918918
    expected_resolution = L - A * np.exp(-k * wavelengths)

    _, wavelength_resolution = event_reduction.compute_wavelength_resolution(wavelengths, resolution_function_str)
    assert len(wavelength_resolution) == len(wavelengths)
    np.testing.assert_array_almost_equal(wavelength_resolution, expected_resolution, decimal=5)
