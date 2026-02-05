from unittest.mock import MagicMock

import numpy as np
import pytest

from lr_reduction import event_reduction


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
