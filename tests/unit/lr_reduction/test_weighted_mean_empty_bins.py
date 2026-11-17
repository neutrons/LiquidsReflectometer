"""
The scale between two overlapping curves is found from the points that hold a measurement.

A Q bin that no neutron fell in has R = 0 and dR = 0. The ratio of the curves there is a division by zero, and one
such point used to turn the weighted mean into NaN: the caller then gave up and left the second curve unscaled.
"""

import warnings

import numpy as np
import pytest

from lr_reduction import nr_tools


def curves(scale=1.25, points=12):
    y2 = np.linspace(2.0, 1.0, points)
    e2 = y2 * 0.02
    return y2 * scale, y2, e2 * scale, e2


def test_an_empty_bin_in_the_second_curve_does_not_spoil_the_scale():
    y1, y2, e1, e2 = curves()
    expected, _ = nr_tools.weighted_mean(y1, y2, e1, e2)
    y2[4], e2[4] = 0.0, 0.0

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        mean, sigma = nr_tools.weighted_mean(y1, y2, e1, e2)

    assert mean == pytest.approx(expected) == pytest.approx(1.25)
    assert np.isfinite(sigma)


def test_an_empty_bin_in_the_first_curve_is_left_out_too():
    y1, y2, e1, e2 = curves()
    y1[7], e1[7] = 0.0, 0.0

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        mean, _ = nr_tools.weighted_mean(y1, y2, e1, e2)

    assert mean == pytest.approx(1.25)


def test_curves_without_empty_bins_give_what_they_gave_before():
    y1, y2, e1, e2 = curves()
    y1[3] *= 1.02  # not all ratios equal, so that the outlier step has something to weigh

    v = y1 / y2
    w = 1 / ((e1 / y2) ** 2 + (e2 * y1 / y2**2) ** 2)
    first = np.sum(v * w) / np.sum(w)
    kept = np.abs(v - first) < 3 * np.std(v - first)
    expected = np.sum(v[kept] * w[kept]) / np.sum(w[kept])

    mean, sigma = nr_tools.weighted_mean(y1, y2, e1, e2)

    assert mean == pytest.approx(expected, rel=1e-12)
    assert sigma == pytest.approx(np.sqrt(1 / np.sum(w[kept])), rel=1e-12)


def test_no_measured_point_at_all_says_so_without_a_warning():
    y1, y2, e1, e2 = curves()

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        mean, sigma = nr_tools.weighted_mean(y1, np.zeros_like(y2), e1, np.zeros_like(e2))

    assert np.isnan(mean)
    assert np.isinf(sigma)
