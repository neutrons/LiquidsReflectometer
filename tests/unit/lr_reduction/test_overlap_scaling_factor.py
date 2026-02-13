"""Test OverlapScalingFactor class methods."""

import numpy as np
import pytest
from pytest import approx

from lr_reduction.scaling_factors import OverlapScalingFactor, ReducedData


# Create some fake linear data for testing
@pytest.fixture
def fake_data():
    q = np.linspace(0.01, 0.2, 100)
    r = np.linspace(2, 2, 100)
    err = 0.05 * np.ones_like(r)  # Constant error

    return ReducedData(q=q, r=r, err=err)


@pytest.fixture
def shifted_fake_data():
    q = np.linspace(0.015, 0.215, 100)  # Slightly shifted q range
    r = np.linspace(0.5, 0.5, 100)  # Different constant reflectivity
    err = 0.05 * np.ones_like(r)  # Constant error

    return ReducedData(q=q, r=r, err=err)


@pytest.fixture
def fake_data_no_overlap():
    q = np.linspace(0.3, 0.5, 100)  # No overlap with fake_data
    r = np.linspace(1, 1, 100)
    err = 0.05 * np.ones_like(r)  # Constant error

    return ReducedData(q=q, r=r, err=err)


@pytest.fixture
def overlap_scaling_factor(fake_data, shifted_fake_data):
    return OverlapScalingFactor(left_data=fake_data, right_data=shifted_fake_data, sf_auto=0.5)


### Test OverlapScalingFactor class methods ###


def test_apply_sf(overlap_scaling_factor):
    left_data_scaled = overlap_scaling_factor.apply_sf(overlap_scaling_factor.left_data)
    assert np.allclose(left_data_scaled.temp_r, overlap_scaling_factor.left_data.r * overlap_scaling_factor.sf_auto)


def test_get_fitting_overlap_range(overlap_scaling_factor):
    min_x = 0.05
    max_x = 0.1
    nbr_points = 5
    fit_range = overlap_scaling_factor.get_fitting_overlap_range(min_x, max_x, nbr_points)
    expected_range = np.array([0.05, 0.06, 0.07, 0.08, 0.09, 0.1])
    assert np.allclose(fit_range, expected_range)


def test_calculate_axis_overlap(overlap_scaling_factor):
    left_x_axis = overlap_scaling_factor.left_data.q
    right_x_axis = overlap_scaling_factor.right_data.q
    overlap = overlap_scaling_factor.calculate_axis_overlap(left_x_axis, right_x_axis)
    assert not overlap.no_overlap
    assert overlap.min_x < overlap.max_x


def test_fit_data(overlap_scaling_factor):
    left_set = overlap_scaling_factor.apply_sf(overlap_scaling_factor.left_data)
    [a_left, b_left] = overlap_scaling_factor.fit_data(left_set, 10, data_type="left")
    assert isinstance(a_left, float)
    assert isinstance(b_left, float)
    assert a_left == approx(0.0, abs=1e-9)  # Since data is constant
    assert b_left == approx(left_set.temp_r[10], abs=1e-9)

    right_set = overlap_scaling_factor.right_data
    [a_right, b_right] = overlap_scaling_factor.fit_data(right_set, 10, data_type="right")
    assert isinstance(a_right, float)
    assert isinstance(b_right, float)
    assert a_right == approx(0.0, abs=1e-9)  # Since data is constant
    assert b_right == approx(right_set.r[10], abs=1e-9)


def test_find_nearest(overlap_scaling_factor):
    axis = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    value = 0.33
    index = overlap_scaling_factor.find_nearest(axis, value)
    assert index == 2  # Nearest to 0.33 is 0.3 at index 2


def test_scale_factor_for_overlap_region(overlap_scaling_factor):
    fit_range = np.array([0.05, 0.06, 0.07, 0.08, 0.09, 0.1])
    a_left = 0.0
    b_left = 2.0
    a_right = 0.0
    b_right = 0.5
    sf = overlap_scaling_factor.scale_factor_for_overlap_region(fit_range, a_left, b_left, a_right, b_right)
    assert sf == approx(4.0, abs=1e-9)  # Since left mean is 2.0 and right mean is 0.5


def test_calculate_mean_over_range(overlap_scaling_factor):
    fit_range = np.array([0.05, 0.06, 0.07, 0.08, 0.09, 0.1])
    a = 0.0
    b = 2.0
    mean_value = overlap_scaling_factor.calculate_mean_over_range(fit_range, a, b)
    assert mean_value == approx(2.0, abs=1e-9)  # Since data is constant at 2.0
