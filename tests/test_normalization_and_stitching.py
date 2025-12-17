import numpy as np
import pytest
from numpy.ma.testutils import approx

from lr_reduction.normalization_and_stitching import OverlapScalingFactor, ReducedData, scaling_factor_critical_edge


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

def test_scaling_factor_critical_edge(fake_data, shifted_fake_data):
    # Define critical edge range
    q_min = 0.05
    q_max = 0.1

    # Calculate scaling factor for fake_data
    sf_linear = scaling_factor_critical_edge(q_min, q_max, [fake_data])
    assert isinstance(sf_linear, float)
    assert approx(sf_linear, 0.5, rtol=1e-1)

    # Calculate scaling factor for shifted_fake_data
    sf_shifted = scaling_factor_critical_edge(q_min, q_max, [shifted_fake_data])
    assert isinstance(sf_shifted, float)
    assert approx(sf_shifted, 2.0, rtol=1e-1)

    # Check that the scaling factors are different due to different data
    assert sf_linear != sf_shifted


def test_overlap_scaling_factor(fake_data, shifted_fake_data):
    # Create OverlapScalingFactor instance
    overlap_sf = OverlapScalingFactor(left_data=fake_data, right_data=shifted_fake_data)

    # Get scaling factor
    sf = overlap_sf.get_scaling_factor()

    assert isinstance(sf, float)
    assert approx(sf, 0.25, rtol=1e-1)  # Since fake_data r=2 and shifted_fake_data r=0.5

def test_overlap_scaling_factor_no_overlap(fake_data, fake_data_no_overlap):
    # Create OverlapScalingFactor instance with no overlap
    overlap_sf = OverlapScalingFactor(left_data=fake_data, right_data=fake_data_no_overlap)

    # Get scaling factor
    sf = overlap_sf.get_scaling_factor()

    assert isinstance(sf, float)
    assert approx(sf, 1.0, rtol=1e-5)  # Should return 1.0 when no overlap