"""Unit tests for RunCollection stitching functionality"""
import numpy as np
import pytest
from pytest import approx

from lr_reduction.output import RunCollection
from lr_reduction.scaling_factors.calculate import StitchingConfiguration, StitchingType


@pytest.fixture
def stitching_configuration_default():
    """Create default stitching configuration"""
    return StitchingConfiguration()


@pytest.fixture
def stitching_configuration_automatic_average():
    """Create a stitching configuration with automatic stitching"""
    stitching_configuration = StitchingConfiguration()
    stitching_configuration.type = StitchingType.AUTOMATIC_AVERAGE
    stitching_configuration.scale_factor_qmin = 0.01
    stitching_configuration.scale_factor_qmax = 0.05
    stitching_configuration.normalize_first_angle = False
    return stitching_configuration


@pytest.fixture
def mock_run_data():
    """Create mock run data for testing"""
    # First run
    q1 = np.linspace(0.01, 0.1, 50)
    r1 = np.ones_like(q1) * 2.0
    dr1 = np.ones_like(q1) * 0.1
    meta1 = {
        "experiment": "TEST",
        "run_number": 1001,
        "sequence_id": 1,
        "scaling_factors": {"a": 1.0, "b": 0.0},
        "two_theta": 1.0,
        "lambda_min": 2.0,
        "lambda_max": 10.0,
        "q_min": 0.01,
        "q_max": 0.1,
        "dq_over_q": 0.027,
    }

    # Second run with overlap
    q2 = np.linspace(0.08, 0.2, 50)
    r2 = np.ones_like(q2) * 0.5
    dr2 = np.ones_like(q2) * 0.05
    meta2 = {
        "experiment": "TEST",
        "run_number": 1002,
        "sequence_id": 2,
        "scaling_factors": {"a": 1.0, "b": 0.0},
        "two_theta": 2.0,
        "lambda_min": 2.0,
        "lambda_max": 10.0,
        "q_min": 0.08,
        "q_max": 0.2,
        "dq_over_q": 0.027,
    }

    return [(q1, r1, dr1, meta1), (q2, r2, dr2, meta2)]


class TestRunCollectionStitching:
    """Test RunCollection stitching functionality"""

    def test_calculate_scale_factors_none(self, stitching_configuration_default, mock_run_data):
        """Test that NONE stitching type returns all 1.0 scale factors"""
        coll = RunCollection(stitching_configuration=stitching_configuration_default)

        for q, r, dr, meta in mock_run_data:
            coll.add(q, r, dr, meta_data=meta)

        coll.calculate_scale_factors()

        assert len(coll.stitching_reflectivity_scale_factors) == 2
        assert all(sf == 1.0 for sf in coll.stitching_reflectivity_scale_factors)

    def test_calculate_scale_factors_automatic(self, stitching_configuration_automatic_average, mock_run_data):
        """Test automatic stitching scale factor calculation"""
        coll = RunCollection(stitching_configuration=stitching_configuration_automatic_average)

        for q, r, dr, meta in mock_run_data:
            coll.add(q, r, dr, meta_data=meta)

        coll.calculate_scale_factors()

        assert len(coll.stitching_reflectivity_scale_factors) == 2
        assert coll.stitching_reflectivity_scale_factors[0] == approx(1.0)
        assert coll.stitching_reflectivity_scale_factors[1] == approx(4.0)

    def test_calculate_scale_factors_automatic_out_of_order(
        self, stitching_configuration_automatic_average, mock_run_data
    ):
        """Test automatic stitching when runs are added out of q-order"""
        coll = RunCollection(stitching_configuration=stitching_configuration_automatic_average)
        # Add higher-q run first, then lower-q run
        q2, r2, dr2, meta2 = mock_run_data[1]
        q1, r1, dr1, meta1 = mock_run_data[0]
        coll.add(q2, r2, dr2, meta_data=meta2)
        coll.add(q1, r1, dr1, meta_data=meta1)
        coll.calculate_scale_factors()
        # Scale factors should be reversed from the base case (test_calculate_scale_factors_automatic)
        assert len(coll.stitching_reflectivity_scale_factors) == 2
        assert coll.stitching_reflectivity_scale_factors[0] == approx(4.0)
        assert coll.stitching_reflectivity_scale_factors[1] == approx(1.0)

    def test_merge_applies_scale_factors(self, stitching_configuration_default, mock_run_data):
        """Test that merge applies scale factors to reflectivity"""
        coll = RunCollection(stitching_configuration=stitching_configuration_default)

        q1, r1, dr1, meta1 = mock_run_data[0]
        coll.add(q1, r1, dr1, meta_data=meta1)

        # Set manual scale factor
        coll.stitching_reflectivity_scale_factors = [2.0]  # Scale by 2
        coll.merge()

        # Check that reflectivity was scaled
        assert coll.refl_all[0] == approx(r1[0] * 2.0)
        assert coll.d_refl_all[0] == approx(dr1[0] * 2.0)

    def test_normalize_first_angle_true(self, stitching_configuration_automatic_average, mock_run_data):
        """Test that setting normalize_first_angle True updates the first edge"""
        stitching_configuration_automatic_average.normalize_first_angle = True
        coll = RunCollection(stitching_configuration=stitching_configuration_automatic_average)

        for q, r, dr, meta in mock_run_data:
            coll.add(q, r, dr, meta_data=meta)

        coll.calculate_scale_factors()

        assert len(coll.stitching_reflectivity_scale_factors) == 2
        assert coll.stitching_reflectivity_scale_factors[0] == approx(0.5)
        assert coll.stitching_reflectivity_scale_factors[1] == approx(2.0)


if __name__ == "__main__":
    pytest.main([__file__])
