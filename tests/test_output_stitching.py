"""Unit tests for RunCollection stitching functionality"""
import numpy as np
import pytest
from pytest import approx

from lr_reduction.output import RunCollection
from lr_reduction.reduction_template_reader import ReductionParameters
from lr_reduction.scaling_factors.calculate import StitchingType


@pytest.fixture
def template_data_none():
    """Create template data with no stitching"""
    template = ReductionParameters()
    template.stitching_type = StitchingType.NONE
    return template


@pytest.fixture
def template_data_automatic():
    """Create template data with automatic stitching"""
    template = ReductionParameters()
    template.stitching_type = StitchingType.AUTOMATIC_AVERAGE
    template.sf_qmin = 0.01
    template.sf_qmax = 0.05
    return template


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

    def test_init_with_template_data(self, template_data_none):
        """Test RunCollection initialization with template data"""
        coll = RunCollection(template_data=template_data_none)
        assert coll.template_data == template_data_none
        assert coll.collection == []

    def test_calculate_scale_factors_none(self, template_data_none, mock_run_data):
        """Test that NONE stitching type returns all 1.0 scale factors"""
        coll = RunCollection(template_data=template_data_none)

        for q, r, dr, meta in mock_run_data:
            coll.add(q, r, dr, meta_data=meta)

        coll.stitching_type = template_data_none.stitching_type
        coll.scale_factors = []
        coll.calculate_scale_factors()

        assert len(coll.scale_factors) == 2
        assert all(sf == 1.0 for sf in coll.scale_factors)

    def test_calculate_scale_factors_automatic(self, template_data_automatic, mock_run_data):
        """Test automatic stitching scale factor calculation"""
        coll = RunCollection(template_data=template_data_automatic)

        for q, r, dr, meta in mock_run_data:
            coll.add(q, r, dr, meta_data=meta)

        coll.stitching_type = template_data_automatic.stitching_type
        coll.calculate_scale_factors()

        assert len(coll.scale_factors) == 2
        # First scale factor should be from critical edge
        assert isinstance(coll.scale_factors[0], float)
        # Second scale factor should be from overlap
        assert isinstance(coll.scale_factors[1], float)

    def test_merge_applies_scale_factors(self, template_data_none, mock_run_data):
        """Test that merge applies scale factors to reflectivity"""
        coll = RunCollection(template_data=template_data_none)

        q1, r1, dr1, meta1 = mock_run_data[0]
        coll.add(q1, r1, dr1, meta_data=meta1)

        # Set manual scale factor
        coll.stitching_type = StitchingType.NONE
        coll.scale_factors = [2.0]  # Scale by 2
        coll.merge()

        # Check that reflectivity was scaled
        assert coll.refl_all[0] == approx(r1[0] * 2.0)
        assert coll.d_refl_all[0] == approx(dr1[0] * 2.0)



if __name__ == "__main__":
    pytest.main([__file__])
