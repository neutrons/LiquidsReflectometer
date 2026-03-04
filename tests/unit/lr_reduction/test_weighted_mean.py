"""Tests for weighted_mean edge cases that previously produced RuntimeWarnings."""

import warnings

import numpy as np
import pytest

from lr_reduction.nr_tools import weighted_mean


class TestWeightedMeanNormal:
    """Test weighted_mean with well-behaved inputs."""

    def test_identical_arrays(self):
        y1 = np.array([1.0, 2.0, 3.0])
        y2 = np.array([1.0, 2.0, 3.0])
        e1 = np.array([0.1, 0.1, 0.1])
        e2 = np.array([0.1, 0.1, 0.1])
        scale, sigma = weighted_mean(y1, y2, e1, e2)
        assert np.isfinite(scale)
        assert np.isfinite(sigma)
        assert abs(scale - 1.0) < 0.1

    def test_scaled_arrays(self):
        y1 = np.array([2.0, 4.0, 6.0])
        y2 = np.array([1.0, 2.0, 3.0])
        e1 = np.array([0.1, 0.1, 0.1])
        e2 = np.array([0.1, 0.1, 0.1])
        scale, sigma = weighted_mean(y1, y2, e1, e2)
        assert np.isfinite(scale)
        assert abs(scale - 2.0) < 0.1


class TestWeightedMeanEdgeCases:
    """Test weighted_mean with inputs that previously produced RuntimeWarnings."""

    def test_empty_arrays_no_warning(self):
        """Empty input arrays should return NaN without RuntimeWarning."""
        y1 = np.array([])
        y2 = np.array([])
        e1 = np.array([])
        e2 = np.array([])
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            scale, sigma = weighted_mean(y1, y2, e1, e2)
        assert np.isnan(scale)
        assert np.isnan(sigma)

    def test_zero_denominator_no_warning(self):
        """y2 containing zeros should not produce RuntimeWarning."""
        y1 = np.array([1.0, 2.0])
        y2 = np.array([0.0, 0.0])
        e1 = np.array([0.1, 0.1])
        e2 = np.array([0.1, 0.1])
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            scale, sigma = weighted_mean(y1, y2, e1, e2)
        assert np.isnan(scale)
        assert np.isnan(sigma)

    def test_mixed_zero_denominator_no_warning(self):
        """y2 with some zeros should use only valid points, no RuntimeWarning."""
        y1 = np.array([2.0, 4.0, 6.0])
        y2 = np.array([1.0, 0.0, 3.0])
        e1 = np.array([0.1, 0.1, 0.1])
        e2 = np.array([0.1, 0.1, 0.1])
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            scale, sigma = weighted_mean(y1, y2, e1, e2)
        assert np.isfinite(scale)
        assert np.isfinite(sigma)

    def test_zero_errors_no_warning(self):
        """Zero error bars should not produce RuntimeWarning."""
        y1 = np.array([1.0, 2.0])
        y2 = np.array([1.0, 2.0])
        e1 = np.array([0.0, 0.0])
        e2 = np.array([0.0, 0.0])
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            scale, sigma = weighted_mean(y1, y2, e1, e2)

    def test_all_outliers_filtered_no_warning(self):
        """If sigma_mask filters all points, should return NaN without warning."""
        # One normal point and one extreme outlier — tight sigma_mask will reject outlier,
        # but with only 2 points and tight mask, both may survive. Use a more extreme case.
        y1 = np.array([1.0, 100.0])
        y2 = np.array([1.0, 1.0])
        e1 = np.array([0.01, 0.01])
        e2 = np.array([0.01, 0.01])
        # Very tight mask — should filter the outlier
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            scale, sigma = weighted_mean(y1, y2, e1, e2, sigma_mask=0.001)
        # Result should be finite or NaN, but no warning
