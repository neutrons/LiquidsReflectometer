import numpy as np
import pytest

from lr_reduction.exceptions import IncompleteDataError, MalformedDataError
from lr_reduction.models.results import ReflectivityCurve


def _four_columns():
    return {
        "q": np.array([0.01, 0.02]),
        "r": np.array([1.0, 0.5]),
        "dr": np.array([0.1, 0.05]),
        "dq": np.array([0.001, 0.002]),
    }


def test_to_array_four_column_order():
    curve = ReflectivityCurve(**_four_columns())
    array = curve.to_array()
    assert array.shape == (2, 4)
    np.testing.assert_array_equal(array[:, 0], curve.q)
    np.testing.assert_array_equal(array[:, 1], curve.r)
    np.testing.assert_array_equal(array[:, 2], curve.dr)
    np.testing.assert_array_equal(array[:, 3], curve.dq)
    assert not curve.is_eight_column


def test_to_array_eight_column_order():
    curve = ReflectivityCurve(
        **_four_columns(),
        theta=np.array([0.6, 0.6]),
        lam=np.array([4.0, 5.0]),
        dtheta=np.array([0.01, 0.01]),
        dlam=np.array([0.04, 0.05]),
    )
    array = curve.to_array()
    assert array.shape == (2, 8)
    np.testing.assert_array_equal(array[:, 4], curve.theta)
    np.testing.assert_array_equal(array[:, 5], curve.lam)
    np.testing.assert_array_equal(array[:, 6], curve.dtheta)
    np.testing.assert_array_equal(array[:, 7], curve.dlam)
    assert curve.is_eight_column


def test_partial_eight_column_rejected():
    with pytest.raises(IncompleteDataError, match="provided together"):
        ReflectivityCurve(**_four_columns(), theta=np.array([0.6, 0.6]))


def test_multidimensional_column_rejected():
    columns = _four_columns()
    columns["r"] = np.ones((2, 2))
    with pytest.raises(MalformedDataError, match="one-dimensional"):
        ReflectivityCurve(**columns)


def test_mismatched_column_length_rejected():
    columns = _four_columns()
    columns["dq"] = np.array([0.001, 0.002, 0.003])
    with pytest.raises(MalformedDataError, match="equal lengths"):
        ReflectivityCurve(**columns)


def test_empty_curve():
    curve = ReflectivityCurve.empty()
    assert curve.to_array().shape == (0, 4)
