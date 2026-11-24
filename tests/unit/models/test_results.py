import numpy as np
import pytest

from lr_reduction.exceptions import IncompleteDataError, MalformedDataError
from lr_reduction.models.results import CompositeDirectBeam, DirectBeamSegment, ReflectivityCurve


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


def _segment(run_number, cd_thickness_um=0.0):
    return DirectBeamSegment(
        run_number=run_number,
        cd_thickness_um=cd_thickness_um,
        slit_scan_factor=1.0,
        dead_time_threshold=1.25,
        tof_range=(12000.0, 45000.0),
        chopper_cut=11000.0,
        dead_time_cut=12000.0,
        intensity_cut=45000.0,
        wavelength_band=(2.5, 9.5),
        db_pixel=141.3,
        db_pixel_err=0.02,
        tthd=0.0,
    )


@pytest.mark.parametrize(
    "run_numbers",
    [(12346, 12345), (12345,), (12345, 12346, 12347)],
    ids=["reordered", "missing segment", "extra run"],
)
def test_composite_direct_beam_rejects_segments_not_matching_run_numbers(run_numbers):
    with pytest.raises(MalformedDataError, match="segments are for runs"):
        CompositeDirectBeam(
            workspace="REF_L_12345_db",
            run_numbers=run_numbers,
            segments=(_segment(12345), _segment(12346)),
            gaps=(),
            weights_workspace=None,
        )
