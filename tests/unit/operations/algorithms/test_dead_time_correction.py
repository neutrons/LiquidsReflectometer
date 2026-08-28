"""Unit tests for the SingleReadoutDeadTimeCorrection Mantid algorithm."""

import os

import mantid.simpleapi as mtd_api
import numpy as np
import pytest

from lr_reduction.operations.algorithms.dead_time_correction import SingleReadoutDeadTimeCorrection


@pytest.fixture(scope="module")
def event_workspace(nexus_dir):
    """Event workspace for a single REF_L run, loaded once for all tests in this module."""
    return mtd_api.LoadEventNexus(
        Filename=os.path.join(nexus_dir, "REF_L_198409.nxs.h5"), OutputWorkspace="dead_time_input"
    )


def run_algorithm(input_workspace, **properties):
    """Run the algorithm on ``input_workspace`` and return the output workspace."""
    algorithm = SingleReadoutDeadTimeCorrection()
    algorithm.PyInit()
    algorithm.setProperty("InputWorkspace", input_workspace)
    algorithm.setProperty("OutputWorkspace", "dead_time_corr")
    for name, value in properties.items():
        algorithm.setProperty(name, value)
    algorithm.PyExec()
    return algorithm.getProperty("OutputWorkspace").value


def compute_correction(input_workspace, **properties):
    """Run the algorithm on ``input_workspace`` and return the correction for each TOF bin."""
    return run_algorithm(input_workspace, **properties).readY(0)


@pytest.mark.datarepo
def test_deadtime(event_workspace):
    """The non-paralyzable correction is a small positive scaling for every TOF bin."""
    correction = compute_correction(event_workspace)

    for value in correction:
        assert value > 0
        assert value < 1.001


@pytest.mark.datarepo
def test_deadtime_paralyzable(event_workspace):
    """The paralyzable correction is a small positive scaling for every TOF bin."""
    correction = compute_correction(event_workspace, Paralyzable=True)

    for value in correction:
        assert value > 0
        assert value < 1.001


@pytest.mark.datarepo
def test_deadtime_threshold(event_workspace):
    """TOF bins requiring a correction above the threshold are zeroed out."""
    # the mean correction is 1.0003, so use that as the test value
    correction = compute_correction(event_workspace, UseDeadTimeThreshold=True, DeadTimeThreshold=1.0003)

    # manual inspection showed 88 corrections were above the threshold
    assert len(correction[correction == 0]) == 88

    for value in correction:
        assert value <= 1.0003


# Correction values recorded from REF_L_198409 at commit 58668e6.
RTOL = 1e-12

CORRECTION_STEP_2000 = [
    1.0004051852254192,
    1.0007312476124917,
    1.0005003675582942,
    1.000421458641246,
    1.0003282635169768,
    1.0002528325613604,
    1.0001749484805473,
    1.0001148184375086,
    1.0000305499712825,
]

PARALYZABLE_CORRECTION_STEP_2000 = [
    1.0004052673573192,
    1.0007315152349066,
    1.0005004928257026,
    1.00042154750487,
    1.000328317419035,
    1.00025286453429,
    1.000174963787603,
    1.0001148250301546,
    1.0000305504379519,
]

CORRECTION_TOF_RANGE_55000_60000 = [
    1.0004880372194696,
    1.000454500235419,
    1.000378557141473,
    1.000373626170375,
    1.0003124861670771,
]

# Every 20th bin of the correction computed with the default TOFStep of 100 microseconds, sampled
# across the full TOF range so the production default is pinned without listing all 167 bins.
CORRECTION_DEFAULT_EVERY_20TH_BIN = [
    1.0000098546255032,
    1.0008185954067546,
    1.0005521583004042,
    1.0005521583004042,
    1.0004041989748087,
    1.0003548889235223,
    1.0001774129808103,
    1.0001379824339913,
    1.0,
]


@pytest.mark.datarepo
def test_correction_values_are_stable(event_workspace):
    """The non-paralyzable correction reproduces its recorded values bin for bin."""
    correction = compute_correction(event_workspace, TOFStep=2000.0)

    np.testing.assert_allclose(correction, CORRECTION_STEP_2000, rtol=RTOL, atol=0)


@pytest.mark.datarepo
def test_paralyzable_correction_values_are_stable(event_workspace):
    """The paralyzable correction reproduces its recorded values bin for bin."""
    correction = compute_correction(event_workspace, TOFStep=2000.0, Paralyzable=True)

    np.testing.assert_allclose(correction, PARALYZABLE_CORRECTION_STEP_2000, rtol=RTOL, atol=0)


@pytest.mark.datarepo
def test_correction_values_with_explicit_tof_range(event_workspace):
    """An explicit TOFRange bins the correction over that range rather than the workspace range."""
    output_workspace = run_algorithm(event_workspace, TOFRange=[55000.0, 60000.0], TOFStep=1000.0)

    np.testing.assert_allclose(output_workspace.readY(0), CORRECTION_TOF_RANGE_55000_60000, rtol=RTOL, atol=0)
    np.testing.assert_allclose(output_workspace.readX(0), np.arange(55000.0, 60001.0, 1000.0), rtol=RTOL, atol=0)


@pytest.mark.datarepo
def test_default_binning_spans_the_workspace_tof_range(event_workspace):
    """Without a TOFRange, the correction is binned in 100 microsecond steps over the run's TOF range."""
    output_workspace = run_algorithm(event_workspace)
    bin_edges = output_workspace.readX(0)

    assert bin_edges[0] == event_workspace.getTofMin()
    assert bin_edges[-1] == event_workspace.getTofMax()
    np.testing.assert_allclose(np.diff(bin_edges[:-1]), 100.0, rtol=RTOL, atol=0)

    correction = output_workspace.readY(0)
    assert len(correction) == 167
    np.testing.assert_allclose(correction[::20], CORRECTION_DEFAULT_EVERY_20TH_BIN, rtol=RTOL, atol=0)


@pytest.mark.datarepo
def test_output_workspace_is_a_distribution_without_errors(event_workspace):
    """The correction is a single distribution spectrum carrying no uncertainty of its own."""
    output_workspace = run_algorithm(event_workspace)

    assert output_workspace.getNumberHistograms() == 1
    assert output_workspace.isDistribution()
    np.testing.assert_array_equal(output_workspace.readE(0), 0.0)
