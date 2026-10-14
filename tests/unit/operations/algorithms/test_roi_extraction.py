"""Unit tests for the ExtractDetectorROI Mantid algorithm."""

import os

import numpy as np
import pytest
from mantid.api import AlgorithmManager, EventType, IEventWorkspace
from mantid.kernel import DateAndTime
from mantid.simpleapi import CreateSampleWorkspace, DeleteWorkspace, LoadEventNexus, SetInstrumentParameter, mtd

from lr_reduction.operations.algorithms.roi_extraction import ExtractDetectorROI, roi_workspace_indices

# A non-square detector, so that swapping x and y in the index rule cannot go unnoticed.
N_X = 3
N_Y = 4

OUTPUT_WORKSPACE = "roi_extracted"


@pytest.fixture(autouse=True)
def clean_output():
    yield
    if mtd.doesExist(OUTPUT_WORKSPACE):
        DeleteWorkspace(OUTPUT_WORKSPACE)


@pytest.fixture
def make_detector():
    """Factory for an ``n_x`` x ``n_y`` event workspace with a single-bin TOF axis.

    Pixel ``index`` holds ``index`` events at distinct TOFs, so every spectrum is distinguishable.
    With ``weighted``, the events carry distinct non-unit weights and errors, as they do after a
    correction such as the Cd attenuation. The instrument's pixel counts are ``instrument_shape``,
    which defaults to the true shape; a dimension of None is left undefined. Deletes each workspace
    on teardown.
    """
    names = []

    def make(n_x=N_X, n_y=N_Y, instrument_shape=None, weighted=False):
        name = mtd.unique_hidden_name()
        CreateSampleWorkspace(
            WorkspaceType="Event",
            Function="Flat background",
            NumBanks=n_x * n_y,
            BankPixelWidth=1,
            NumEvents=0,
            XMin=0.0,
            XMax=1000.0,
            BinWidth=1000.0,
            OutputWorkspace=name,
        )
        workspace = mtd[name]
        for index in range(workspace.getNumberHistograms()):
            event_list = workspace.getSpectrum(index)
            if weighted:
                event_list.switchTo(EventType.WEIGHTED)
            for event in range(index):
                tof = 10.0 * event + index / 100.0
                if weighted:
                    event_list.addWeightedEventQuickly(tof, 1.0 + tof / 1000.0, (tof / 1000.0) ** 2, DateAndTime(0))
                else:
                    event_list.addEventQuickly(tof, DateAndTime(0))
        for parameter, value in zip(("number-of-x-pixels", "number-of-y-pixels"), instrument_shape or (n_x, n_y)):
            if value is not None:
                SetInstrumentParameter(
                    Workspace=name, ParameterName=parameter, ParameterType="Number", Value=f"{float(value)}"
                )
        names.append(name)
        return name

    yield make
    for name in names:
        if mtd.doesExist(name):
            DeleteWorkspace(name)


def run_algorithm(input_workspace, x_range, y_range, output_workspace=OUTPUT_WORKSPACE):
    """Run the algorithm as a managed, top-level algorithm and return the output workspace."""
    algorithm = AlgorithmManager.create("ExtractDetectorROI")
    algorithm.setProperty("InputWorkspace", input_workspace)
    algorithm.setProperty("XRange", x_range)
    algorithm.setProperty("YRange", y_range)
    algorithm.setProperty("OutputWorkspace", output_workspace)
    algorithm.execute()
    return mtd[output_workspace]


def assert_same_pixel(output_events, input_events):
    """The two spectra are the same detector pixel with the same events, weights and errors."""
    assert set(output_events.getDetectorIDs()) == set(input_events.getDetectorIDs())
    assert output_events.getEventType() == input_events.getEventType()
    assert output_events.getNumberEvents() == input_events.getNumberEvents()
    np.testing.assert_array_equal(output_events.getTofs(), input_events.getTofs())
    np.testing.assert_array_equal(output_events.getWeights(), input_events.getWeights())
    np.testing.assert_array_equal(output_events.getWeightErrors(), input_events.getWeightErrors())


def test_algorithm_metadata():
    algorithm = ExtractDetectorROI()

    assert algorithm.name() == "ExtractDetectorROI"
    assert algorithm.category() == "Reflectometry\\SNS"
    assert algorithm.version() == 1


def test_index_rule_matches_the_legacy_convention():
    """Legacy BL4B code reads pixel (x, y) from spectrum ``x * n_y + y``, y running fastest."""
    assert roi_workspace_indices([1, 2], [0, 1], n_y=N_Y) == [4, 5, 8, 9]
    assert roi_workspace_indices([0, 0], [3, 3], n_y=304) == [3]
    assert roi_workspace_indices([255, 255], [303, 303], n_y=304) == [256 * 304 - 1]


@pytest.mark.parametrize("weighted", [False, True], ids=["unweighted", "weighted"])
def test_output_holds_exactly_the_pixels_inside_the_rectangle(make_detector, weighted):
    """x in [1, 2], y in [1, 3] on the 3 x 4 detector is spectra 5-7 and 9-11, in that order."""
    input_workspace = make_detector(weighted=weighted)
    if weighted:
        assert mtd[input_workspace].getSpectrum(5).getEventType() == EventType.WEIGHTED
        assert not np.allclose(mtd[input_workspace].getSpectrum(5).getWeights(), 1.0)

    output = run_algorithm(input_workspace, [1, 2], [1, 3])

    expected = [5, 6, 7, 9, 10, 11]
    assert output.getNumberHistograms() == len(expected)
    for output_index, input_index in enumerate(expected):
        assert_same_pixel(output.getSpectrum(output_index), mtd[input_workspace].getSpectrum(input_index))


def test_single_pixel_roi(make_detector):
    """Pixel (2, 3) is the last spectrum, 2 * 4 + 3 = 11."""
    input_workspace = make_detector()

    output = run_algorithm(input_workspace, [2, 2], [3, 3])

    assert output.getNumberHistograms() == 1
    assert_same_pixel(output.getSpectrum(0), mtd[input_workspace].getSpectrum(11))


def test_full_detector_keeps_every_pixel(make_detector):
    input_workspace = make_detector()

    output = run_algorithm(input_workspace, [0, N_X - 1], [0, N_Y - 1])

    assert output.getNumberHistograms() == N_X * N_Y
    assert output.getNumberEvents() == mtd[input_workspace].getNumberEvents()


def test_input_is_unchanged(make_detector):
    input_workspace = make_detector()
    before = [mtd[input_workspace].getSpectrum(i).getTofs().copy() for i in range(N_X * N_Y)]

    run_algorithm(input_workspace, [0, 1], [2, 3])

    workspace = mtd[input_workspace]
    assert workspace.getNumberHistograms() == N_X * N_Y
    for index, tofs in enumerate(before):
        np.testing.assert_array_equal(workspace.getSpectrum(index).getTofs(), tofs)


def test_output_can_replace_the_input(make_detector):
    input_workspace = make_detector()
    expected = [mtd[input_workspace].getSpectrum(index).getTofs().copy() for index in (1, 2, 5, 6)]

    output = run_algorithm(input_workspace, [0, 1], [1, 2], output_workspace=input_workspace)

    assert mtd[input_workspace].getNumberHistograms() == 4
    for output_index, tofs in enumerate(expected):
        np.testing.assert_array_equal(output.getSpectrum(output_index).getTofs(), tofs)


def test_output_stays_an_event_workspace_with_the_placeholder_axis(make_detector):
    """ExtractSpectra must not histogram: the output keeps its events and the input's single bin."""
    input_workspace = make_detector()

    output = run_algorithm(input_workspace, [0, 1], [0, 1])

    assert isinstance(output, IEventWorkspace)
    assert output.blocksize() == 1
    for index in range(output.getNumberHistograms()):
        np.testing.assert_array_equal(output.readX(index), mtd[input_workspace].readX(0))


@pytest.mark.parametrize(
    ("x_range", "y_range", "match"),
    [
        ([0, 3], [0, 1], r"XRange \[0, 3\] runs off the 3 x 4 detector, whose x pixels are 0 to 2"),
        ([0, 1], [2, 4], r"YRange \[2, 4\] runs off the 3 x 4 detector, whose y pixels are 0 to 3"),
    ],
)
def test_invalid_range_raises(make_detector, x_range, y_range, match):
    with pytest.raises(RuntimeError, match=match):
        run_algorithm(make_detector(), x_range, y_range)


def test_every_invalid_range_is_reported(make_detector):
    with pytest.raises(RuntimeError) as error:
        run_algorithm(make_detector(), [0, 3], [2, 4])

    assert "XRange [0, 3] runs off the 3 x 4 detector" in str(error.value)
    assert "YRange [2, 4] runs off the 3 x 4 detector" in str(error.value)


def test_instrument_size_disagreeing_with_the_spectrum_count_raises(make_detector):
    input_workspace = make_detector(instrument_shape=(3, 3))

    with pytest.raises(RuntimeError, match="3 x 3 detector has 9 pixels, but the workspace has 12 spectra"):
        run_algorithm(input_workspace, [0, 1], [0, 1])


@pytest.mark.parametrize(
    ("instrument_shape", "found"), [((None, 4), "None x 4"), ((3, None), "3 x None")], ids=["no_x", "no_y"]
)
def test_instrument_without_pixel_counts_raises(make_detector, instrument_shape, found):
    input_workspace = make_detector(instrument_shape=instrument_shape)

    with pytest.raises(RuntimeError, match=rf"must define number-of-x-pixels and number-of-y-pixels \(found {found}\)"):
        run_algorithm(input_workspace, [0, 1], [0, 1])


@pytest.mark.parametrize(
    ("value", "match"),
    [
        ([0], "Incorrect size"),
        ([0, 1, 2], "Incorrect size"),
        ([-1, 1], r"Selected value -1 is < the lower bound \(0\)"),
        ([2, 1], r"Pair \(2, 1\) is not ordered"),
    ],
    ids=["too_short", "too_long", "negative", "reversed"],
)
@pytest.mark.parametrize("name", ["XRange", "YRange"])
def test_property_validators_reject_malformed_ranges(name, value, match):
    """Each malformed range is caught by its own validator: length, lower bound and ordering."""
    algorithm = ExtractDetectorROI()
    algorithm.PyInit()

    with pytest.raises(ValueError, match=match):
        algorithm.setProperty(name, value)


def test_history_records_the_ranges(make_detector):
    output = run_algorithm(make_detector(), [0, 2], [1, 3])

    history = output.getHistory().lastAlgorithm()
    reference = ExtractDetectorROI()
    reference.PyInit()
    reference.setProperty("XRange", [0, 2])
    reference.setProperty("YRange", [1, 3])
    assert history.name() == "ExtractDetectorROI"
    assert history.getPropertyValue("XRange") == reference.getPropertyValue("XRange")
    assert history.getPropertyValue("YRange") == reference.getPropertyValue("YRange")


@pytest.fixture(scope="module")
def measured_run(nexus_dir):
    """A REF_L run, loaded once for this module."""
    name = "roi_measured_run"
    LoadEventNexus(Filename=os.path.join(nexus_dir, "REF_L_198388.nxs.h5"), OutputWorkspace=name)
    yield name
    DeleteWorkspace(name)


@pytest.mark.datarepo
def test_roi_on_a_measured_run(measured_run):
    """The REF_L instrument definition gives the 256 x 304 detector, read with the legacy index rule.

    The ROI covers the run's beam spot, so it holds events. On REF_L, each spectrum's detector ID
    equals its workspace index, ``x * 304 + y``.
    """
    input_workspace = mtd[measured_run]
    assert input_workspace.getNumberHistograms() == 256 * 304

    output = run_algorithm(measured_run, [100, 159], [130, 159])

    expected = [x * 304 + y for x in range(100, 160) for y in range(130, 160)]
    assert isinstance(output, IEventWorkspace)
    assert output.getNumberHistograms() == len(expected)
    assert output.getNumberEvents() > 0
    for output_index, input_index in enumerate(expected):
        assert set(output.getSpectrum(output_index).getDetectorIDs()) == {input_index}
        assert_same_pixel(output.getSpectrum(output_index), input_workspace.getSpectrum(input_index))
