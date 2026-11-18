"""Unit tests for the FitBeamPeak Mantid algorithm and the SuperGaussian fit function."""

import itertools
import os

import numpy as np
import pytest
from mantid.api import AlgorithmManager, EventType, FunctionFactory, IEventWorkspace
from mantid.fitfunctions import FunctionWrapper
from mantid.kernel import DateAndTime
from mantid.simpleapi import (
    ConvertUnits,
    CreateSampleWorkspace,
    DeleteWorkspace,
    ExtractSpectra,
    LoadEventNexus,
    Rebin,
    Scale,
    SetInstrumentParameter,
    mtd,
)

from lr_reduction.legacy.peak_finding import fit_signal_flat_bck
from lr_reduction.operations.algorithms.fit_beam_peak import CONVERGED, FitBeamPeak, SuperGaussian, y_profile

# A non-square detector, so that swapping x and y in the pixel rule cannot go unnoticed.
N_X = 2
N_Y = 40

# The column-summed truth: a super-Gaussian between two pixels on a sloped background.
CENTRE = 21.5
HEIGHT = 3000.0
SIGMA = 4.0
EXPONENT = 4.0

PROFILE_WORKSPACE = "beam_peak_profile"


def truth(y):
    return HEIGHT * np.exp(-(np.abs((y - CENTRE) / SIGMA) ** EXPONENT)) + 30.0 + 0.5 * y


@pytest.fixture(autouse=True)
def clean_output():
    yield
    if mtd.doesExist(PROFILE_WORKSPACE):
        DeleteWorkspace(PROFILE_WORKSPACE)


@pytest.fixture
def make_detector():
    """Factory for an ``N_X`` x ``n_y`` event workspace whose column sums follow `profile(y)`.

    The pixel of detector ID ``d`` is ``(d // n_y, d % n_y)``, the rule the algorithm reads; each
    pixel gets a Poisson draw of ``profile(y) / N_X`` events at random TOFs. With ``weighted``,
    every event has weight 2 and error² 4. A dimension of ``instrument_shape`` that is None is left
    undefined. Deletes each workspace on teardown.
    """
    names = []
    rng = np.random.default_rng(2026)

    def make(profile=truth, n_y=N_Y, weighted=False, instrument_shape=None):
        name = mtd.unique_hidden_name()
        CreateSampleWorkspace(
            WorkspaceType="Event",
            Function="Flat background",
            NumBanks=N_X * n_y,
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
            y = event_list.getDetectorIDs()[0] % n_y
            if weighted:
                event_list.switchTo(EventType.WEIGHTED)
            for tof in rng.uniform(1.0, 999.0, rng.poisson(profile(y) / N_X)):
                if weighted:
                    event_list.addWeightedEventQuickly(tof, 2.0, 4.0, DateAndTime(0))
                else:
                    event_list.addEventQuickly(tof, DateAndTime(0))
        for parameter, value in zip(("number-of-x-pixels", "number-of-y-pixels"), instrument_shape or (N_X, n_y)):
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


def run_algorithm(input_workspace, peak="supergauss", background="linear", output_workspace=None):
    """Run the algorithm as a managed, top-level algorithm and return it."""
    algorithm = AlgorithmManager.create("FitBeamPeak")
    algorithm.setRethrows(True)
    algorithm.setProperty("InputWorkspace", input_workspace)
    algorithm.setProperty("PeakFunction", peak)
    algorithm.setProperty("BackgroundFunction", background)
    if output_workspace is not None:
        algorithm.setProperty("OutputWorkspace", output_workspace)
    algorithm.execute()
    return algorithm


def centre(algorithm):
    return algorithm.getProperty("PeakCentre").value


def test_algorithm_metadata():
    algorithm = FitBeamPeak()

    assert algorithm.name() == "FitBeamPeak"
    assert algorithm.category() == "Reflectometry\\SNS"
    assert algorithm.version() == 1


@pytest.mark.parametrize("weight", [1.0, 1e-3], ids=["counts", "normalized"])
def test_y_profile_sums_over_x_in_quadrature_and_gives_empty_pixels_one_event_error(weight):
    """IDs 1, 5 and 9 are y = 1 on a 4-pixel-high detector; y = 2 has no counts.

    Every event weighs `weight`, so a pixel of n events has counts ``n w`` and error ``√n w``.
    """
    n_events = np.array([1.0, 2.0, 3.0, 0.0, 4.0])
    y_pixels, counts, errors = y_profile([1, 5, 9, 2, 3], weight * n_events, weight * np.sqrt(n_events), 4)

    np.testing.assert_array_equal(y_pixels, [1, 2, 3])
    np.testing.assert_allclose(counts, weight * np.array([6.0, 0.0, 4.0]))
    np.testing.assert_allclose(errors, weight * np.array([np.sqrt(6.0), 1.0, 2.0]))


@pytest.mark.parametrize("weighted", [False, True], ids=["unweighted", "weighted"])
def test_centroid_recovery(make_detector, weighted):
    algorithm = run_algorithm(make_detector(weighted=weighted))

    error = algorithm.getProperty("PeakCentreError").value
    assert 0 < error < 0.1
    assert abs(centre(algorithm) - CENTRE) < max(3 * error, 0.02)
    assert abs(centre(algorithm) - CENTRE) < 0.1
    assert 0.5 < algorithm.getProperty("ChiSquaredOverDoF").value < 2.0
    assert algorithm.getProperty("PeakWidth").value == pytest.approx(SIGMA, rel=0.05)


def test_input_is_unchanged(make_detector):
    input_workspace = make_detector(weighted=True)
    n_events = mtd[input_workspace].getNumberEvents()

    run_algorithm(input_workspace)

    workspace = mtd[input_workspace]
    assert isinstance(workspace, IEventWorkspace)
    assert workspace.getNumberEvents() == n_events
    assert workspace.getSpectrum(0).getEventType() == EventType.WEIGHTED
    assert workspace.getAxis(0).getUnit().unitID() == "TOF"


def test_centroid_does_not_depend_on_the_x_unit_or_on_events(make_detector):
    tof = make_detector()
    wavelength = ConvertUnits(InputWorkspace=tof, Target="Wavelength", StoreInADS=False)
    histogram = Rebin(InputWorkspace=tof, Params=[0.0, 100.0, 1000.0], PreserveEvents=False, StoreInADS=False)

    centres = [centre(run_algorithm(workspace)) for workspace in (tof, wavelength, histogram)]

    assert centres == pytest.approx([centres[0]] * 3, abs=1e-6)


def test_profile_counts_every_event(make_detector):
    """Including events outside the X axis and on its edges, which Integration's default range drops."""
    input_workspace = make_detector()
    for tof in (0.0, 1000.0, 1500.0, -5.0):
        mtd[input_workspace].getSpectrum(0).addEventQuickly(tof, DateAndTime(0))

    run_algorithm(input_workspace, output_workspace=PROFILE_WORKSPACE)

    assert mtd[PROFILE_WORKSPACE].readY(0).sum() == mtd[input_workspace].getNumberEvents()


@pytest.mark.parametrize(
    ("peak", "background"), itertools.product(["supergauss", "gauss"], ["linear", "constant", "none"])
)
def test_every_function_choice_fits(make_detector, peak, background):
    algorithm = run_algorithm(make_detector(), peak, background)

    assert abs(centre(algorithm) - CENTRE) < 0.5
    function = algorithm.getProperty("Function").value
    assert function.startswith({"supergauss": "name=SuperGaussian", "gauss": "name=Gaussian"}[peak])
    assert ("Background" in function) == (background != "none")


def test_super_gaussian_sigma_is_root_two_gaussian_sigma_at_exponent_two():
    x = np.linspace(-10.0, 10.0, 41)
    gaussian = FunctionWrapper("Gaussian", Height=5.0, PeakCentre=1.0, Sigma=2.0)
    # Evaluated directly: FunctionWrapper passes parameter values through a six-digit string.
    super_gaussian = SuperGaussian()
    super_gaussian.initialize()
    for name, value in (("Height", 5.0), ("PeakCentre", 1.0), ("Sigma", 2.0 * np.sqrt(2)), ("Exponent", 2.0)):
        super_gaussian.setParameter(name, value)

    np.testing.assert_allclose(super_gaussian.function1D(x), gaussian(x), rtol=1e-12)


def test_exponent_bound_is_enforced(make_detector):
    """A cusp-like peak fits below Exponent 2 unless the bound is enforced."""

    def cusp(y):
        return HEIGHT * np.exp(-(np.abs((y - CENTRE) / SIGMA) ** 1.2)) + 30.0

    algorithm = run_algorithm(make_detector(profile=cusp))
    function = FunctionFactory.createInitialized(algorithm.getProperty("Function").value)

    assert function.getFunction(0).getParameterValue("Exponent") == 2.0


def test_function_string_reproduces_the_fitted_profile(make_detector):
    algorithm = run_algorithm(make_detector(), output_workspace=PROFILE_WORKSPACE)

    profile = mtd[PROFILE_WORKSPACE]
    assert [profile.getAxis(1).label(index) for index in range(3)] == ["Data", "Calc", "Diff"]
    function = FunctionFactory.createInitialized(algorithm.getProperty("Function").value)
    # The function string carries six significant digits.
    np.testing.assert_allclose(FunctionWrapper(function)(profile.readX(1)), profile.readY(1), atol=1e-3 * HEIGHT)


def test_empty_y_pixel_has_unit_error(make_detector):
    run_algorithm(
        make_detector(profile=lambda y: 0.0 if y == 5 else truth(y)),
        output_workspace=PROFILE_WORKSPACE,
    )

    data = mtd[PROFILE_WORKSPACE]
    index = list(data.readX(0)).index(5.0)
    assert data.readY(0)[index] == 0.0
    assert data.readE(0)[index] == 1.0


def test_centroid_does_not_depend_on_the_weight_scale(make_detector):
    """Proton-charge normalization rescales every weight; empty pixels must not change the fit then."""
    # No background: the tails are empty pixels.
    sparse = make_detector(profile=lambda y: HEIGHT / 20 * np.exp(-(np.abs((y - CENTRE) / SIGMA) ** EXPONENT)))
    assert 0.0 in mtd[sparse].extractY()
    centres = [
        centre(run_algorithm(Scale(InputWorkspace=sparse, Factor=factor, StoreInADS=False)))
        for factor in (1.0, 1e-3, 1e3)
    ]

    # Scaled weights are float32, and the minimizer stops within its tolerance.
    assert centres == pytest.approx([centres[0]] * 3, abs=1e-4)


def test_history_records_the_function_choices(make_detector):
    run_algorithm(make_detector(), "gauss", "constant", output_workspace=PROFILE_WORKSPACE)

    history = mtd[PROFILE_WORKSPACE].getHistory().lastAlgorithm()
    assert history.name() == "FitBeamPeak"
    assert history.getPropertyValue("PeakFunction") == "gauss"
    assert history.getPropertyValue("BackgroundFunction") == "constant"


def test_nothing_is_registered_without_an_output_workspace(make_detector):
    before = set(mtd.getObjectNames())

    run_algorithm(make_detector())

    assert set(mtd.getObjectNames()) == before


@pytest.mark.parametrize(
    ("peak", "background", "match"),
    [
        ("supergauss", "linear", "holds 7 y pixels, but fitting supergauss \\+ linear needs at least 8"),
        ("gauss", "none", None),
    ],
)
def test_too_few_y_pixels(make_detector, peak, background, match):
    # Spectra 17-23 hold detector IDs 18-24: y pixels 18 to 24 of the first column, around the peak.
    window = ExtractSpectra(InputWorkspace=make_detector(), WorkspaceIndexList=list(range(17, 24)), StoreInADS=False)

    if match is None:
        run_algorithm(window, peak, background)  # 5 parameters: 7 pixels are enough
    else:
        with pytest.raises(RuntimeError, match=match):
            run_algorithm(window, peak, background)


def test_missing_number_of_y_pixels(make_detector):
    with pytest.raises(RuntimeError, match="number-of-y-pixels"):
        run_algorithm(make_detector(instrument_shape=(N_X, None)))


def test_grouped_spectrum(make_detector):
    input_workspace = make_detector()
    mtd[input_workspace].getSpectrum(0).addDetectorID(2)

    with pytest.raises(RuntimeError, match="exactly one detector pixel"):
        run_algorithm(input_workspace)


def test_no_events(make_detector):
    with pytest.raises(RuntimeError, match="no events in InputWorkspace"):
        run_algorithm(make_detector(profile=lambda _: 0.0))


def test_all_zero_profile(make_detector):
    """Events whose weights cancel: a histogram with no counts."""
    histogram = Rebin(
        InputWorkspace=make_detector(), Params=[0.0, 1000.0, 1000.0], PreserveEvents=False, StoreInADS=False
    )
    for index in range(histogram.getNumberHistograms()):
        histogram.dataY(index)[:] = 0.0

    with pytest.raises(RuntimeError, match="no counts in the 40 y pixels of InputWorkspace"):
        run_algorithm(histogram)


@pytest.mark.parametrize(
    ("status", "converged"),
    [
        ("success", True),
        ("Changes in function value are too small", True),
        ("Changes in parameter value are too small", True),  # what the measured direct beam 198403 ends with
        ("Failed to converge after 500 iterations", False),
        ("Failed to converge, maximum mu reached.", False),
    ],
)
def test_converged_statuses(status, converged):
    assert (status in CONVERGED) == converged


# Template sequence 7 of tests/data/template.xml: (peak y range, background y range, low-resolution x range).
REFERENCE_RUNS = {
    198388: ((134, 148), (131, 151), (70, 182)),  # reflected
    198403: ((138, 147), (135, 150), (107, 133)),  # its direct beam
}


@pytest.mark.datarepo
@pytest.mark.parametrize("run_number", list(REFERENCE_RUNS), ids=["reflected", "direct_beam"])
def test_measured_run(nexus_dir, run_number):
    """Centroid inside the template's peak range and within a pixel of the legacy Gaussian fit."""
    (peak_low, peak_high), (background_low, background_high), low_res = REFERENCE_RUNS[run_number]
    name = f"fit_beam_peak_{run_number}"
    LoadEventNexus(Filename=os.path.join(nexus_dir, f"REF_L_{run_number}.nxs.h5"), OutputWorkspace=name)
    try:
        # The window is the hull of the peak range widened by the default pad of 1 and the background strips.
        window = AlgorithmManager.create("ExtractDetectorROI")
        window.setChild(True)
        window.setProperty("InputWorkspace", name)
        window.setProperty("XRange", list(low_res))
        window.setProperty("YRange", [min(peak_low - 1, background_low), max(peak_high + 1, background_high)])
        window.setProperty("OutputWorkspace", "_")
        window.execute()

        algorithm = run_algorithm(window.getProperty("OutputWorkspace").value, output_workspace=PROFILE_WORKSPACE)
    finally:
        DeleteWorkspace(name)

    assert peak_low <= centre(algorithm) <= peak_high
    data = mtd[PROFILE_WORKSPACE]
    x, y = np.array(data.readX(0)), np.array(data.readY(0))
    legacy_centre, _, _ = fit_signal_flat_bck(x, y, x_min=0, x_max=x.size, center=x[np.argmax(y)])
    assert centre(algorithm) == pytest.approx(legacy_centre.value, abs=1.0)
