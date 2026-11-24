"""Unit tests for the SubtractStripBackground Mantid algorithm and its fit functions."""

import os

import numpy as np
import pytest
from mantid.api import AlgorithmManager
from mantid.simpleapi import (
    ConvertUnits,
    CreateSampleWorkspace,
    DeleteWorkspace,
    GroupDetectors,
    LoadEventNexus,
    Rebin,
    SetInstrumentParameter,
    mtd,
)

from lr_reduction.legacy.background import functional_background
from lr_reduction.operations.algorithms.subtract_strip_background import (
    SubtractStripBackground,
    background_rows,
    evaluate_strip_background,
    fit_strip_background,
    spectrum_rows,
)

# A non-square detector, so that swapping x and y in the row rule cannot go unnoticed.
N_X = 3
N_Y = 30
N_BINS = 50
PEAK = [12, 16]
STRIPS = [3, 9, 19, 25]

OUTPUT_WORKSPACE = "bkg_subtracted"
BACKGROUND_WORKSPACE = "bkg_fitted"


@pytest.fixture(autouse=True)
def clean_output():
    yield
    for name in (OUTPUT_WORKSPACE, BACKGROUND_WORKSPACE):
        if mtd.doesExist(name):
            DeleteWorkspace(name)


@pytest.fixture
def make_rows():
    """Factory for a histogram with one spectrum per detector row, each grouping the row's ``N_X`` pixels.

    ``counts`` and ``errors`` are ``(n_rows, N_BINS)`` arrays for ``rows`` (all rows by default), row ``rows[i]`` taking
    ``counts[i]``. Rows are identified by detector ID, ``y = id mod N_Y``. GroupDetectors orders the spectra by their
    first workspace index, so row 0 comes last and the spectra are not in ascending ``y``. The instrument's
    ``number-of-y-pixels`` is ``y_pixels`` (default ``N_Y``; None leaves it undefined). Deletes each workspace
    on teardown.
    """
    names = []

    def make(counts, errors, rows=None, unit="TOF", y_pixels=N_Y, n_y=N_Y):
        rows = list(range(n_y)) if rows is None else list(rows)
        name = mtd.unique_hidden_name()
        CreateSampleWorkspace(
            Function="Flat background",
            NumBanks=N_X * n_y,
            BankPixelWidth=1,
            XMin=0.0,
            XMax=float(N_BINS),
            BinWidth=1.0,
            OutputWorkspace=name,
        )
        # CreateSampleWorkspace numbers detector IDs from 1, at workspace index id - 1
        ids = range(1, N_X * n_y + 1)
        pattern = ",".join("+".join(str(d - 1) for d in ids if d % n_y == row) for row in rows)
        GroupDetectors(InputWorkspace=name, OutputWorkspace=name, GroupingPattern=pattern)
        if y_pixels is not None:
            SetInstrumentParameter(
                Workspace=name, ParameterName="number-of-y-pixels", ParameterType="Number", Value=f"{float(y_pixels)}"
            )
        workspace = mtd[name]
        workspace.getAxis(0).setUnit(unit)
        position = {row: i for i, row in enumerate(rows)}
        for index, row in enumerate(spectrum_rows(workspace, n_y)):
            workspace.setSharedY(index, np.asarray(counts[position[row]], dtype=float))
            workspace.setSharedE(index, np.asarray(errors[position[row]], dtype=float))
        names.append(name)
        return name

    yield make
    for name in names:
        if mtd.doesExist(name):
            DeleteWorkspace(name)


def run_algorithm(input_workspace, peak=PEAK, strips=STRIPS, model="Linear", output=OUTPUT_WORKSPACE, background=None):
    algorithm = AlgorithmManager.create("SubtractStripBackground")
    algorithm.setProperty("InputWorkspace", input_workspace)
    algorithm.setProperty("PeakRange", peak)
    algorithm.setProperty("BackgroundRanges", strips)
    algorithm.setProperty("BackgroundModel", model)
    algorithm.setProperty("OutputWorkspace", output)
    if background:
        algorithm.setProperty("BackgroundWorkspace", background)
    algorithm.execute()
    return mtd[output], algorithm


def validation_issues(input_workspace, peak=PEAK, strips=STRIPS, model="Linear"):
    algorithm = SubtractStripBackground()
    algorithm.initialize()
    algorithm.setProperty("InputWorkspace", input_workspace)
    algorithm.setProperty("PeakRange", peak)
    algorithm.setProperty("BackgroundRanges", strips)
    algorithm.setProperty("BackgroundModel", model)
    algorithm.setProperty("OutputWorkspace", OUTPUT_WORKSPACE)
    return algorithm.validateInputs()


def linear_background(rows, slope=0.5, intercept=20.0):
    """A background that varies across y and λ: ``(slope·y + intercept)·(1 + j/N_BINS)``."""
    rows = np.asarray(rows, dtype=float)[:, None]
    return (slope * rows + intercept) * (1.0 + np.arange(N_BINS) / N_BINS)


def peak_signal(rows):
    rows = np.asarray(rows, dtype=float)[:, None]
    profile = np.where((rows >= PEAK[0]) & (rows <= PEAK[1]), 1000.0 * np.exp(-0.5 * (rows - 14.0) ** 2 / 2.0), 0.0)
    return profile * np.ones(N_BINS)


class TestFitFunctions:
    @pytest.mark.parametrize("model", ["Linear", "Constant"])
    def test_closed_form_matches_polyfit_and_matrix_covariance(self, model):
        rng = np.random.default_rng(1)
        y = np.array([2, 3, 4, 7, 20, 21, 22, 25])
        r = rng.normal(50.0, 5.0, (y.size, 7))
        sigma = rng.uniform(1.0, 6.0, r.shape)
        a, c, cov = fit_strip_background(y, r, sigma, model)

        degree = 1 if model == "Linear" else 0
        design = np.vander(y.astype(float), degree + 1)
        xtx_inv = np.linalg.inv(design.T @ design)
        for j in range(r.shape[1]):
            coefficients = np.polyfit(y, r[:, j], degree)
            expected_cov = xtx_inv @ design.T @ np.diag(sigma[:, j] ** 2) @ design @ xtx_inv
            if model == "Linear":
                np.testing.assert_allclose([a[j], c[j]], coefficients)
                np.testing.assert_allclose(cov[j], expected_cov)
            else:
                assert a[j] == 0.0
                np.testing.assert_allclose(c[j], coefficients[0])
                np.testing.assert_allclose(cov[j, 1, 1], expected_cov[0, 0])
                np.testing.assert_allclose([cov[j, 0, 0], cov[j, 0, 1]], 0.0)

    def test_unknown_model_rejected(self):
        with pytest.raises(ValueError, match="'linear'"):
            fit_strip_background([1, 2], np.ones((2, 3)), np.ones((2, 3)), "linear")

    def test_evaluate_includes_covariance(self):
        cov = np.array([[[0.04, -0.3], [-0.3, 4.0]]])
        bkg, var = evaluate_strip_background((np.array([2.0]), np.array([1.0]), cov), [0.0, 10.0])
        np.testing.assert_allclose(bkg[:, 0], [1.0, 21.0])
        np.testing.assert_allclose(var[:, 0], [4.0, 100 * 0.04 - 2 * 10 * 0.3 + 4.0])

    @pytest.mark.parametrize("level", [0.3, 1.0])
    def test_low_counts_unbiased_with_correct_variance(self, level):
        rng = np.random.default_rng(2)
        y = np.array(background_rows(PEAK, [STRIPS[:2], STRIPS[2:]]))
        trials = 4000
        r = rng.poisson(level, (y.size, trials)).astype(float)
        fit = fit_strip_background(y, r, np.sqrt(r), "Linear")
        bkg, var = evaluate_strip_background(fit, [14])
        mean, spread = bkg.mean(), bkg.std()
        assert abs(mean - level) < 4 * spread / np.sqrt(trials)
        assert var.mean() == pytest.approx(spread**2, rel=0.1)


class TestBackgroundRows:
    @pytest.mark.parametrize(
        ("strips", "expected"),
        [
            ([[3, 9], [19, 25]], [*range(3, 10), *range(19, 26)]),  # one strip on each side
            ([[3, 9]], list(range(3, 10))),  # one side only
            ([[5, 22]], [*range(5, 12), *range(17, 23)]),  # one range spanning the peak
            ([[8, 14]], list(range(8, 12))),  # overlapping the peak on one side
        ],
    )
    def test_pixel_set(self, make_rows, strips, expected):
        counts = linear_background(range(N_Y))
        _, algorithm = run_algorithm(make_rows(counts, np.sqrt(counts)), strips=[v for s in strips for v in s])
        pixels = list(algorithm.getProperty("BackgroundPixels").value)
        assert pixels == expected
        assert PEAK[0] not in pixels and PEAK[1] not in pixels


class TestSubtractStripBackground:
    def test_noiseless_linear_background_removed_exactly(self, make_rows):
        rows = np.arange(N_Y)
        counts = linear_background(rows) + peak_signal(rows)
        output, _ = run_algorithm(make_rows(counts, np.sqrt(counts)))
        np.testing.assert_allclose(output.extractY(), peak_signal(rows)[PEAK[0] : PEAK[1] + 1], atol=1e-9)

    def test_constant_model_is_mean_of_background_rows(self, make_rows):
        rng = np.random.default_rng(3)
        counts = rng.uniform(5.0, 15.0, (N_Y, N_BINS))
        name = make_rows(counts, np.sqrt(counts))
        run_algorithm(name, model="Constant", background=BACKGROUND_WORKSPACE)
        bkg_rows = background_rows(PEAK, [STRIPS[:2], STRIPS[2:]])
        np.testing.assert_allclose(mtd[BACKGROUND_WORKSPACE].extractY(), np.tile(counts[bkg_rows].mean(axis=0), (5, 1)))

    def test_recovery_pulls_have_unit_spread(self, make_rows):
        rng = np.random.default_rng(4)
        rows = np.arange(N_Y)
        truth = linear_background(rows, slope=20.0, intercept=500.0) + peak_signal(rows)
        counts = rng.poisson(truth).astype(float)
        output, _ = run_algorithm(make_rows(counts, np.sqrt(counts)))
        expected = peak_signal(rows)[PEAK[0] : PEAK[1] + 1]
        pulls = (output.extractY() - expected) / output.extractE()
        assert abs(pulls.mean()) < 0.15
        assert pulls.std() == pytest.approx(1.0, abs=0.1)

    def test_unit_agnostic(self, make_rows):
        rng = np.random.default_rng(5)
        counts = rng.poisson(linear_background(range(N_Y))).astype(float)
        tof, _ = run_algorithm(make_rows(counts, np.sqrt(counts), unit="TOF"), output=OUTPUT_WORKSPACE)
        lam, _ = run_algorithm(make_rows(counts, np.sqrt(counts), unit="Wavelength"), output=BACKGROUND_WORKSPACE)
        np.testing.assert_array_equal(tof.extractY(), lam.extractY())
        np.testing.assert_array_equal(tof.extractE(), lam.extractE())
        assert lam.getAxis(0).getUnit().unitID() == "Wavelength"

    def test_output_shape_identity_and_background_workspace(self, make_rows):
        rng = np.random.default_rng(6)
        rows = np.arange(N_Y)
        counts = rng.poisson(linear_background(rows) + peak_signal(rows)).astype(float)
        peak = [0, 2]  # row 0 is the input's last spectrum, so the output order is [1, 2, 0]
        name = make_rows(counts, np.sqrt(counts), unit="Wavelength")
        output, _ = run_algorithm(name, peak=peak, strips=[3, 9], background=BACKGROUND_WORKSPACE)
        input_workspace = mtd[name]
        input_rows = spectrum_rows(input_workspace, N_Y)
        peak_indices = [i for i, row in enumerate(input_rows) if peak[0] <= row <= peak[1]]
        assert spectrum_rows(output, N_Y) == [1, 2, 0]

        assert output.getNumberHistograms() == len(peak_indices)
        for k, index in enumerate(peak_indices):
            assert list(output.getSpectrum(k).getDetectorIDs()) == list(
                input_workspace.getSpectrum(index).getDetectorIDs()
            )
            np.testing.assert_array_equal(output.x(k), input_workspace.x(index))
        assert output.getAxis(0).getUnit().unitID() == "Wavelength"
        assert list(output.getInstrument().getNumberParameter("number-of-y-pixels")) == [float(N_Y)]

        background = mtd[BACKGROUND_WORKSPACE]
        counts_in = input_workspace.extractY()[peak_indices]
        np.testing.assert_allclose(counts_in - background.extractY(), output.extractY())
        np.testing.assert_allclose(np.sqrt(counts_in + background.extractE() ** 2), output.extractE())

    def test_input_untouched_and_in_place(self, make_rows):
        counts = linear_background(range(N_Y)) + peak_signal(range(N_Y))
        name = make_rows(counts, np.sqrt(counts))
        run_algorithm(name)
        rows = spectrum_rows(mtd[name], N_Y)
        np.testing.assert_array_equal(mtd[name].extractY(), counts[rows])
        np.testing.assert_array_equal(mtd[name].extractE(), np.sqrt(counts)[rows])

        output, _ = run_algorithm(name, output=name)
        assert output.getNumberHistograms() == PEAK[1] - PEAK[0] + 1


class TestValidation:
    def test_non_finite_values(self, make_rows):
        counts = linear_background(range(N_Y))
        counts[3, 4] = np.nan
        assert "finite" in validation_issues(make_rows(counts, np.ones_like(counts)))["InputWorkspace"]

    def test_negative_errors(self, make_rows):
        counts = linear_background(range(N_Y))
        assert "non-negative" in validation_issues(make_rows(counts, -np.ones_like(counts)))["InputWorkspace"]

    def test_missing_y_pixel_parameter(self, make_rows):
        counts = linear_background(range(N_Y))
        name = make_rows(counts, np.sqrt(counts), y_pixels=None)
        assert "number-of-y-pixels" in validation_issues(name)["InputWorkspace"]

    def test_spectrum_spanning_two_rows(self, make_rows):
        counts = linear_background(range(N_Y))
        # a 25-row detector setting makes each 3-pixel group of the 30-row detector straddle rows
        name = make_rows(counts, np.sqrt(counts), y_pixels=25)
        assert "single y row" in validation_issues(name, peak=[1, 2], strips=[3, 5])["InputWorkspace"]

    def test_two_spectra_on_one_row(self, make_rows):
        counts = linear_background(range(N_Y))
        # a 15-row detector setting maps rows r and r + 15 of the 30-row detector to the same y
        name = make_rows(counts, np.sqrt(counts), y_pixels=15)
        assert "y rows [0, 1, 2" in validation_issues(name, peak=[1, 2], strips=[3, 5])["InputWorkspace"]

    def test_missing_rows(self, make_rows):
        rows = range(10, N_Y)
        counts = linear_background(rows)
        issues = validation_issues(make_rows(counts, np.sqrt(counts), rows=rows))
        assert "absent" in issues["BackgroundRanges"] and "[3, 9]" in issues["BackgroundRanges"]

    def test_bad_strip_count(self, make_rows):
        counts = linear_background(range(N_Y))
        assert (
            "one or two"
            in validation_issues(make_rows(counts, np.sqrt(counts)), strips=[1, 2, 3, 4, 5, 6])["BackgroundRanges"]
        )

    @pytest.mark.parametrize(("model", "strips"), [("Linear", [11, 13]), ("Constant", [12, 16])])
    def test_too_few_background_rows(self, make_rows, model, strips):
        counts = linear_background(range(N_Y))
        issues = validation_issues(make_rows(counts, np.sqrt(counts)), strips=strips, model=model)
        assert f"the {model} model needs" in issues["BackgroundRanges"]

    def test_event_workspace_rejected(self):
        name = mtd.unique_hidden_name()
        CreateSampleWorkspace(WorkspaceType="Event", NumBanks=1, BankPixelWidth=1, OutputWorkspace=name)
        try:
            algorithm = SubtractStripBackground()
            algorithm.initialize()
            algorithm.setProperty("InputWorkspace", name)
            algorithm.setProperty("PeakRange", PEAK)
            algorithm.setProperty("BackgroundRanges", STRIPS)
            algorithm.setProperty("OutputWorkspace", OUTPUT_WORKSPACE)
            assert "event workspace" in algorithm.validateInputs()["InputWorkspace"]
        finally:
            DeleteWorkspace(name)


class LegacyReflectivityStub:
    """Stands in for the legacy ``EventReflectivity``: returns the per-row counts of a row range, per unit charge."""

    theta = 0.0

    def __init__(self, workspace, first_row):
        self.charge = workspace.getRun().getProtonCharge()
        self.counts, self.errors, self.first_row = workspace.extractY(), workspace.extractE(), first_row

    def _reflectivity(self, _workspace, peak, **_kwargs):
        rows = slice(peak[0] - self.first_row, peak[1] - self.first_row + 1)
        return self.counts[rows] / self.charge, self.errors[rows] / self.charge


@pytest.mark.datarepo
def test_matches_legacy_functional_background(nexus_dir):
    """On a measured run, the peak-summed background matches the legacy ``functional_background``
    on the same λ bins, within the legacy fit's error, where the strips hold many counts per cell."""
    peak, strip, x_range = (
        [134, 148],
        [131, 151],
        (70, 182),
    )  # REF_L_198388_auto_template.xml, one range spanning the peak
    name = mtd.unique_hidden_name()
    LoadEventNexus(Filename=os.path.join(nexus_dir, "REF_L_198388.nxs.h5"), OutputWorkspace=name)
    try:
        rows = range(strip[0], strip[1] + 1)
        pattern = ",".join("+".join(str(x * 304 + y) for x in range(x_range[0], x_range[1] + 1)) for y in rows)
        GroupDetectors(InputWorkspace=name, OutputWorkspace=name, GroupingPattern=pattern)
        ConvertUnits(InputWorkspace=name, OutputWorkspace=name, Target="Wavelength")
        Rebin(InputWorkspace=name, OutputWorkspace=name, Params="2.5,0.1,6.5", PreserveEvents=False)
        run_algorithm(name, peak=peak, strips=strip, background=BACKGROUND_WORKSPACE)

        stub = LegacyReflectivityStub(mtd[name], strip[0])
        legacy, legacy_error = functional_background(mtd[name], stub, peak, [*strip, 0, 0], x_range)
        legacy, legacy_error = legacy * stub.charge, legacy_error * stub.charge
        summed = mtd[BACKGROUND_WORKSPACE].extractY().sum(axis=0)

        strip_rows = [y - strip[0] for y in background_rows(peak, [strip])]
        high_counts = stub.counts[strip_rows].min(axis=0) >= 20
        assert high_counts.sum() >= 10
        np.testing.assert_array_less(np.abs(summed - legacy)[high_counts], legacy_error[high_counts])
    finally:
        DeleteWorkspace(name)
