"""Unit tests for the MergeDirectBeamSegments Mantid algorithm and its kernels."""

from collections.abc import Callable, Iterator

import numpy as np
import pytest
from mantid.api import EventType
from mantid.kernel import DateAndTime
from mantid.simpleapi import (
    CompressEvents,
    ConvertUnits,
    CreateSampleWorkspace,
    DeleteWorkspace,
    Rebin,
    SumSpectra,
    mtd,
)

from lr_reduction.operations.algorithms.merge_direct_beam_segments import (
    MergeDirectBeamSegments,
    coverage_gaps,
    event_weight,
    in_band_mask,
    overlap_weights,
    summed_weights_and_counts,
    weighting_grid,
)

OUTPUT_WORKSPACE = "merged_direct_beam"
WEIGHTS_WORKSPACE = "merged_direct_beam_weights"


# --------------------------------------------------------------------------------------------
# Kernels


def test_coverage_gaps_finds_gap_between_bands():
    assert coverage_gaps([(4.0, 6.0), (1.0, 3.0)]) == [(3.0, 4.0, 1, 0)]


@pytest.mark.parametrize("bands", [[(1.0, 3.0), (2.0, 5.0)], [(1.0, 5.0), (2.0, 3.0), (3.0, 4.0)], [(1.0, 2.0)]])
def test_coverage_gaps_none_for_overlapping_or_nested_bands(bands):
    assert coverage_gaps(bands) == []


def test_weighting_grid_contains_every_band_edge():
    bands = [(1.0, 3.013), (2.0071, 5.0)]
    edges = weighting_grid(bands, 0.1)
    assert np.all(np.diff(edges) > 0)
    for low, high in bands:
        assert low in edges
        assert np.nextafter(high, np.inf) in edges
    assert edges[0] == 1.0
    assert edges[-1] == np.nextafter(5.0, np.inf)


def test_in_band_mask_includes_upper_edge_bin():
    bands = [(1.0, 2.0), (1.5, 3.0)]
    edges = weighting_grid(bands, 0.25)
    mask = in_band_mask(edges, bands[0])
    centers = 0.5 * (edges[:-1] + edges[1:])
    assert np.all(centers[mask] <= 2.0)
    assert np.all(centers[~mask] > 2.0)
    assert edges[1:][mask][-1] == np.nextafter(2.0, np.inf)


def test_event_weight_is_exact_regardless_of_counts():
    c = 0.25
    counts = np.array([3, 1, 7, 2])
    in_band = np.ones(4, dtype=bool)
    result = event_weight(counts * c, counts, in_band, np.arange(4.0))
    np.testing.assert_allclose(result, c)


def test_event_weight_is_mean_weight_of_the_bin_events():
    # Events of weight 1 and 3 in one bin: the event weight is their mean, 2.
    result = event_weight(np.array([4.0]), np.array([2]), np.array([True]), np.array([0.0]))
    np.testing.assert_allclose(result, [2.0])


def test_event_weight_interpolates_inverse_across_empty_bins_and_is_infinite_outside_band():
    # Event weights 1 and 1/3 on either side of an empty bin: 1/c interpolates to 2, so c = 1/2.
    y = np.array([2.0, 0.0, 1.0, 0.0])
    n = np.array([2, 0, 3, 0])
    in_band = np.array([True, True, True, False])
    result = event_weight(y, n, in_band, np.arange(4.0))
    np.testing.assert_allclose(result[:3], [1.0, 0.5, 1 / 3])
    assert result[3] == np.inf


def test_event_weight_is_infinite_on_zero_weight_bins_and_empty_bins_next_to_them():
    # Bin 1 is empty between valid bins; bin 3 holds only zero-weight events; bin 4 is empty beside it.
    y = np.array([2.0, 0.0, 3.0, 0.0, 0.0, 1.0])
    n = np.array([2, 0, 3, 4, 0, 1])
    in_band = np.ones(6, dtype=bool)
    result = event_weight(y, n, in_band, np.arange(6.0))
    np.testing.assert_allclose(result[[0, 1, 2, 5]], 1.0)
    assert np.all(result[[3, 4]] == np.inf)


def test_overlap_weights_finite_where_one_segment_has_zero_weight_events():
    in_band = np.ones(3, dtype=bool)
    centers = np.arange(3.0)
    a = event_weight(np.array([1.0, 0.0, 1.0]), np.array([1, 5, 1]), in_band, centers)
    b = event_weight(np.array([1.0, 1.0, 1.0]), np.array([1, 1, 1]), in_band, centers)
    weights = overlap_weights([a, b])
    np.testing.assert_allclose(weights, [[0.5, 0.0, 0.5], [0.5, 1.0, 0.5]])


def test_overlap_weights_sum_to_one_where_covered():
    event_weights = np.array([[0.5, 0.5, np.inf, np.inf], [np.inf, 1 / 6, 1.0, np.inf]])
    weights = overlap_weights(event_weights)
    np.testing.assert_allclose(weights, [[1.0, 0.25, 0.0, 0.0], [0.0, 0.75, 1.0, 0.0]])


def test_summed_weights_and_counts_over_all_spectra(make_segment):
    name = make_segment(2.0, 4.0, c=0.5, density=1.0e3, seed=14)
    workspace = mtd[name]
    edges = np.array([2.0, 3.0, np.nextafter(workspace.getTofMax(), np.inf)])
    summed, counts = summed_weights_and_counts(workspace, edges)
    tofs = np.concatenate([workspace.getSpectrum(i).getTofs() for i in range(workspace.getNumberHistograms())])
    assert counts.sum() == workspace.getNumberEvents() == len(tofs)
    np.testing.assert_array_equal(counts, [np.sum(tofs < 3.0), np.sum(tofs >= 3.0)])
    np.testing.assert_allclose(summed, 0.5 * counts)


# --------------------------------------------------------------------------------------------
# Algorithm


@pytest.fixture
def workspaces() -> Iterator[list[str]]:
    """Collects workspace names to delete on teardown, starting with the outputs."""
    names = [OUTPUT_WORKSPACE, WEIGHTS_WORKSPACE]
    yield names
    for name in names:
        if mtd.doesExist(name):
            DeleteWorkspace(name)


@pytest.fixture
def make_segment(workspaces: list[str]) -> Callable[..., str]:
    """Factory for a four-pixel event workspace in wavelength with events of weight `c`.

    Events are drawn from a flat flux: the expected raw count density is `density / c` per Å over
    `[low, high]`, so each segment is an unbiased estimate of the same flux `density`.
    """

    def make(
        low: float,
        high: float,
        c: float,
        density: float = 2.0e4,
        seed: int = 0,
        run_number: int = 1,
        zero: tuple[float, float] | None = None,
        name: str | None = None,
    ) -> str:
        """`zero`, a ``(min, max)`` range, gives its events weight 0, as the dead time threshold does."""
        rng = np.random.default_rng(seed)
        count = rng.poisson(density / c * (high - low))
        wavelengths = np.sort(rng.uniform(low, high, count))
        name = name or mtd.unique_hidden_name()
        CreateSampleWorkspace(
            WorkspaceType="Event",
            NumBanks=1,
            BankPixelWidth=2,
            NumEvents=0,
            XMin=0.0,
            XMax=20000.0,
            BinWidth=20000.0,
            OutputWorkspace=name,
        )
        ConvertUnits(InputWorkspace=name, Target="Wavelength", OutputWorkspace=name)
        workspace = mtd[name]
        workspace.getRun().addProperty("run_number", str(run_number), True)
        spectra = [workspace.getSpectrum(index) for index in range(workspace.getNumberHistograms())]
        for spectrum in spectra:
            spectrum.switchTo(EventType.WEIGHTED)
        for index, wavelength in enumerate(wavelengths):
            weight = 0.0 if zero is not None and zero[0] <= wavelength <= zero[1] else c
            spectra[index % len(spectra)].addWeightedEventQuickly(wavelength, weight, weight**2, DateAndTime(0))
        workspaces.append(name)
        return name

    return make


def run_merge(names: list[str], **properties: object) -> MergeDirectBeamSegments:
    algorithm = MergeDirectBeamSegments()
    algorithm.initialize()
    algorithm.setRethrows(True)
    algorithm.setProperty("InputWorkspaces", names)
    algorithm.setProperty("OutputWorkspace", OUTPUT_WORKSPACE)
    algorithm.setProperty("OutputWeightsWorkspace", WEIGHTS_WORKSPACE)
    for name, value in properties.items():
        algorithm.setProperty(name, value)
    algorithm.execute()
    return algorithm


def histogram(name: str, params: list[float]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    aligned = Rebin(InputWorkspace=name, Params=params, PreserveEvents=True, StoreInADS=False)
    summed = SumSpectra(InputWorkspace=aligned, StoreInADS=False)
    rebinned = Rebin(InputWorkspace=summed, Params=params, PreserveEvents=False, StoreInADS=False)
    return rebinned.readX(0), rebinned.readY(0), rebinned.readE(0)


def test_overlap_is_averaged_not_summed(make_segment):
    density = 2.0e4
    a = make_segment(2.0, 6.0, c=0.5, density=density, seed=1)
    b = make_segment(4.0, 8.0, c=2.0, density=density, seed=2)
    run_merge([a, b])

    edges, y, e = histogram(OUTPUT_WORKSPACE, [2.0, 0.5, 8.0])
    flux = y / np.diff(edges)
    sigma = e / np.diff(edges)
    # Flat at the common flux across A-only, overlap, and B-only bins; summing would give 2× in the overlap.
    assert np.all(np.abs(flux - density) < 5 * sigma), flux
    assert mtd[OUTPUT_WORKSPACE].getNumberEvents() == mtd[a].getNumberEvents() + mtd[b].getNumberEvents()


def test_overlap_bin_is_pooled_counts_with_combined_event_weight(make_segment):
    c_a, c_b = 0.5, 2.0
    a = make_segment(2.0, 6.0, c=c_a, seed=3)
    b = make_segment(4.0, 8.0, c=c_b, seed=4)
    run_merge([a, b])

    params = [4.5, 0.5, 5.5]
    _, y_a, _ = histogram(a, params)
    _, y_b, _ = histogram(b, params)
    _, y, e = histogram(OUTPUT_WORKSPACE, params)
    pooled_counts = y_a / c_a + y_b / c_b
    combined_event_weight = 1 / (1 / c_a + 1 / c_b)
    np.testing.assert_allclose(y, combined_event_weight * pooled_counts, rtol=1e-5)
    np.testing.assert_allclose(e**2, combined_event_weight**2 * pooled_counts, rtol=1e-5)


def test_gap_is_reported_and_holds_no_flux(make_segment):
    a = make_segment(2.0, 4.0, c=1.0, seed=5)
    b = make_segment(5.0, 7.0, c=1.0, seed=6)
    algorithm = run_merge([a, b])

    gap_min, gap_max = algorithm.getProperty("GapMin").value, algorithm.getProperty("GapMax").value
    assert len(gap_min) == 1
    assert gap_min[0] == pytest.approx(4.0, abs=1e-3)
    assert gap_max[0] == pytest.approx(5.0, abs=1e-3)
    _, y, _ = histogram(OUTPUT_WORKSPACE, [4.1, 0.8, 4.9])
    assert y[0] == 0.0
    # Alone everywhere, so neither segment is scaled.
    np.testing.assert_array_equal(np.unique(mtd[WEIGHTS_WORKSPACE].extractY()), [0.0, 1.0])


def test_band_edge_events_are_weighted(make_segment):
    a = make_segment(2.0, 6.0, c=1.0, seed=7)
    b = make_segment(4.0, 8.0, c=1.0, seed=8)
    high_a = mtd[a].getTofMax()
    run_merge([a, b])

    # The last event of A sits in the overlap and must carry weight 1/2, not 1.
    weights = [event_weight for event_weight in _weights_at(OUTPUT_WORKSPACE, high_a)]
    assert weights
    np.testing.assert_allclose(weights, 0.5, rtol=0.05)


def _weights_at(name: str, wavelength: float) -> Iterator[float]:
    workspace = mtd[name]
    for index in range(workspace.getNumberHistograms()):
        spectrum = workspace.getSpectrum(index)
        tofs, weights = spectrum.getTofs(), spectrum.getWeights()
        yield from weights[tofs == wavelength]


def _all_weights(name: str) -> np.ndarray:
    workspace = mtd[name]
    return np.concatenate([workspace.getSpectrum(i).getWeights() for i in range(workspace.getNumberHistograms())])


def test_zero_weight_events_in_overlap_leave_other_segment_at_full_weight(make_segment):
    density = 2.0e4
    a = make_segment(2.0, 6.0, c=1.0, density=density, seed=16, zero=(4.5, 5.0))
    b = make_segment(4.0, 8.0, c=1.0, density=density, seed=17)
    run_merge([a, b])

    assert np.all(np.isfinite(_all_weights(OUTPUT_WORKSPACE)))
    assert np.all(np.isfinite(mtd[WEIGHTS_WORKSPACE].extractY()))
    edges, y, e = histogram(OUTPUT_WORKSPACE, [4.6, 0.3, 4.9])
    flux, sigma = y / np.diff(edges), e / np.diff(edges)
    assert np.all(np.abs(flux - density) < 5 * sigma), flux


def test_zero_weight_events_where_one_segment_is_alone_give_zero_flux(make_segment):
    a = make_segment(2.0, 6.0, c=1.0, seed=18, zero=(2.5, 3.0))
    b = make_segment(4.0, 8.0, c=1.0, seed=19)
    run_merge([a, b])

    assert np.all(np.isfinite(_all_weights(OUTPUT_WORKSPACE)))
    assert np.all(np.isfinite(mtd[WEIGHTS_WORKSPACE].extractY()))
    _, y, _ = histogram(OUTPUT_WORKSPACE, [2.6, 0.3, 2.9])
    assert y[0] == 0.0


def test_input_named_like_a_former_temporary_is_kept(make_segment):
    a = make_segment(2.0, 6.0, c=0.5, seed=20)
    b = make_segment(4.0, 8.0, c=2.0, seed=21, name=f"__{OUTPUT_WORKSPACE}_weighted_0")
    events = mtd[b].getNumberEvents()
    run_merge([a, b])
    assert mtd.doesExist(b)
    assert mtd[b].getNumberEvents() == events
    assert mtd[OUTPUT_WORKSPACE].getNumberEvents() == mtd[a].getNumberEvents() + events


def test_one_segment_is_cloned(make_segment):
    a = make_segment(2.0, 6.0, c=0.5, seed=9)
    run_merge([a])
    np.testing.assert_array_equal(histogram(OUTPUT_WORKSPACE, [2.0, 0.1, 6.0])[1], histogram(a, [2.0, 0.1, 6.0])[1])


def test_inputs_unchanged_and_no_intermediates_left(make_segment):
    a = make_segment(2.0, 6.0, c=0.5, seed=10)
    b = make_segment(4.0, 8.0, c=2.0, seed=11)
    before = histogram(a, [2.0, 0.1, 8.0])[1].copy()
    names_before = set(mtd.getObjectNames())
    run_merge([a, b])
    np.testing.assert_array_equal(histogram(a, [2.0, 0.1, 8.0])[1], before)
    assert set(mtd.getObjectNames()) - names_before == {OUTPUT_WORKSPACE, WEIGHTS_WORKSPACE}


def test_output_axis_is_single_bin(make_segment):
    a = make_segment(2.0, 6.0, c=0.5, seed=12)
    b = make_segment(4.0, 8.0, c=2.0, seed=13)
    run_merge([a, b])
    x = mtd[OUTPUT_WORKSPACE].readX(0)
    assert len(x) == 2
    assert x[0] == mtd[a].getTofMin()
    assert x[1] == np.nextafter(mtd[b].getTofMax(), np.inf)


def test_rejects_compressed_events(make_segment):
    name = make_segment(2.0, 6.0, c=0.5, seed=15)
    CompressEvents(InputWorkspace=name, Tolerance=0.01, OutputWorkspace=name)
    with pytest.raises(RuntimeError, match="compressed"):
        run_merge([name])


def test_rejects_tof_input(workspaces):
    name = mtd.unique_hidden_name()
    CreateSampleWorkspace(WorkspaceType="Event", NumBanks=1, BankPixelWidth=2, OutputWorkspace=name)
    workspaces.append(name)
    with pytest.raises(RuntimeError, match="Wavelength"):
        run_merge([name])
