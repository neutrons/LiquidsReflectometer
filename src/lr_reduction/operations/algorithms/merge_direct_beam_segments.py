"""
Direct-beam segment merging: combine the corrected direct-beam segments of one composite into a
single event workspace, averaging their overlapping wavelength bands.

The event weight ``c_i(lambda)`` is the factor the per-run corrections put on each raw event of
segment ``i``. Every event of segment ``i`` is scaled by an overlap weight inversely proportional to
its event weight, ``alpha_i(lambda) = (1 / c_i(lambda)) / sum_k (1 / c_k(lambda))``, and the
weighted segments are appended with `MergeRuns`. The event weight is measured from the data as the
mean weight of a bin's events, ``c_i = Y / N`` on the segment's ROI-summed spectrum, which is exact
where ``c_i`` is constant across a bin.
"""

import numpy as np
from mantid.api import (
    AlgorithmFactory,
    AnalysisDataService,
    EventType,
    IEventWorkspace,
    IEventWorkspaceProperty,
    MatrixWorkspaceProperty,
    PropertyMode,
    PythonAlgorithm,
)
from mantid.kernel import Direction, FloatArrayProperty, FloatBoundedValidator, StringArrayProperty
from mantid.simpleapi import (
    CloneWorkspace,
    CreateWorkspace,
    DeleteWorkspace,
    MergeRuns,
    Multiply,
    Rebin,
    RebinToWorkspace,
    SumSpectra,
    logger,
)


def coverage_gaps(bands):
    """Ranges between the union's minimum and maximum that no band contains.

    Parameters
    ----------
    bands
        ``(lambda_min, lambda_max)`` per segment, inclusive at both ends.

    Returns
    -------
        ``(gap_min, gap_max, left, right)`` per gap, in increasing wavelength: the gap is the open
        range ``(gap_min, gap_max)``, and ``left``/``right`` are the indices of the bands that end
        and begin at its edges.
    """
    order = sorted(range(len(bands)), key=lambda index: bands[index][0])
    gaps = []
    reach, reach_index = bands[order[0]][1], order[0]
    for index in order[1:]:
        low, high = bands[index]
        if low > reach:
            gaps.append((reach, low, reach_index, index))
        if high > reach:
            reach, reach_index = high, index
    return gaps


def weighting_grid(bands, step):
    """Bin edges of step `step` that also contain every band edge.

    Each upper band edge is moved to the next float up, so the band's last event falls inside a
    half-open bin of its own band. Every bin then lies wholly inside or wholly outside each band.
    """
    lows = np.array([low for low, _ in bands], dtype=float)
    highs = np.nextafter(np.array([high for _, high in bands], dtype=float), np.inf)
    start, stop = lows.min(), highs.max()
    regular = start + step * np.arange(int(np.ceil((stop - start) / step)))
    return np.unique(np.concatenate([regular, lows, highs]))


def in_band_mask(edges, band):
    """Bins of `edges` that lie inside `band`, whose upper edge is inclusive (see `weighting_grid`)."""
    low, high = band[0], np.nextafter(band[1], np.inf)
    return (edges[:-1] >= low) & (edges[1:] <= high)


def summed_weights_and_counts(ws, edges):
    """Summed event weight ``Y`` and number of events ``N`` per bin of `edges`, over all spectra of `ws`.

    The spectra are summed on events into one spectrum, whose events are then histogrammed onto
    `edges`. `SumSpectra` needs a common axis, so the events are first given the single bin
    `[edges[0], edges[-1]]` with `Rebin(PreserveEvents=True)`, which keeps every event. The last
    bin is closed on the right, as in `numpy.histogram`.
    """
    low, high = edges[0], edges[-1]
    common = Rebin(InputWorkspace=ws, Params=[low, high - low, high], PreserveEvents=True, StoreInADS=False)
    spectrum = SumSpectra(InputWorkspace=common, StoreInADS=False).getSpectrum(0)
    wavelengths = spectrum.getTofs()
    summed = np.histogram(wavelengths, edges, weights=spectrum.getWeights())[0].astype(float)
    counts = np.histogram(wavelengths, edges)[0].astype(float)
    return summed, counts


def event_weight(y, n, in_band, centers):
    """Event weight ``Y / N``, the mean weight of a bin's events, per bin of a segment's band.

    In-band bins without events take a ``1 / c`` interpolated linearly in wavelength from the
    neighboring bins. Outside the band the event weight is infinite: the segment contributes
    ``1 / c = 0`` there. Bins whose events all carry weight 0, as the dead time correction's
    threshold leaves them, hold no valid data and are treated like bins outside the band, as are
    empty bins next to them.

    Parameters
    ----------
    y, n
        Summed event weight and number of events, per bin, of the segment's ROI-summed spectrum.
    in_band
        Boolean mask of the bins inside the segment's band.
    centers
        Bin centers, the abscissa of the interpolation.
    """
    y = np.asarray(y, dtype=float)
    n = np.asarray(n, dtype=float)
    has_events = in_band & (n > 0)
    measured = has_events & (y > 0)
    empty = in_band & ~has_events
    inverse = np.zeros_like(y)
    if measured.any():
        inverse[measured] = n[measured] / y[measured]
        discarded = has_events & ~measured
        next_to_discarded = np.interp(centers[empty], centers[has_events], discarded[has_events]) > 0
        interpolated = empty.copy()
        interpolated[empty] = ~next_to_discarded
        inverse[interpolated] = np.interp(centers[interpolated], centers[measured], inverse[measured])
    weights = np.full_like(y, np.inf)
    np.divide(1.0, inverse, out=weights, where=inverse > 0)
    return weights


def overlap_weights(event_weights):
    """Overlap weights ``(1 / c_i) / sum_k (1 / c_k)`` per bin; 0 where no segment covers the bin.

    Parameters
    ----------
    event_weights
        Array of shape ``(segments, bins)``, infinite where a segment does not cover the bin.
    """
    inverse = 1.0 / np.asarray(event_weights, dtype=float)
    total = inverse.sum(axis=0)
    weights = np.zeros_like(inverse)
    np.divide(inverse, total, out=weights, where=total > 0)
    return weights


def zero_weight_ranges(edges, masks, weights):
    """Ranges inside some band where no segment holds valid data, so every overlap weight is 0.

    Returns
    -------
        ``(lambda_min, lambda_max)`` per contiguous run of such bins of `edges`, in increasing
        wavelength.
    """
    lost = np.any(masks, axis=0) & ~(np.asarray(weights).sum(axis=0) > 0)
    changes = np.flatnonzero(np.diff(np.concatenate([[0], lost.astype(int), [0]])))
    return [(edges[start], edges[stop]) for start, stop in zip(changes[::2], changes[1::2], strict=True)]


class MergeDirectBeamSegments(PythonAlgorithm):
    def category(self):
        return "Reflectometry\\SNS"

    def name(self):
        return "MergeDirectBeamSegments"

    def version(self):
        return 1

    def summary(self):
        return "Merge corrected direct-beam segments into one event workspace, weighting overlaps by event weight."

    def seeAlso(self):  # noqa: N802
        return ["MergeRuns"]

    def PyInit(self):  # noqa: N802
        self.declareProperty(
            StringArrayProperty("InputWorkspaces", direction=Direction.Input),
            "Corrected direct-beam event workspaces in wavelength, one per segment; not modified",
        )
        self.declareProperty(
            "WeightingBinWidth",
            0.01,
            FloatBoundedValidator(lower=0.0, exclusive=True),
            "Width (Å) of the grid the overlap weights are evaluated on",
        )
        self.declareProperty(
            IEventWorkspaceProperty("OutputWorkspace", "", Direction.Output), "The composite event workspace"
        )
        self.declareProperty(
            MatrixWorkspaceProperty("OutputWeightsWorkspace", "", Direction.Output, PropertyMode.Optional),
            "One spectrum per segment holding its overlap weight on the weighting grid",
        )
        for name, doc in (
            ("WavelengthMin", "Each segment's lower band edge (Å), in input order"),
            ("WavelengthMax", "Each segment's upper band edge (Å), in input order"),
            ("GapMin", "Lower end (Å) of each uncovered range between the bands"),
            ("GapMax", "Upper end (Å) of each uncovered range between the bands"),
        ):
            self.declareProperty(FloatArrayProperty(name, [], direction=Direction.Output), doc)

    def validateInputs(self):  # noqa: N802
        names = self.getProperty("InputWorkspaces").value
        if len(names) == 0:
            return {"InputWorkspaces": "At least one segment is required"}

        issues = []
        reference = None
        for name in names:
            if not AnalysisDataService.doesExist(name):
                issues.append(f"{name}: no such workspace")
                continue
            ws = AnalysisDataService.retrieve(name)
            label = f"{name} (run {ws.getRunNumber()})"
            if not isinstance(ws, IEventWorkspace):
                issues.append(f"{label}: not an event workspace")
                continue
            if ws.getAxis(0).getUnit().unitID() != "Wavelength":
                issues.append(f"{label}: X axis must be in Wavelength, not {ws.getAxis(0).getUnit().unitID()!r}")
            if ws.getNumberEvents() == 0:
                issues.append(f"{label}: contains no events")
            spectra = [ws.getSpectrum(i) for i in range(ws.getNumberHistograms())]
            if any(spectrum.getEventType() == EventType.WEIGHTED_NOTIME for spectrum in spectra):
                issues.append(f"{label}: holds compressed events; the event weight needs one event per neutron")
            detectors = [tuple(spectrum.getDetectorIDs()) for spectrum in spectra]
            if reference is None:
                reference = (label, detectors)
            elif detectors != reference[1]:
                issues.append(f"{label}: spectra or detector IDs differ from {reference[0]}")
        return {"InputWorkspaces": "; ".join(issues)} if issues else {}

    def PyExec(self):  # noqa: N802
        names = list(self.getProperty("InputWorkspaces").value)
        step = self.getProperty("WeightingBinWidth").value
        segments = [AnalysisDataService.retrieve(name) for name in names]
        runs = [ws.getRunNumber() for ws in segments]

        # 1. Bands, from where each segment's events lie.
        bands = [(ws.getTofMin(), ws.getTofMax()) for ws in segments]
        self.setProperty("WavelengthMin", [low for low, _ in bands])
        self.setProperty("WavelengthMax", [high for _, high in bands])

        # 2. Coverage.
        gaps = coverage_gaps(bands)
        for gap_min, gap_max, left, right in gaps:
            logger.warning(
                f"No direct-beam segment covers {gap_min:.4f}-{gap_max:.4f} Å, between run {runs[left]} "
                f"and run {runs[right]}; the composite flux is zero there"
            )
        self.setProperty("GapMin", [gap[0] for gap in gaps])
        self.setProperty("GapMax", [gap[1] for gap in gaps])

        # 3. Grid.
        edges = weighting_grid(bands, step)
        centers = 0.5 * (edges[:-1] + edges[1:])

        # 4. Event weight, on each segment's ROI-summed spectrum.
        masks = [in_band_mask(edges, band) for band in bands]
        event_weights = []
        for ws, mask in zip(segments, masks, strict=True):
            summed, counts = summed_weights_and_counts(ws, edges)
            event_weights.append(event_weight(summed, counts, mask, centers))

        # 5. Overlap weights.
        weights = overlap_weights(event_weights)
        for low, high in zero_weight_ranges(edges, masks, weights):
            first = np.searchsorted(edges, low)
            covering = [str(run) for run, mask in zip(runs, masks, strict=True) if mask[first]]
            logger.warning(
                f"Every event of run {', '.join(covering)} in {low:.4f}-{high:.4f} Å has weight 0; "
                "the composite flux is zero there"
            )
        if not self.getProperty("OutputWeightsWorkspace").isDefault:
            weights_ws = CreateWorkspace(
                DataX=np.tile(edges, len(segments)),
                DataY=weights.ravel(),
                DataE=np.zeros(weights.size),
                NSpec=len(segments),
                UnitX="Wavelength",
                VerticalAxisUnit="Label",
                VerticalAxisValues=[str(run) for run in runs],
                StoreInADS=False,
            )
            self.setProperty("OutputWeightsWorkspace", weights_ws)

        if len(segments) == 1:
            merged = CloneWorkspace(InputWorkspace=segments[0], StoreInADS=False)
        else:
            merged = self._apply_and_merge(names, segments, masks, edges, weights)

        # 8. Axis: the single placeholder bin spanning every band.
        low, high = edges[0], edges[-1]
        output_ws = Rebin(InputWorkspace=merged, Params=[low, high - low, high], PreserveEvents=True, StoreInADS=False)
        self.setProperty("OutputWorkspace", output_ws)

    def _apply_and_merge(self, names, segments, masks, edges, weights):
        """Steps 6 and 7: scale each segment's events by its weights and append the segments."""
        registered = []
        try:
            merge_inputs = []
            for index, (name, ws, mask) in enumerate(zip(names, segments, masks, strict=True)):
                if np.all(weights[index][mask] == 1.0):
                    merge_inputs.append(name)
                    continue
                factors = CreateWorkspace(
                    DataX=edges,
                    DataY=weights[index],
                    DataE=np.zeros(len(edges) - 1),
                    UnitX="Wavelength",
                    StoreInADS=False,
                    EnableLogging=False,
                )
                aligned = RebinToWorkspace(
                    WorkspaceToRebin=ws, WorkspaceToMatch=factors, PreserveEvents=True, StoreInADS=False
                )
                weighted = Multiply(LHSWorkspace=aligned, RHSWorkspace=factors, StoreInADS=False)
                weighted_name = AnalysisDataService.unique_hidden_name()
                RebinToWorkspace(
                    WorkspaceToRebin=weighted, WorkspaceToMatch=ws, PreserveEvents=True, OutputWorkspace=weighted_name
                )
                registered.append(weighted_name)
                merge_inputs.append(weighted_name)
            return MergeRuns(InputWorkspaces=merge_inputs, StoreInADS=False)
        finally:
            for name in registered:
                if AnalysisDataService.doesExist(name):
                    DeleteWorkspace(Workspace=name, EnableLogging=False)


AlgorithmFactory.subscribe(MergeDirectBeamSegments)
