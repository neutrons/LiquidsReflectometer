"""
Background subtraction under the specular peak of a reflected-beam histogram ``R(y, λ)``.

In each λ bin, a line ``a·y + c`` (or a constant ``c``) is fitted by unweighted least squares over
the rows of the background strips, subtracted from every peak row, and its variance, including the
slope-intercept covariance, is added to each peak row's variance.
"""

from collections import Counter

import numpy as np
from mantid.api import (
    CommonBinsValidator,
    HistogramValidator,
    IEventWorkspace,
    MatrixWorkspaceProperty,
    PropertyMode,
    PythonAlgorithm,
)
from mantid.kernel import (
    CompositeValidator,
    Direction,
    IntArrayBoundedValidator,
    IntArrayLengthValidator,
    IntArrayOrderedPairsValidator,
    IntArrayProperty,
    StringListValidator,
)

from lr_reduction.operations.algorithms.roi_extraction import detector_shape
from lr_reduction.utils.registration import register

MODELS = ("Linear", "Constant")
MIN_BACKGROUND_ROWS = {"Linear": 2, "Constant": 1}


def background_rows(peak_range, strips):
    """The sorted y rows in the union of the inclusive `strips`, minus the rows of the inclusive `peak_range`."""
    rows = set()
    for low, high in strips:
        rows.update(range(low, high + 1))
    rows.difference_update(range(peak_range[0], peak_range[1] + 1))
    return sorted(rows)


def fit_strip_background(y, r, sigma, model="Linear"):
    """Fit ``bkg = a·y + c`` independently in each bin by unweighted least squares.

    Parameters
    ----------
    y
        ``(n,)`` row coordinates of the background rows.
    r, sigma
        ``(n, n_bins)`` values and errors of the background rows.
    model
        ``"Linear"`` or ``"Constant"`` (``a = 0``).

    Returns
    -------
    tuple
        ``(a, c, cov)``: ``(n_bins,)`` slope and intercept, and the ``(n_bins, 2, 2)`` covariance of ``(a, c)``.
    """
    y = np.asarray(y, dtype=float)
    n = y.size
    if model == "Linear":
        s_y, s_yy = y.sum(), (y * y).sum()
        d = n * s_yy - s_y**2
        alpha = (n * y - s_y) / d
        gamma = (s_yy - s_y * y) / d
    elif model == "Constant":
        alpha = np.zeros(n)
        gamma = np.full(n, 1.0 / n)
    else:
        raise ValueError(f"model must be one of {MODELS}, got {model!r}")
    var = np.asarray(sigma, dtype=float) ** 2
    a = alpha @ r
    c = gamma @ r
    cov = np.empty((a.size, 2, 2))
    cov[:, 0, 0] = (alpha**2) @ var
    cov[:, 1, 1] = (gamma**2) @ var
    cov[:, 0, 1] = cov[:, 1, 0] = (alpha * gamma) @ var
    return a, c, cov


def evaluate_strip_background(fit, y):
    """Evaluate a `fit_strip_background` result on rows `y`; returns ``(bkg, var)``, each ``(len(y), n_bins)``."""
    a, c, cov = fit
    y = np.asarray(y, dtype=float)[:, np.newaxis]
    bkg = a * y + c
    var = y**2 * cov[:, 0, 0] + 2 * y * cov[:, 0, 1] + cov[:, 1, 1]
    return bkg, np.maximum(var, 0.0)  # cancellation can leave about -1e-17 where the variance is zero


def spectrum_rows(workspace, n_y):
    """The y row of each spectrum, or None for a spectrum whose detector IDs give no single row."""
    rows = []
    for index in range(workspace.getNumberHistograms()):
        ids = workspace.getSpectrum(index).getDetectorIDs()
        row_set = {detector_id % n_y for detector_id in ids}
        rows.append(row_set.pop() if len(row_set) == 1 else None)
    return rows


class SubtractStripBackground(PythonAlgorithm):
    def category(self):
        return "Reflectometry\\SNS"

    def name(self):
        return "SubtractStripBackground"

    def version(self):
        return 1

    def summary(self):
        return (
            "Fit a per-bin linear or constant background across detector rows over the background strips, "
            "and subtract it, with its uncertainty, from the rows of the peak."
        )

    def PyInit(self):  # noqa: N802
        self.declareProperty(
            MatrixWorkspaceProperty(
                "InputWorkspace",
                "",
                Direction.Input,
                CompositeValidator([HistogramValidator(), CommonBinsValidator()]),
            ),
            "Histogram with one spectrum per detector row (x summed), on common bins",
        )

        def range_validator():
            validator = CompositeValidator()
            validator.add(IntArrayBoundedValidator(lower=0))
            validator.add(IntArrayOrderedPairsValidator())
            return validator

        peak_validator = range_validator()
        peak_validator.add(IntArrayLengthValidator(2))
        self.declareProperty(
            IntArrayProperty("PeakRange", [], peak_validator, direction=Direction.Input),
            "Inclusive [y_lo, y_hi] of the peak",
        )
        self.declareProperty(
            IntArrayProperty("BackgroundRanges", [], range_validator(), direction=Direction.Input),
            "One or two inclusive [y_lo, y_hi] background strips, flattened",
        )
        self.declareProperty(
            "BackgroundModel", "Linear", StringListValidator(list(MODELS)), "a·y + c (Linear) or c (Constant)"
        )
        self.declareProperty(
            MatrixWorkspaceProperty("OutputWorkspace", "", Direction.Output),
            "The peak rows of the input, background subtracted",
        )
        self.declareProperty(
            MatrixWorkspaceProperty("BackgroundWorkspace", "", Direction.Output, PropertyMode.Optional),
            "The fitted background and its error on the peak rows",
        )
        self.declareProperty(
            IntArrayProperty("BackgroundPixels", [], direction=Direction.Output), "The y rows used by the fit"
        )

    def validateInputs(self):  # noqa: N802
        issues = {}
        workspace = self.getProperty("InputWorkspace").value
        if workspace is None:
            return issues
        if isinstance(workspace, IEventWorkspace):
            issues["InputWorkspace"] = "InputWorkspace is an event workspace; rebin it to a histogram first"
            return issues
        counts, errors = workspace.extractY(), workspace.extractE()
        if not (np.all(np.isfinite(counts)) and np.all(np.isfinite(errors)) and np.all(errors >= 0)):
            issues["InputWorkspace"] = "InputWorkspace must hold finite values and finite, non-negative errors"
            return issues

        _, n_y = detector_shape(workspace)
        if n_y is None:
            issues["InputWorkspace"] = "The instrument must define number-of-y-pixels"
            return issues
        rows = spectrum_rows(workspace, n_y)
        if None in rows:
            index = rows.index(None)
            issues["InputWorkspace"] = f"The detector IDs of workspace index {index} do not belong to a single y row"
            return issues
        shared = sorted(row for row, count in Counter(rows).items() if count > 1)
        if shared:
            issues["InputWorkspace"] = f"More than one spectrum maps to y rows {shared}"
            return issues

        peak = self.getProperty("PeakRange").value.tolist()
        flat = self.getProperty("BackgroundRanges").value.tolist()
        if len(flat) not in (2, 4):
            issues["BackgroundRanges"] = f"BackgroundRanges must hold one or two [min, max] pairs, got {flat}"
            return issues
        strips = [flat[i : i + 2] for i in range(0, len(flat), 2)]

        present = set(rows)
        for name, ranges in (("PeakRange", [peak]), ("BackgroundRanges", strips)):
            for low, high in ranges:
                missing = sorted(set(range(low, high + 1)) - present)
                if missing:
                    issues[name] = f"{name} [{low}, {high}] includes rows absent from the input: {missing}"
        if issues:
            return issues

        model = self.getProperty("BackgroundModel").value
        remaining = background_rows(peak, strips)
        if len(remaining) < MIN_BACKGROUND_ROWS[model]:
            issues["BackgroundRanges"] = (
                f"Strips {strips} minus peak {peak} leave rows {remaining}; "
                f"the {model} model needs at least {MIN_BACKGROUND_ROWS[model]}"
            )
        return issues

    def PyExec(self):  # noqa: N802
        workspace = self.getProperty("InputWorkspace").value
        peak = self.getProperty("PeakRange").value.tolist()
        flat = self.getProperty("BackgroundRanges").value.tolist()
        model = self.getProperty("BackgroundModel").value

        _, n_y = detector_shape(workspace)
        row_of_index = np.array(spectrum_rows(workspace, n_y))
        index_of_row = {row: index for index, row in enumerate(row_of_index)}
        bkg_rows = background_rows(peak, [flat[i : i + 2] for i in range(0, len(flat), 2)])
        bkg_indices = [index_of_row[row] for row in bkg_rows]
        peak_indices = [i for i, row in enumerate(row_of_index) if peak[0] <= row <= peak[1]]  # input order

        counts, errors = workspace.extractY(), workspace.extractE()
        fit = fit_strip_background(bkg_rows, counts[bkg_indices], errors[bkg_indices], model)
        bkg, var = evaluate_strip_background(fit, row_of_index[peak_indices])
        self.log().information(f"Background from rows {bkg_rows} ({model}) under peak rows {peak}")

        output = self._extract(workspace, peak_indices)
        for k, index in enumerate(peak_indices):
            output.setSharedY(k, counts[index] - bkg[k])
            output.setSharedE(k, np.sqrt(errors[index] ** 2 + var[k]))
        self.setProperty("OutputWorkspace", output)
        self.setProperty("BackgroundPixels", bkg_rows)

        if not self.getProperty("BackgroundWorkspace").isDefault:
            background = self._extract(workspace, peak_indices)
            for k in range(len(peak_indices)):
                background.setSharedY(k, bkg[k])
                background.setSharedE(k, np.sqrt(var[k]))
            self.setProperty("BackgroundWorkspace", background)

    def _extract(self, workspace, indices):
        extract = self.createChildAlgorithm("ExtractSpectra")
        extract.setProperty("InputWorkspace", workspace)
        extract.setProperty("WorkspaceIndexList", indices)
        extract.execute()
        return extract.getProperty("OutputWorkspace").value


register(SubtractStripBackground)
