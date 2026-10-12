"""
Detector region-of-interest extraction for the BL4B area detector.
"""

from mantid.api import AlgorithmFactory, IEventWorkspaceProperty, PythonAlgorithm
from mantid.kernel import (
    CompositeValidator,
    Direction,
    IntArrayBoundedValidator,
    IntArrayLengthValidator,
    IntArrayOrderedPairsValidator,
    IntArrayProperty,
)
from mantid.simpleapi import logger


def detector_shape(workspace):
    """Return ``(n_x, n_y)``, the detector's pixel counts from its instrument parameters.

    Returns None for a dimension the instrument does not define.
    """
    instrument = workspace.getInstrument()
    shape = []
    for parameter in ("number-of-x-pixels", "number-of-y-pixels"):
        values = instrument.getNumberParameter(parameter)
        shape.append(int(values[0]) if values else None)
    return tuple(shape)


def roi_workspace_indices(x_range, y_range, n_y):
    """Workspace indices of the inclusive pixel rectangle, using the BL4B rule ``index = x * n_y + y``."""
    x_min, x_max = x_range
    y_min, y_max = y_range
    return [x * n_y + y for x in range(x_min, x_max + 1) for y in range(y_min, y_max + 1)]


class ExtractDetectorROI(PythonAlgorithm):
    def category(self):
        return "Reflectometry\\SNS"

    def name(self):
        return "ExtractDetectorROI"

    def version(self):
        return 1

    def summary(self):
        return (
            "Keep only the detector pixels inside a rectangular region of interest (ROI), "
            "given as inclusive x and y pixel ranges. Events are kept as they are."
        )

    def PyInit(self):  # noqa: N802
        self.declareProperty(
            IEventWorkspaceProperty("InputWorkspace", "", Direction.Input),
            "Event workspace with one spectrum per detector pixel, in the order ``index = x * n_y + y``",
        )
        self.declareProperty(
            IEventWorkspaceProperty("OutputWorkspace", "", Direction.Output),
            "Event workspace holding only the spectra of the pixels inside the ROI",
        )
        for name in ("XRange", "YRange"):
            range_validator = CompositeValidator()
            range_validator.add(IntArrayLengthValidator(2))
            range_validator.add(IntArrayBoundedValidator(lower=0))
            range_validator.add(IntArrayOrderedPairsValidator())
            self.declareProperty(
                IntArrayProperty(name, [0, 0], range_validator, direction=Direction.Input),
                f"Inclusive [min, max] range of {name[0].lower()} pixels to keep",
            )

    def validateInputs(self):  # noqa: N802
        issues = {}
        input_workspace = self.getProperty("InputWorkspace").value
        if input_workspace is None:
            return issues

        n_x, n_y = detector_shape(input_workspace)
        n_spectra = input_workspace.getNumberHistograms()
        if n_x is None or n_y is None:
            issues["InputWorkspace"] = (
                f"The instrument must define number-of-x-pixels and number-of-y-pixels (found {n_x} x {n_y})"
            )
            return issues
        if n_x * n_y != n_spectra:
            issues["InputWorkspace"] = (
                f"The instrument's {n_x} x {n_y} detector has {n_x * n_y} pixels, "
                f"but the workspace has {n_spectra} spectra"
            )
            return issues

        for name, n_pixels in (("XRange", n_x), ("YRange", n_y)):
            low, high = self.getProperty(name).value
            if high >= n_pixels:
                issues[name] = (
                    f"{name} [{low}, {high}] runs off the {n_x} x {n_y} detector, "
                    f"whose {name[0].lower()} pixels are 0 to {n_pixels - 1}"
                )
        return issues

    def PyExec(self):  # noqa: N802
        input_workspace = self.getProperty("InputWorkspace").value
        x_range = self.getProperty("XRange").value
        y_range = self.getProperty("YRange").value

        # Get the detector shape and compute the workspace indices of the ROI pixels
        _, n_y = detector_shape(input_workspace)
        indices = roi_workspace_indices(x_range, y_range, n_y)
        logger.information(f"Extracting a ROI of {len(indices)} pixels: x {list(x_range)}, y {list(y_range)}")

        # Run ExtractSpectra as a child algorithm to keep only the spectra of the pixels in the ROI
        extract = self.createChildAlgorithm("ExtractSpectra")
        extract.setProperty("InputWorkspace", input_workspace)
        extract.setProperty("WorkspaceIndexList", indices)
        extract.execute()
        self.setProperty("OutputWorkspace", extract.getProperty("OutputWorkspace").value)


AlgorithmFactory.subscribe(ExtractDetectorROI)
