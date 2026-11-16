"""
Beam-peak localisation: the fitted y-pixel centroid of a beam on the BL4B area detector.

The ROI-extracted events are collapsed into a transient y-profile, summed over x and over the
whole X range, and the profile is fitted with Mantid's ``Fit``: a super-Gaussian (or Gaussian)
peak plus a linear (or flat, or no) background. The super-Gaussian is not in Mantid's function
library, so this module also registers it, as `SuperGaussian`.
"""

import numpy as np
from mantid.api import (
    AlgorithmFactory,
    CompositeFunction,
    FunctionFactory,
    IEventWorkspace,
    IFunction1D,
    MatrixWorkspaceProperty,
    PropertyMode,
    PythonAlgorithm,
)
from mantid.kernel import Direction, StringListValidator
from mantid.simpleapi import logger

from lr_reduction.operations.algorithms.roi_extraction import detector_shape

#: Mantid peak function for each `PeakFunction` choice.
PEAK_FUNCTIONS = {"supergauss": "SuperGaussian", "gauss": "Gaussian"}

#: Mantid background function for each `BackgroundFunction` choice, and its number of parameters.
BACKGROUND_FUNCTIONS = {"linear": ("LinearBackground", 2), "flat": ("FlatBackground", 1), "none": (None, 0)}

#: Number of parameters of each peak function.
PEAK_PARAMETERS = {"supergauss": 4, "gauss": 3}

EXPONENT_BOUNDS = (2.0, 10.0)

#: Penalty factor that makes Mantid's (penalty-based) boundary constraints effectively hard.
CONSTRAINT_PENALTY = 1e10

#: Fit statuses of the Levenberg-Marquardt minimizer that mean it converged.
CONVERGED = ("success", "Changes in function value are too small", "Changes in parameter value are too small")


class SuperGaussian(IFunction1D):
    """``Height * exp(-|(x - PeakCentre) / Sigma| ** Exponent)``, with ``2 <= Exponent <= 10``.

    No factor of two in the exponent: at ``Exponent = 2`` its ``Sigma`` is √2 times the
    ``Sigma`` of Mantid's ``Gaussian``.
    """

    def category(self):
        return "Peak"

    def init(self):
        self.declareParameter("Height", 1.0, "Peak height")
        self.declareParameter("PeakCentre", 0.0, "Peak centre")
        self.declareParameter("Sigma", 1.0, "Peak width")
        self.declareParameter("Exponent", 2.0, "Shape exponent; 2 is a Gaussian, larger is flatter-topped")
        self.addConstraints(f"{EXPONENT_BOUNDS[0]} < Exponent < {EXPONENT_BOUNDS[1]}")
        self.addConstraints("0 < Sigma")
        # Constraints are penalties; at the default factor a measured peak fitted Exponent 1.70.
        self.setConstraintPenaltyFactor("Exponent", CONSTRAINT_PENALTY)
        self.setConstraintPenaltyFactor("Sigma", CONSTRAINT_PENALTY)

    def function1D(self, xvals):  # noqa: N802
        height = self.getParameterValue("Height")
        centre = self.getParameterValue("PeakCentre")
        sigma = self.getParameterValue("Sigma")
        exponent = self.getParameterValue("Exponent")
        return height * np.exp(-(np.abs((xvals - centre) / sigma) ** exponent))


def y_profile(detector_ids, counts, errors, n_y):
    """Sum per-pixel counts over x into a y-profile, errors in quadrature.

    Pixels map to ``y = detector_id mod n_y``. A y pixel with no counts gets an error of 1, so it
    neither dominates nor drops out of a weighted fit.

    Returns
    -------
    ``(y_pixels, counts, errors)``, with `y_pixels` sorted.
    """
    y_pixels, inverse = np.unique(np.asarray(detector_ids) % n_y, return_inverse=True)
    summed = np.bincount(inverse, weights=counts, minlength=y_pixels.size)
    summed_errors = np.sqrt(np.bincount(inverse, weights=np.square(errors), minlength=y_pixels.size))
    summed_errors[summed == 0] = 1.0
    return y_pixels.astype(float), summed, summed_errors


def starting_function(y_pixels, counts, peak_function, background_function):
    """Mantid function string for the fit, with its starting values estimated from the profile."""
    floor = counts.min()
    centre = y_pixels[np.argmax(counts)]
    above = counts - floor
    sigma = 0.5
    if above.any():
        mean = np.average(y_pixels, weights=above)
        sigma = max(np.sqrt(np.average((y_pixels - mean) ** 2, weights=above)), 0.5)

    peak = f"name={PEAK_FUNCTIONS[peak_function]},Height={counts.max() - floor},PeakCentre={centre},Sigma={sigma}"
    if peak_function == "supergauss":
        peak += ",Exponent=3"
    peak += f",constraints=({y_pixels[0]}<PeakCentre<{y_pixels[-1]},0<Height)"

    background, _ = BACKGROUND_FUNCTIONS[background_function]
    if background == "LinearBackground":
        return f"{peak};name=LinearBackground,A0={floor},A1=0"
    if background == "FlatBackground":
        return f"{peak};name=FlatBackground,A0={floor}"
    return peak


class FitBeamPeak(PythonAlgorithm):
    def category(self):
        return "Reflectometry\\SNS"

    def name(self):
        return "FitBeamPeak"

    def version(self):
        return 1

    def summary(self):
        return (
            "Fit the y-profile of an ROI-extracted beam, summed over x and the whole X range, with a peak "
            "plus background, and return the fitted centroid in pixels."
        )

    def PyInit(self):  # noqa: N802
        self.declareProperty(
            MatrixWorkspaceProperty("InputWorkspace", "", Direction.Input),
            "ROI-extracted workspace with one detector pixel per spectrum; every y pixel in it is fitted",
        )
        self.declareProperty("PeakFunction", "supergauss", StringListValidator(list(PEAK_FUNCTIONS)), "Peak model")
        self.declareProperty(
            "BackgroundFunction", "linear", StringListValidator(list(BACKGROUND_FUNCTIONS)), "Background model"
        )
        for name, doc in (
            ("PeakCentre", "Fitted centroid, in y pixels"),
            ("PeakCentreError", "Fit uncertainty of PeakCentre"),
            ("PeakWidth", "Fitted Sigma of the peak function"),
            ("ChiSquaredOverDoF", "Goodness of fit"),
        ):
            self.declareProperty(name, 0.0, doc, direction=Direction.Output)
        self.declareProperty("Function", "", "The fitted function with its final parameters", Direction.Output)
        self.declareProperty(
            MatrixWorkspaceProperty("OutputWorkspace", "", Direction.Output, PropertyMode.Optional),
            "The fitted profile: spectra Data, Calc and Diff, X in y pixels",
        )

    def validateInputs(self):  # noqa: N802
        issues = {}
        input_ws = self.getProperty("InputWorkspace").value
        if input_ws is None:
            return issues

        _, n_y = detector_shape(input_ws)
        if n_y is None:
            issues["InputWorkspace"] = "The instrument must define number-of-y-pixels"
            return issues

        detector_ids = self._detector_ids(input_ws)
        if detector_ids is None:
            issues["InputWorkspace"] = "Every spectrum must hold exactly one detector pixel"
            return issues

        peak = self.getProperty("PeakFunction").value
        background = self.getProperty("BackgroundFunction").value
        n_parameters = PEAK_PARAMETERS[peak] + BACKGROUND_FUNCTIONS[background][1]
        n_pixels = np.unique(detector_ids % n_y).size
        if n_pixels < n_parameters + 2:
            issues["InputWorkspace"] = (
                f"InputWorkspace holds {n_pixels} y pixels, but fitting {peak} + {background} "
                f"needs at least {n_parameters + 2}"
            )
        return issues

    def PyExec(self):  # noqa: N802
        input_ws = self.getProperty("InputWorkspace").value
        peak = self.getProperty("PeakFunction").value
        background = self.getProperty("BackgroundFunction").value
        _, n_y = detector_shape(input_ws)
        detector_ids = self._detector_ids(input_ws)
        no_counts = f"no counts in the {np.unique(detector_ids % n_y).size} y pixels of InputWorkspace"
        # Integration cannot bin an event workspace without events.
        if isinstance(input_ws, IEventWorkspace) and input_ws.getNumberEvents() == 0:
            raise RuntimeError(no_counts)

        # Transient histogram: one bin per pixel over the whole X range; the input is untouched.
        integration = self.createChildAlgorithm("Integration")
        integration.setProperty("InputWorkspace", input_ws)
        # Without an explicit range, Integration drops some events of an event workspace (Mantid 6.16).
        integration.setProperty("RangeLower", -1e300)
        integration.setProperty("RangeUpper", 1e300)
        integration.execute()
        integrated = integration.getProperty("OutputWorkspace").value

        y_pixels, counts, errors = y_profile(
            detector_ids, integrated.extractY()[:, 0], integrated.extractE()[:, 0], n_y
        )
        if not np.any(counts):
            raise RuntimeError(no_counts)

        create = self.createChildAlgorithm("CreateWorkspace")
        create.setProperty("DataX", y_pixels)
        create.setProperty("DataY", counts)
        create.setProperty("DataE", errors)
        create.setProperty("UnitX", "Empty")
        create.execute()
        profile = create.getProperty("OutputWorkspace").value

        start = starting_function(y_pixels, counts, peak, background)
        fit = self.createChildAlgorithm("Fit")
        fit.setProperty("Function", start)
        fit.setProperty("InputWorkspace", profile)
        fit.setProperty("Minimizer", "Levenberg-Marquardt")
        fit.setProperty("CreateOutput", True)
        fit.setProperty("Output", "beam_peak_fit")
        fit.execute()
        status = fit.getProperty("OutputStatus").value
        if status not in CONVERGED:
            raise RuntimeError(f"The beam-peak fit failed: {status}. Starting function: {start}")

        function = fit.getProperty("Function").value
        peak_function = function.getFunction(0) if isinstance(function, CompositeFunction) else function
        centre = peak_function.getParameterValue("PeakCentre")
        self.setProperty("PeakCentre", centre)
        self.setProperty("PeakCentreError", peak_function.getError("PeakCentre"))
        self.setProperty("PeakWidth", peak_function.getParameterValue("Sigma"))
        self.setProperty("ChiSquaredOverDoF", fit.getProperty("OutputChi2overDoF").value)
        self.setProperty("Function", str(function))
        if not self.getProperty("OutputWorkspace").isDefault:
            self.setProperty("OutputWorkspace", fit.getProperty("OutputWorkspace").value)

        if peak == "supergauss":
            exponent = peak_function.getParameterValue("Exponent")
            # Mantid constraints are penalties, so a fit can end slightly outside the bounds.
            if not EXPONENT_BOUNDS[0] + 1e-3 < exponent < EXPONENT_BOUNDS[1] - 1e-3:
                logger.warning(f"FitBeamPeak: the fitted Exponent {exponent} is at a bound of {EXPONENT_BOUNDS}")
        if centre - y_pixels[0] < 1 or y_pixels[-1] - centre < 1:
            logger.warning(
                f"FitBeamPeak: the fitted centre {centre} is within a pixel of the window "
                f"[{y_pixels[0]}, {y_pixels[-1]}]"
            )
        logger.information(f"FitBeamPeak: centre {centre}, function {function}")

    @staticmethod
    def _detector_ids(workspace):
        """The single detector ID of each spectrum, or None when a spectrum holds other than one."""
        ids = [workspace.getSpectrum(index).getDetectorIDs() for index in range(workspace.getNumberHistograms())]
        if any(len(spectrum_ids) != 1 for spectrum_ids in ids):
            return None
        return np.array([spectrum_ids[0] for spectrum_ids in ids])


AlgorithmFactory.subscribe(FitBeamPeak)
FunctionFactory.subscribe(SuperGaussian)
