"""
Direct-beam TOF trimming: remove the unreliable time-of-flight regions of a direct-beam run
before the composite direct beam is assembled.

Three judgments each resolve to a bound on the TOF axis:

- frame overlap: a lower bound from the chopper-2 phase, `chopper_cut`;
- dead-time reliability: a lower bound from the ``DTC(TOF)`` curve, `dead_time_cut`;
- intensity floor: an upper bound from a transient histogram of the run's events, `intensity_cut`.

The events are cropped once, to ``[max(chopper_cut, dead_time_cut), intensity_cut)``.
"""

import numpy as np
from mantid.api import AlgorithmFactory, IEventWorkspaceProperty, MatrixWorkspaceProperty, PythonAlgorithm
from mantid.kernel import (
    CompositeValidator,
    Direction,
    FloatArrayLengthValidator,
    FloatArrayProperty,
    FloatBoundedValidator,
    FloatMandatoryValidator,
    Property,
)
from mantid.simpleapi import CropWorkspace, Rebin, SumSpectra, logger

from lr_reduction.operations.algorithms.apply_dead_time_correction import dead_time_curve_issues
from lr_reduction.utils.sample_logs import SampleLogs

#: Sample log holding the chopper-2 phase, averaged over the run for the frame-overlap cut.
CHOPPER2_PHASE_LOG = "BL4B:Chop:Skf2:PhaseAccuracy"


def chopper_cut(phase, coefficients, offset_ms):
    """Lower TOF bound (µs) below which chopper-2 frame overlap contaminates the run.

    ``t_cut [ms] = (a1 * phase - a0) / 1000 + offset``.

    Parameters
    ----------
    phase
        Time-averaged chopper-2 phase.
    coefficients
        ``(a1, a0)`` of the linear fit of the cut time to the phase.
    offset_ms
        Offset added to the fitted cut time, in ms.
    """
    a1, a0 = coefficients
    cut_ms = (a1 * phase - a0) / 1000.0 + offset_ms
    return cut_ms * 1000.0


def dead_time_cut(edges, factors, threshold):
    """Lower TOF bound (µs) after which the dead time correction stays at or below `threshold`.

    The upper edge of the last bin whose factor exceeds `threshold`, or None when no bin does.

    Parameters
    ----------
    edges
        Bin edges of the ``DTC(TOF)`` curve.
    factors
        The curve's correction factors, one per bin.
    threshold
        Largest correction factor considered reliable.
    """
    above = np.flatnonzero(np.asarray(factors) > threshold)
    if above.size == 0:
        return None
    return float(edges[above[-1] + 1])


def intensity_cut(edges, counts, threshold):
    """Upper TOF bound (µs) after which the intensity stays below `threshold`.

    The upper edge of the last bin at or above `threshold`, or None when no bin reaches it.
    Bins below `threshold` before that one are kept: only the tail is trimmed.

    Parameters
    ----------
    edges
        Bin edges of the intensity histogram.
    counts
        The histogram's summed event weights, one per bin.
    threshold
        Smallest intensity per bin considered reliable.
    """
    reaching = np.flatnonzero(np.asarray(counts) >= threshold)
    if reaching.size == 0:
        return None
    return float(edges[reaching[-1] + 1])


class TrimDirectBeamTof(PythonAlgorithm):
    def category(self):
        return "Reflectometry\\SNS"

    def name(self):
        return "TrimDirectBeamTof"

    def version(self):
        return 1

    def summary(self):
        return (
            "Remove the TOF regions of a direct-beam run affected by chopper-2 frame overlap, "
            "unreliable dead time correction, or low intensity."
        )

    def PyInit(self):  # noqa: N802
        self.declareProperty(
            IEventWorkspaceProperty("InputWorkspace", "", Direction.Input),
            "ROI-extracted, dead-time-corrected direct-beam event workspace in TOF",
        )
        self.declareProperty(
            MatrixWorkspaceProperty("CorrectionWorkspace", "", Direction.Input),
            "Single-spectrum DTC(TOF) curve from SingleReadoutDeadTimeCorrection, computed from the full detector",
        )
        self.declareProperty(
            "DeadTimeThreshold",
            Property.EMPTY_DBL,
            CompositeValidator([FloatBoundedValidator(lower=1.0), FloatMandatoryValidator()]),
            "Largest dead time correction factor considered reliable",
        )
        self.declareProperty(
            "IntensityThreshold",
            Property.EMPTY_DBL,
            CompositeValidator([FloatBoundedValidator(lower=0.0), FloatMandatoryValidator()]),
            "Smallest intensity considered reliable, in the input's event weights per IntensityTOFStep bin",
        )
        self.declareProperty(
            "IntensityTOFStep",
            50.0,
            FloatBoundedValidator(lower=0.0, exclusive=True),
            "Bin width (µs) of the histogram the intensity floor is judged on",
        )
        self.declareProperty(
            FloatArrayProperty("ChopperCutCoefficients", [2.077, 16818.0], FloatArrayLengthValidator(2)),
            "(a1, a0) of the chopper-2 cut, t_cut [ms] = (a1 * phase - a0) / 1000 + ChopperCutOffset",
        )
        self.declareProperty("ChopperCutOffset", 1.0, "Offset (ms) added to the chopper-2 cut time")
        self.declareProperty(
            IEventWorkspaceProperty("OutputWorkspace", "", Direction.Output),
            "Event workspace holding only the events in the surviving TOF range",
        )
        for name, doc in (
            ("TOFMin", "Lower end (µs) of the applied TOF range"),
            ("TOFMax", "Upper end (µs, excluded) of the applied TOF range"),
            ("ChopperCut", "Frame-overlap lower bound (µs)"),
            ("DeadTimeCut", "Dead-time lower bound (µs); 0 when no bin exceeds the threshold"),
            ("IntensityCut", "Intensity-floor upper bound (µs)"),
        ):
            self.declareProperty(name, 0.0, doc, direction=Direction.Output)

    def validateInputs(self):  # noqa: N802
        issues = {}
        input_ws = self.getProperty("InputWorkspace").value
        curve_ws = self.getProperty("CorrectionWorkspace").value
        run = input_ws.getRunNumber()

        input_issues = []
        tof_range = None
        if input_ws.getAxis(0).getUnit().unitID() != "TOF":
            input_issues.append("the X axis must be in TOF")
        elif input_ws.getNumberEvents() == 0:
            input_issues.append("it contains no events")
        else:
            tof_range = (input_ws.getTofMin(), input_ws.getTofMax())

        logs = SampleLogs(input_ws)
        if CHOPPER2_PHASE_LOG not in logs:
            input_issues.append(f"frame overlap: it has no {CHOPPER2_PHASE_LOG} sample log")
        elif np.size(logs.property(CHOPPER2_PHASE_LOG).value) == 0:
            input_issues.append(f"frame overlap: its {CHOPPER2_PHASE_LOG} sample log has no entries")
        if input_issues:
            issues["InputWorkspace"] = f"Run {run}: invalid input workspace: " + "; ".join(input_issues)

        curve_issues = dead_time_curve_issues(curve_ws, tof_range)
        if curve_issues:
            issues["CorrectionWorkspace"] = f"Run {run}: invalid dead time correction curve: " + "; ".join(curve_issues)

        return issues

    def PyExec(self):  # noqa: N802
        """
        Trim the input to the TOF range all three judgments agree is reliable.

        1. Frame overlap: `chopper_cut` from the time-averaged chopper-2 phase.
        2. Dead-time reliability: `dead_time_cut` from the correction curve.
        3. Intensity floor: `intensity_cut` from a transient histogram of the ROI-summed events.
        4. Crop the events to ``[max(chopper_cut, dead_time_cut), intensity_cut)``.
        """
        input_ws = self.getProperty("InputWorkspace").value
        curve_ws = self.getProperty("CorrectionWorkspace").value
        run = input_ws.getRunNumber()

        phase = SampleLogs(input_ws).time_average(CHOPPER2_PHASE_LOG)
        t_chop = chopper_cut(
            phase, self.getProperty("ChopperCutCoefficients").value, self.getProperty("ChopperCutOffset").value
        )
        t_dtc = dead_time_cut(curve_ws.readX(0), curve_ws.readY(0), self.getProperty("DeadTimeThreshold").value)
        t_floor = self._intensity_cut(input_ws, run)

        tof_min = t_chop if t_dtc is None else max(t_chop, t_dtc)
        tof_max = t_floor
        bounds = f"chopper cut {t_chop}, dead time cut {t_dtc}, intensity cut {t_floor}"
        if tof_min >= tof_max:
            raise RuntimeError(f"Run {run}: no TOF range survives the trim ({bounds})")

        output_ws = self._crop(input_ws, tof_min, tof_max)
        if output_ws.getNumberEvents() == 0:
            raise RuntimeError(f"Run {run}: no events remain in the TOF range [{tof_min}, {tof_max}) ({bounds})")

        self.setProperty("OutputWorkspace", output_ws)
        self.setProperty("TOFMin", tof_min)
        self.setProperty("TOFMax", tof_max)
        self.setProperty("ChopperCut", t_chop)
        self.setProperty("DeadTimeCut", 0.0 if t_dtc is None else t_dtc)
        self.setProperty("IntensityCut", t_floor)
        logger.notice(f"Trimmed run {run} to TOF [{tof_min}, {tof_max}) µs ({bounds})")

    def _intensity_cut(self, input_ws, run):
        """Judge the intensity floor on a transient histogram of the summed events."""
        threshold = self.getProperty("IntensityThreshold").value
        step = self.getProperty("IntensityTOFStep").value

        summed = SumSpectra(InputWorkspace=input_ws, StoreInADS=False)
        # Rebin drops an event exactly on the last edge, so the grid ends just past the last event.
        histogram = Rebin(
            InputWorkspace=summed,
            Params=[input_ws.getTofMin(), step, float(np.nextafter(input_ws.getTofMax(), np.inf))],
            PreserveEvents=False,
            StoreInADS=False,
        )
        cut = intensity_cut(histogram.readX(0), histogram.readY(0), threshold)
        if cut is None:
            raise RuntimeError(
                f"Run {run}: intensity floor: no {step} µs bin reaches IntensityThreshold={threshold} "
                f"(histogram maximum {histogram.readY(0).max()})"
            )
        return cut

    @staticmethod
    def _crop(input_ws, tof_min, tof_max):
        """Keep the events with ``tof_min <= TOF < tof_max``, on the single-bin axis ``[tof_min, tof_max]``.

        CropWorkspace keeps events exactly at XMax, so the crop stops just below `tof_max`, and the
        rebin, which keeps every event, then restores the exact axis.
        """
        cropped = CropWorkspace(
            InputWorkspace=input_ws,
            XMin=tof_min,
            XMax=float(np.nextafter(tof_max, -np.inf)),
            StoreInADS=False,
        )
        return Rebin(
            InputWorkspace=cropped,
            Params=[tof_min, tof_max - tof_min, tof_max],
            PreserveEvents=True,
            StoreInADS=False,
        )


AlgorithmFactory.subscribe(TrimDirectBeamTof)
