"""
Slit-scan intensity rescaling for direct-beam runs.

A slit-scan direct beam is measured in ``N`` smaller sub-openings, so its intensity is low by
``N`` per slit dimension. Every event's weight is multiplied by ``N²``.
"""

import math

from mantid.api import AlgorithmFactory, IEventWorkspaceProperty, PythonAlgorithm
from mantid.kernel import Direction
from mantid.simpleapi import CloneWorkspace, Scale, logger

from lr_reduction.exceptions import AmbiguousLogError, LrValidationError
from lr_reduction.utils.sample_logs import SampleLogs

#: Sample log holding the number of slit-scan sub-openings, ``N``.
SCALE_MULTIPLIER_LOG = "BL4B:CS:Autoreduce:ScaleMultiplier"

#: Sample log recording the factor ``N²`` applied, for provenance.
FACTOR_LOG = "slit_scan_factor"


def get_scale_multiplier(workspace) -> float:
    """``N`` from the run's `SCALE_MULTIPLIER_LOG`, or 1 for a run that predates the log.

    Raises
    ------
    LrValidationError
        The log changes during the run, or ``N`` is zero, negative or undefined.
    """
    run_number = workspace.getRunNumber()
    logs = SampleLogs(workspace)
    if SCALE_MULTIPLIER_LOG not in logs:
        logger.notice(f"Run {run_number} has no {SCALE_MULTIPLIER_LOG} log; treating it as N = 1")
        return 1.0
    try:
        n = logs[SCALE_MULTIPLIER_LOG]
    except AmbiguousLogError as exc:
        raise LrValidationError(f"Run {run_number}: {exc}") from exc
    if isinstance(n, bool) or not isinstance(n, (int, float)) or not math.isfinite(n) or n <= 0:
        raise LrValidationError(f"Run {run_number}: {SCALE_MULTIPLIER_LOG} must be a positive number, got {n!r}")
    return float(n)


class SlitScanCorrection(PythonAlgorithm):
    def category(self):
        return "Reflectometry\\SNS"

    def name(self):
        return "SlitScanCorrection"

    def version(self):
        return 1

    def summary(self):
        return (
            "Rescale a slit-scan direct-beam run by N squared, where N is the number of sub-openings "
            f"logged as {SCALE_MULTIPLIER_LOG}."
        )

    def PyInit(self):  # noqa: N802
        self.declareProperty(
            IEventWorkspaceProperty("InputWorkspace", "", Direction.Input),
            f"Direct-beam event workspace carrying a {SCALE_MULTIPLIER_LOG} log; never modified",
        )
        self.declareProperty(
            IEventWorkspaceProperty("OutputWorkspace", "", Direction.Output),
            f"Event workspace with weights multiplied by N squared, carrying a {FACTOR_LOG} log",
        )

    def PyExec(self):  # noqa: N802
        input_workspace = self.getProperty("InputWorkspace").value
        output_name = self.getPropertyValue("OutputWorkspace")

        factor = get_scale_multiplier(input_workspace) ** 2
        if factor == 1.0:
            # Scale by 1 would still turn plain events into weighted ones, costing memory for nothing.
            output = CloneWorkspace(InputWorkspace=input_workspace, OutputWorkspace=output_name)
        else:
            logger.notice(f"Slit-scan factor for run {input_workspace.getRunNumber()}: {factor:g}")
            output = Scale(
                InputWorkspace=input_workspace, Factor=factor, Operation="Multiply", OutputWorkspace=output_name
            )
        SampleLogs(output).insert(FACTOR_LOG, factor)
        self.setProperty("OutputWorkspace", output)


AlgorithmFactory.subscribe(SlitScanCorrection)
