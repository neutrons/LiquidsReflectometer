"""
Apply dead time correction to the data, as calculated by the SingleReadoutDeadTimeCorrection algorithm,
in order to correct for the effects of dead time in the measurement system.
"""

import numpy as np
from mantid.api import AlgorithmFactory, IEventWorkspaceProperty, MatrixWorkspaceProperty, PythonAlgorithm
from mantid.kernel import Direction
from mantid.simpleapi import Multiply, RebinToWorkspace, logger


class ApplyDeadTimeCorrection(PythonAlgorithm):
    def category(self):
        return "Reflectometry\\SNS"

    def name(self):
        return "ApplyDeadTimeCorrection"

    def version(self):
        return 1

    def summary(self):
        return "Apply dead time correction to the data."

    def PyInit(self):  # noqa: N802
        self.declareProperty(
            IEventWorkspaceProperty("InputWorkspace", "", Direction.Input),
            "Input workspace to which the dead time correction will be applied",
        )
        self.declareProperty(
            MatrixWorkspaceProperty("DeadTimeCorrectionWorkspace", "", Direction.Input),
            "Workspace containing the dead time correction factors from the SingleReadoutDeadTimeCorrection algorithm",
        )
        self.declareProperty(
            IEventWorkspaceProperty("OutputWorkspace", "", Direction.Output), "Corrected event workspace"
        )

    def validateInputs(self):  # noqa: N802
        issues = {}
        input_ws = self.getProperty("InputWorkspace").value
        curve_ws = self.getProperty("DeadTimeCorrectionWorkspace").value

        if input_ws.getAxis(0).getUnit().unitID() != "TOF":
            issues["InputWorkspace"] = "The X axis must be in TOF"

        curve_issues = []
        if curve_ws.getAxis(0).getUnit().unitID() != "TOF":
            curve_issues.append("the X axis must be in TOF")
        if curve_ws.getNumberHistograms() != 1:
            curve_issues.append(f"it must have exactly one spectrum, found {curve_ws.getNumberHistograms()}")
        else:
            tof_edges = curve_ws.readX(0)
            if "InputWorkspace" not in issues and input_ws.getNumberEvents() > 0:
                tof_min, tof_max = input_ws.getTofMin(), input_ws.getTofMax()
                if tof_edges[0] > tof_min or tof_edges[-1] < tof_max:
                    curve_issues.append(
                        f"its TOF range [{tof_edges[0]}, {tof_edges[-1]}] does not cover "
                        f"the events' TOF range [{tof_min}, {tof_max}]"
                    )
            factors = curve_ws.readY(0)
            if not np.all(np.isfinite(factors)):
                curve_issues.append("it contains undefined (NaN or infinite) values")
            if np.any(factors < 0):
                curve_issues.append("it contains negative values")
        if curve_issues:
            issues["DeadTimeCorrectionWorkspace"] = "Invalid dead time correction curve: " + "; ".join(curve_issues)

        return issues

    def PyExec(self):  # noqa: N802
        """
        Apply dead time correction to the input workspace using the provided dead time correction workspace.

        1. Line up the event workspace's placeholder axis with the curve's bins (`RebinToWorkspace`).
        2. Multiply the event workspace by the rebinned dead time correction factors (`Multiply`).
        3. Put the event workspace's original placeholder axis back (`RebinToWorkspace`, keeping all events),
           so that the output workspace has the same axis shape as the input.
        """
        input_ws = self.getProperty("InputWorkspace").value
        dead_time_correction_ws = self.getProperty("DeadTimeCorrectionWorkspace").value

        aligned_ws = RebinToWorkspace(
            WorkspaceToRebin=input_ws,
            WorkspaceToMatch=dead_time_correction_ws,
            PreserveEvents=True,
            StoreInADS=False,
        )
        corrected_ws = Multiply(LHSWorkspace=aligned_ws, RHSWorkspace=dead_time_correction_ws, StoreInADS=False)
        output_ws = RebinToWorkspace(
            WorkspaceToRebin=corrected_ws,
            WorkspaceToMatch=input_ws,
            PreserveEvents=True,
            StoreInADS=False,
        )

        self.setProperty("OutputWorkspace", output_ws)
        output_ws_name = self.getPropertyValue("OutputWorkspace")
        logger.information(f"Applied dead time correction to {input_ws.name()} and saved as {output_ws_name}.")


AlgorithmFactory.subscribe(ApplyDeadTimeCorrection)
