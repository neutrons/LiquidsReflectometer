"""
Cd attenuator transmission correction for direct-beam runs.
"""

from mantid.api import AlgorithmFactory, IEventWorkspaceProperty, PythonAlgorithm
from mantid.kernel import Direction
from mantid.simpleapi import logger

from lr_reduction.exceptions import LogTypeError
from lr_reduction.properties.cd_attenuation import apply_correction, compute_cd_thickness
from lr_reduction.utils.sample_logs import SampleLogs


class CdAttenuationCorrection(PythonAlgorithm):
    def category(self):
        return "Reflectometry\\SNS"

    def name(self):
        return "CdAttenuationCorrection"

    def version(self):
        return 1

    def summary(self):
        return (
            "Apply transmission correction for cadmium attenuators along the beam path. "
            "The total Cd thickness is determined from the Atten log."
        )

    def PyInit(self):  # noqa: N802
        self.declareProperty(
            IEventWorkspaceProperty("InputWorkspace", "", Direction.Input),
            "Direct-beam event workspace in wavelength, carrying an ``Atten`` log of per-stage flags",
        )
        self.declareProperty(
            "FlipAtten",
            False,
            doc=(
                "If True, swap the 0/1 flags in the ``Atten`` log before use: it was recorded with "
                "inverted polarity, so 0 means the attenuator stage is in the beam"
            ),
        )
        self.declareProperty(
            IEventWorkspaceProperty("OutputWorkspace", "", Direction.Output), "Corrected event workspace"
        )

    def PyExec(self):  # noqa: N802
        input_workspace = self.getProperty("InputWorkspace").value
        flip_atten = self.getProperty("FlipAtten").value

        atten = SampleLogs(input_workspace).property("Atten")
        if hasattr(atten, "times"):
            raise LogTypeError("The Atten log must hold one flag per attenuator stage, not a time series")
        cd_thickness = compute_cd_thickness(atten.value, flip_atten=flip_atten)
        logger.notice(f"Cd attenuator thickness for run {input_workspace.getRunNumber()}: {cd_thickness:.5f} cm")

        corrected = apply_correction(input_workspace, cd_thickness, self.getPropertyValue("OutputWorkspace"))
        self.setProperty("OutputWorkspace", corrected)


AlgorithmFactory.subscribe(CdAttenuationCorrection)
