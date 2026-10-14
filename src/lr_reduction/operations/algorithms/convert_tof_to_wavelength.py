"""
Emission-time-corrected conversion of an event workspace from time of flight to wavelength.

Requirements §2.1.3 (moderator emission-time TOF correction) and §2.1.4 (TOF -> wavelength),
folded into one step as specified in ``design/toftowavelength.pdf``.
"""

import math

from mantid.api import AlgorithmFactory, IEventWorkspaceProperty, PythonAlgorithm
from mantid.kernel import Direction, PhysicalConstants

from lr_reduction.utils.sample_logs import SampleLogs

#: Name of the instrument parameter `ModeratorTzero` reads the emission-delay formula from.
FORMULA_PARAMETER = "t0_formula"

#: Sample logs recording the coefficients that were applied to the output.
T0_LOG = "lr_emission_t0"
T1_LOG = "lr_emission_t1"

# `ModeratorTzero` stopping criteria. Each iteration shrinks the error against the closed form
# by roughly the ratio t1 / (m_n * L / h), about 1/1000 on REF_L. Measured on Mantid 6.16.1
# against the closed form: Niter=1 leaves ~1e-6 relative error (the two-pass result), Niter=2
# with tolTOF=0.1 us leaves ~1e-9. The unit tests pin these values at a 1e-8 tolerance.
_NITER = 2
_TOL_TOF = 0.1  # microseconds

# h^2 / (2 m_n) in meV * Angstrom^2, so that lambda = sqrt(_LAMBDA_SQ_ENERGY / E). Built from
# Mantid's own constants rather than a textbook value: ModeratorTzero converts each event's
# velocity to `incidentEnergy` with these constants, and a formula written with a different
# CODATA release would put a ~1e-9 relative floor under the result.
_LAMBDA_SQ_ENERGY = PhysicalConstants.h**2 / (2 * PhysicalConstants.NeutronMass * PhysicalConstants.meV) * 1e20


def emission_time_formula(t0: float, t1: float) -> str:
    """The linear emission-delay model ``t_emit(lambda) = t0 + t1 * lambda``, written in
    `ModeratorTzero`'s formula language.

    `ModeratorTzero` evaluates its formula as a function of ``incidentEnergy`` (meV), not of
    wavelength, so lambda is spelled ``sqrt(h^2 / (2 m_n) / incidentEnergy)``.

    Parameters
    ----------
    t0
        Emission-time offset, in microseconds.
    t1
        Emission-time slope, in microseconds per Angstrom.
    """
    return f"{t0:.17g} + {t1:.17g} * sqrt({_LAMBDA_SQ_ENERGY:.17g} / incidentEnergy)"


class ConvertTofToWavelength(PythonAlgorithm):
    """Convert an event workspace from TOF to wavelength, removing the moderator emission delay.

    The delay is modeled as linear in wavelength, ``t_emit(lambda) = t0 + t1 * lambda``. Since
    the delay depends on the wavelength and the wavelength on the corrected TOF, the two are
    solved together; for this model the exact solution is

        lambda = h * (TOF - t0) / (m_n * L + h * t1)

    with ``L`` the full flight path of the pixel that recorded the event.

    The arithmetic is delegated to two Mantid algorithms run on a copy of the input:
    `ModeratorTzero`, which reads the delay formula from an instrument parameter this
    algorithm writes onto the copy and iterates to the solution above; and `ConvertUnits`,
    which rescales each event by its own pixel's flight path and labels the axis Wavelength.
    Both read the geometry from the workspace's instrument, so any override applied earlier
    by the instrument-settings step reaches the conversion without extra code.

    With both coefficients zero the correction is skipped and the result is exactly a plain
    `ConvertUnits` to wavelength. The input is never modified.
    """

    def category(self):
        return "Reflectometry\\SNS"

    def name(self):
        return "ConvertTofToWavelength"

    def version(self):
        return 1

    def summary(self):
        return "Convert TOF to wavelength with the moderator emission-time correction folded in"

    def seeAlso(self):  # noqa: N802
        return ["ModeratorTzero", "ConvertUnits"]

    def PyInit(self):  # noqa: N802
        self.declareProperty(
            IEventWorkspaceProperty("InputWorkspace", "", Direction.Input),
            "Event workspace in TOF, carrying the instrument geometry in force. Not modified",
        )
        self.declareProperty("T0", 0.0, doc="Emission-time offset t0, in microseconds")
        self.declareProperty("T1", 0.0, doc="Emission-time slope t1, in microseconds per Angstrom")
        self.declareProperty(
            IEventWorkspaceProperty("OutputWorkspace", "", Direction.Output),
            "New event workspace in wavelength (Angstrom)",
        )

    def validateInputs(self):  # noqa: N802
        issues = {}

        workspace = self.getProperty("InputWorkspace").value
        if workspace is not None:
            unit = workspace.getAxis(0).getUnit().unitID()
            if unit != "TOF":
                issues["InputWorkspace"] = f"X axis must be in TOF, not {unit!r}"
            else:
                instrument = workspace.getInstrument()
                if instrument.getSource() is None or instrument.getSample() is None:
                    issues["InputWorkspace"] = (
                        "instrument has no source or sample position, so the flight path cannot be computed"
                    )

        for name in ("T0", "T1"):
            value = self.getProperty(name).value
            if not math.isfinite(value):
                issues[name] = f"must be a finite number, got {value}"

        return issues

    def PyExec(self):  # noqa: N802
        t0 = self.getProperty("T0").value
        t1 = self.getProperty("T1").value

        # Work on a copy: the formula parameter is written onto the workspace's instrument
        # parameter map, and ModeratorTzero shifts events in place.
        workspace = self._run_child("CloneWorkspace", InputWorkspace=self.getProperty("InputWorkspace").value)

        if t0 != 0.0 or t1 != 0.0:
            self._run_child(
                "SetInstrumentParameter",
                output_property=None,
                Workspace=workspace,
                ParameterName=FORMULA_PARAMETER,
                ParameterType="String",
                Value=emission_time_formula(t0, t1),
            )
            # EMode must be Elastic. ModeratorTzero's default (Indirect) assumes an
            # indirect-geometry spectrometer and silently produces wrong TOFs on REF_L.
            workspace = self._run_child(
                "ModeratorTzero",
                InputWorkspace=workspace,
                EMode="Elastic",
                Niter=_NITER,
                tolTOF=_TOL_TOF,
            )

        workspace = self._run_child("ConvertUnits", InputWorkspace=workspace, Target="Wavelength", EMode="Elastic")

        logs = SampleLogs(workspace)
        logs.insert(T0_LOG, t0, unit="microsecond")
        logs.insert(T1_LOG, t1, unit="microsecond/Angstrom")

        self.setProperty("OutputWorkspace", workspace)

    def _run_child(self, name, output_property="OutputWorkspace", **properties):
        """Run a child algorithm and return its output workspace (None if `output_property` is None).

        Child algorithms keep their intermediates out of the analysis data service, and their
        calls are nested under this algorithm in the output's history rather than appearing
        as siblings of it.
        """
        algorithm = self.createChildAlgorithm(name)
        for key, value in properties.items():
            algorithm.setProperty(key, value)
        if output_property is not None:
            algorithm.setPropertyValue(output_property, "_unused")
        algorithm.execute()
        return algorithm.getProperty(output_property).value if output_property is not None else None


AlgorithmFactory.subscribe(ConvertTofToWavelength)
