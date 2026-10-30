"""
Emission-time-corrected conversion of an event workspace from time of flight to wavelength.

Requirements §2.1.3 (moderator emission-time TOF correction) and §2.1.4 (TOF -> wavelength),
folded into one step.
"""

import math

from mantid.api import AlgorithmFactory, IEventWorkspaceProperty, PythonAlgorithm
from mantid.kernel import Direction, PhysicalConstants, Property
from mantid.simpleapi import CloneWorkspace, ConvertUnits, ModeratorTzero, SetInstrumentParameter

from lr_reduction.exceptions import SampleLogsError
from lr_reduction.utils.sample_logs import SampleLogs

#: Chopper PV holding the emission-time offset t0, in microseconds. Default for `T0`.
OFFSET_LOG = "BL4B:Chop:Skf2:ChopperOffset"
#: Chopper PV holding the emission-time slope t1, in microseconds per Angstrom. Default for `T1`.
MULTIPLIER_LOG = "BL4B:Chop:Skf2:ChopperMultiplier"

#: Name of the instrument parameter `ModeratorTzero` reads the emission-delay formula from.
FORMULA_PARAMETER = "t0_formula"

#: Sample logs recording the coefficients that were applied to the output.
T0_LOG = "lr_emission_t0"
T1_LOG = "lr_emission_t1"

# `ModeratorTzero` stopping criteria. Each iteration shrinks the error against the closed form
# by the ratio t1 / (m_n * L / h), only about 1/120 on REF_L for the production slope of
# ~29.5 us/A. Measured on Mantid 6.16.1 with t0 = 114 us, t1 = 29.5 us/A (max relative error):
# Niter=2, tolTOF=0.1 us gives 1.9e-6; tolTOF >= 0.01 us halts the iteration at 1.5e-8 however
# many passes are allowed; Niter=10, tolTOF=1e-6 us gives 9e-13. Iteration stops as soon as
# the TOF stops changing, so the generous cap costs little (~0.1 s per million events). The
# unit tests pin the result to the closed form at 1e-8. The error after `_NITER` passes is about
# ratio ** (_NITER + 1), so slopes with ratio above `_MAX_SLOPE_RATIO` are rejected (~1e-11 at the cap).
_NITER = 10
_TOL_TOF = 1e-6  # microseconds
_MAX_SLOPE_RATIO = 0.1

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


def _us_per_angstrom(path_m: float) -> float:
    """Flight time per Angstrom of wavelength, m_n * L / h, in microseconds per Angstrom."""
    return PhysicalConstants.NeutronMass * path_m / PhysicalConstants.h * 1e-4


class ConvertTofToWavelength(PythonAlgorithm):
    """Convert an event workspace from TOF to wavelength, removing the moderator emission delay.

    The delay is modeled as linear in wavelength, ``t_emit(lambda) = t0 + t1 * lambda``. Since
    the delay depends on the wavelength and the wavelength on the corrected TOF, the two are
    solved together; for this model the exact solution is

        lambda = h * (TOF - t0) / (m_n * L + h * t1)

    with ``L`` the full flight path of the pixel that recorded the event.

    The coefficients default to the run's own chopper logs, ``ChopperOffset`` for `T0` and
    ``ChopperMultiplier`` for `T1`; either property, when set, overrides its log. Setting both
    to zero, or setting `ApplyEmissionTimeCorrection` to false, turns the correction off, leaving a
    plain `ConvertUnits` to wavelength; the logs are then not needed. The ``ChopperModerator``
    flag is not consulted, matching the prototype workflow. The applied values are recorded on
    the output in the ``lr_emission_t0`` and ``lr_emission_t1`` logs.

    The arithmetic is delegated to two Mantid algorithms run on a copy of the input:
    `ModeratorTzero`, which reads the delay formula from an instrument parameter this
    algorithm writes onto the copy and iterates to the solution above; and `ConvertUnits`,
    which rescales each event by its own pixel's flight path and labels the axis Wavelength.
    Both read the geometry from the workspace's instrument, so any override applied earlier
    by the instrument-settings step reaches the conversion without extra code. The input is
    never modified.
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
        self.declareProperty(
            "ApplyEmissionTimeCorrection",
            True,
            doc="Remove the moderator emission delay. If false, T0 and T1 are ignored and a plain conversion is done",
        )
        self.declareProperty(
            "T0",
            Property.EMPTY_DBL,
            doc=f"Emission-time offset t0, in microseconds. Default: the run's ``{OFFSET_LOG}`` log",
        )
        self.declareProperty(
            "T1",
            Property.EMPTY_DBL,
            doc=f"Emission-time slope t1, in microseconds per Angstrom. Default: the run's ``{MULTIPLIER_LOG}`` log",
        )
        self.declareProperty(
            IEventWorkspaceProperty("OutputWorkspace", "", Direction.Output),
            "New event workspace in wavelength (Angstrom)",
        )

    def validateInputs(self):  # noqa: N802
        issues = {}

        workspace = self.getProperty("InputWorkspace").value
        if workspace is None:
            return issues

        unit = workspace.getAxis(0).getUnit().unitID()
        if unit != "TOF":
            issues["InputWorkspace"] = f"X axis must be in TOF, not {unit!r}"
        else:
            instrument = workspace.getInstrument()
            if instrument.getSource() is None or instrument.getSample() is None:
                issues["InputWorkspace"] = (
                    "instrument has no source or sample position, so the flight path cannot be computed"
                )

        if not self.getProperty("ApplyEmissionTimeCorrection").value:
            return issues

        coefficients = {}
        for property_name, log_name in (("T0", OFFSET_LOG), ("T1", MULTIPLIER_LOG)):
            try:
                coefficients[property_name] = self._coefficient(property_name, log_name)
            except SampleLogsError as exc:
                issues[property_name] = f"not set, and the run's log cannot supply it: {exc}"
            except ValueError as exc:
                issues[property_name] = str(exc)

        if "T1" in coefficients and "InputWorkspace" not in issues:
            # shortcut: bound uses L1 only, stricter than min(L1 + L2); use that if a real slope ever hits it.
            max_slope = _MAX_SLOPE_RATIO * _us_per_angstrom(workspace.spectrumInfo().l1())
            if abs(coefficients["T1"]) > max_slope:
                issues["T1"] = (
                    f"{coefficients['T1']} us/A is too large for ModeratorTzero to converge on this "
                    f"instrument; |T1| must not exceed {max_slope:.1f} us/A"
                )

        return issues

    def PyExec(self):  # noqa: N802
        if self.getProperty("ApplyEmissionTimeCorrection").value:
            t0 = self._coefficient("T0", OFFSET_LOG)
            t1 = self._coefficient("T1", MULTIPLIER_LOG)
        else:
            t0 = t1 = 0.0

        # Work on a copy: the formula parameter is written onto the workspace's instrument
        # parameter map, and ModeratorTzero shifts events in place.
        workspace = CloneWorkspace(InputWorkspace=self.getProperty("InputWorkspace").value, StoreInADS=False)

        if t0 != 0.0 or t1 != 0.0:
            SetInstrumentParameter(
                Workspace=workspace,
                ParameterName=FORMULA_PARAMETER,
                ParameterType="String",
                Value=emission_time_formula(t0, t1),
                StoreInADS=False,
            )
            # EMode must be Elastic. ModeratorTzero's default (Indirect) assumes an
            # indirect-geometry spectrometer and silently produces wrong TOFs on REF_L.
            workspace = ModeratorTzero(
                InputWorkspace=workspace,
                EMode="Elastic",
                Niter=_NITER,
                tolTOF=_TOL_TOF,
                StoreInADS=False,
            )

        workspace = ConvertUnits(InputWorkspace=workspace, Target="Wavelength", EMode="Elastic", StoreInADS=False)

        logs = SampleLogs(workspace)
        logs.insert(T0_LOG, t0, unit="microsecond")
        logs.insert(T1_LOG, t1, unit="microsecond/Angstrom")

        self.setProperty("OutputWorkspace", workspace)

    def _coefficient(self, property_name, log_name):
        """The value of `property_name` if it was set, otherwise the input run's `log_name` log.

        Raises
        ------
        SampleLogsError
            The property is unset and the log is missing, or varies across the run, so it has
            no single value to apply.
        ValueError
            The value, from either source, is not a finite number.
        """
        prop = self.getProperty(property_name)
        if prop.isDefault:
            value = SampleLogs(self.getProperty("InputWorkspace").value)[log_name]
            source = f"log {log_name!r}"
        else:
            value, source = prop.value, "the given value"
        # bool is an int, but a flag is not a coefficient.
        if isinstance(value, bool) or not isinstance(value, int | float) or not math.isfinite(value):
            raise ValueError(f"{source} is not a finite number: {value!r}")
        return float(value)


AlgorithmFactory.subscribe(ConvertTofToWavelength)
