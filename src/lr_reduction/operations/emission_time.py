"""Where the moderator emission-time coefficients for a run come from (§2.1.3).

`ConvertTofToWavelength` takes the two coefficients as plain numbers; deciding them is the
caller's job and is done here, once, for both operations that reduce raw runs.
"""

from lr_reduction.exceptions import EmissionTimeCorrectionError, LogNotFoundError
from lr_reduction.models.config import EmissionTimeConfig
from lr_reduction.types import MantidWorkspace
from lr_reduction.utils.sample_logs import SampleLogs

#: Chopper PV flagging that the chopper-2 phasing includes the moderator emission delay.
CHOPPER_MODERATOR_LOG = "BL4B:Chop:Skf2:ChopperModerator"
#: Chopper PV holding the emission-time offset t0, in microseconds.
CHOPPER_OFFSET_LOG = "BL4B:Chop:Skf2:ChopperOffset"
#: Chopper PV holding the emission-time slope t1, in microseconds per Angstrom.
CHOPPER_MULTIPLIER_LOG = "BL4B:Chop:Skf2:ChopperMultiplier"


def resolve_emission_coefficients(config: EmissionTimeConfig, reflected_run: MantidWorkspace) -> tuple[float, float]:
    """The emission-time coefficients ``(t0, t1)`` to apply when converting a run to wavelength.

    In order of precedence:

    1. ``config.apply`` is False: ``(0, 0)``, which makes the conversion a plain TOF ->
       wavelength.
    2. ``config.coefficients`` is set: those values. The configuration wins over the logs.
    3. Otherwise the run's chopper logs: when the ``ChopperModerator`` flag is set, the
       ``ChopperOffset`` and ``ChopperMultiplier`` values; when it is clear, ``(0, 0)``.

    The reflected beam and every direct beam it is divided by must use the same coefficients
    (§2.1.4.1), so callers resolve them once, against the *reflected* run, and pass the same
    pair to every direct-beam run, even when a direct-beam run's own logs would give
    different numbers.

    Parameters
    ----------
    config
        The ``corrections.emission_time`` block of the reduction configuration.
    reflected_run
        The reflected-beam run workspace (object or name) whose chopper logs are read when the
        configuration supplies no coefficients.

    Returns
    -------
        ``(t0, t1)`` in microseconds and microseconds per Angstrom.

    Raises
    ------
    EmissionTimeCorrectionError
        The correction is enabled, no coefficients are configured, and the run lacks the
        chopper logs needed to decide them. There is no silent fallback to "no correction".
    """
    if not config.apply:
        return 0.0, 0.0
    if config.coefficients is not None:
        t0, t1 = config.coefficients
        return float(t0), float(t1)

    logs = SampleLogs(reflected_run)
    flag = _read_chopper_log(logs, CHOPPER_MODERATOR_LOG)
    if flag != 1:
        return 0.0, 0.0
    return float(_read_chopper_log(logs, CHOPPER_OFFSET_LOG)), float(_read_chopper_log(logs, CHOPPER_MULTIPLIER_LOG))


def _read_chopper_log(logs: SampleLogs, name: str) -> float:
    """A chopper PV's value, or an error explaining how to proceed without it.

    Chopper PVs are recorded as time series; the first entry is the setting the run was
    acquired with, the same value the production workflow reads.
    """
    try:
        return logs.single_value(name, operation=lambda values: values[0])
    except LogNotFoundError as exc:
        raise EmissionTimeCorrectionError(
            f"Cannot resolve emission-time coefficients: the run has no {name!r} log. "
            "Set `corrections.emission_time.coefficients: [t0, t1]` explicitly, or disable the "
            "correction with `corrections.emission_time.apply: false`"
        ) from exc
