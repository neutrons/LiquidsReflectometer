import logging

from mantid.api import IEventWorkspace

from lr_reduction.exceptions import ConfigError, OperationInputError, SampleLogsError
from lr_reduction.models.config import PeakFitConfig, ReductionConfig, RunParameters
from lr_reduction.models.results import CompositeDirectBeam
from lr_reduction.models.run_data import RunData
from lr_reduction.operations.interfaces import OperationInterface
from lr_reduction.properties.cd_attenuation import compute_cd_thickness, get_attenuator_log
from lr_reduction.utils.workspace import workspace_exists, workspace_handle

logger = logging.getLogger(__name__)


def cd_order(runs: list[RunData], flip_atten: bool) -> list[tuple[RunData, float]]:
    """The direct-beam runs in ascending order of the Cd thickness in their beam.

    The composite uses this one order everywhere: the thickest-Cd run, last, gets the
    lowest-wavelength dead-time threshold, and the segments are merged and logged in this order.

    Parameters
    ----------
    runs
        The composite's direct-beam runs.
    flip_atten
        Whether the runs' attenuator logs were recorded with inverted polarity.

    Returns
    -------
    list[tuple[RunData, float]]
        Each run with its Cd thickness in centimeters, thinnest first. Runs of equal
        thickness keep their input order.

    Raises
    ------
    SampleLogsError
        A run has no usable attenuator log.
    """
    thicknesses = [compute_cd_thickness(get_attenuator_log(run.workspace), flip_atten=flip_atten) for run in runs]
    return sorted(zip(runs, thicknesses, strict=True), key=lambda pair: pair[1])


def direct_beam_roi(params: RunParameters) -> tuple[tuple[int, int], tuple[int, int]]:
    """The detector box each direct-beam run is cut to before it is corrected.

    The box is the direct beam's own characterization window, with the peak range widened by
    `peak_fit.pad` on both sides so the beam-peak fit sees the peak's shoulders. A direct beam
    has no background subtraction, so no background strips are included.

    Parameters
    ----------
    params
        The direct beam's effective parameters, its own overrides merged with the defaults.

    Returns
    -------
    tuple[tuple[int, int], tuple[int, int]]
        ``(x_range, y_range)``, each an inclusive ``(min, max)`` pixel range: x from `low_res`,
        y from the padded `peak`.

    Raises
    ------
    ConfigError
        `peak` or `low_res` is unset, or the padded peak range starts below pixel 0.
    """
    if params.peak is None:
        raise ConfigError("the direct beam has no `peak` range, in its own block or in `defaults`")
    if params.low_res is None:
        raise ConfigError("the direct beam has no `low_res` range, in its own block or in `defaults`")
    pad = (params.peak_fit or PeakFitConfig()).pad
    y_min = params.peak.min - pad
    if y_min < 0:
        raise ConfigError(
            f"the direct beam's `peak` range starts at pixel {params.peak.min}, so padding it by "
            f"`peak_fit.pad` = {pad} would start the box at pixel {y_min}"
        )
    return (params.low_res.min, params.low_res.max), (y_min, params.peak.max + pad)


def _run_label(run: RunData) -> str:
    return "run " + "+".join(str(run_number) for run_number in run.run_numbers)


def _workspace_issues(name: str) -> list[str]:
    """Why the workspace `name` cannot enter the composition, or nothing when it can."""
    if not workspace_exists(name):
        return [f"no workspace named {name!r} in the analysis data service"]
    workspace = workspace_handle(name)
    if not isinstance(workspace, IEventWorkspace):
        return [f"workspace {name!r} is not an event workspace"]
    return []


class DirectBeamCompositionOperation(OperationInterface[list[RunData], ReductionConfig, CompositeDirectBeam]):
    """Combine multiple direct beam runs into a single composite direct beam result.

    Attributes
    ----------
    data : list[RunData]
        The direct-beam runs to be combined.
    config : ReductionConfig
        The full reduction configuration. Besides the direct beam's own block, the composition
        reads the correction settings under `corrections` and the instrument settings under
        `defaults`.
    direct_beam_name : str
        Name of the composite to build, selecting `config.direct_beams[direct_beam_name]`.
    """

    def __init__(self, data: list[RunData], config: ReductionConfig, direct_beam_name: str):
        super().__init__(data, config)
        self.direct_beam_name = direct_beam_name

    def validate_input(self) -> None:
        """Check what the operation itself reads before the first algorithm runs.

        That is the configuration, the run set, each run's workspace and its attenuator log, which
        orders the runs. Each algorithm validates its own inputs, the unit and the logs it reads
        included, when the composition calls it. Every problem found here is reported together,
        each naming its run.

        Raises
        ------
        OperationInputError
            `direct_beam_name` is not configured; `data` is empty or its runs differ from the
            configured ones; or a run's workspace is missing or not an event workspace, or its
            attenuator log cannot be read.
        """
        if self.direct_beam_name not in self.config.direct_beams:
            raise OperationInputError(
                f"No direct beam named {self.direct_beam_name!r} in the configuration; "
                f"defined: {sorted(self.config.direct_beams)}"
            )
        context = f"Cannot build the composite direct beam {self.direct_beam_name!r}"
        if not self.data:
            raise OperationInputError(f"{context}: no direct-beam runs were given")

        direct_beam = self.config.direct_beams[self.direct_beam_name]
        given = sorted(run_number for run in self.data for run_number in run.run_numbers)
        if given != sorted(direct_beam.direct_beam_run_numbers):
            raise OperationInputError(
                f"{context}: runs {given} were given, but the configuration lists "
                f"{sorted(direct_beam.direct_beam_run_numbers)}"
            )

        issues = []
        for run in self.data:
            label = _run_label(run)
            run_issues = _workspace_issues(run.workspace)
            if run_issues:
                issues.extend(f"{label}: {issue}" for issue in run_issues)
                continue
            try:
                get_attenuator_log(run.workspace)
            except SampleLogsError as error:
                issues.append(f"{label}: cannot read the attenuator log: {error}")
            if self.config.corrections.dead_time.apply and run.error_events_workspace is None:
                logger.warning(
                    f"{label} has no error-events workspace, so its dead-time correction is computed without "
                    "the error events"
                )
        if issues:
            raise OperationInputError(f"{context}:\n" + "\n".join(f"- {issue}" for issue in issues))

    def process(self) -> CompositeDirectBeam:
        """
        Perform the direct beam composition operations and return the result.
        """
        # Implement processing logic here
        ...

    def cleanup(self) -> None:
        """
        Perform any necessary cleanup after processing.
        """
        # Implement cleanup logic here
        pass
