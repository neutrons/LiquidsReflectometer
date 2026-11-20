import logging
from secrets import token_hex

import numpy as np
from mantid import simpleapi as api
from mantid.api import IEventWorkspace

# Importing the algorithms package registers every reduction algorithm with `mantid.simpleapi`.
import lr_reduction.operations.algorithms  # noqa: F401
from lr_reduction import app_config
from lr_reduction.exceptions import ConfigError, OperationInputError, SampleLogsError
from lr_reduction.models.config import DirectBeamConfig, PeakFitConfig, ReductionConfig, RunParameters
from lr_reduction.models.results import CompositeDirectBeam, DirectBeamSegment
from lr_reduction.models.run_data import RunData
from lr_reduction.operations.algorithms.slit_scan_correction import FACTOR_LOG
from lr_reduction.operations.interfaces import OperationInterface
from lr_reduction.properties.cd_attenuation import compute_cd_thickness, get_attenuator_log
from lr_reduction.types import ID, MantidWorkspaceName
from lr_reduction.utils.sample_logs import SampleLogs
from lr_reduction.utils.workspace import workspace_exists, workspace_handle

logger = logging.getLogger(__name__)

#: The log of the detector arm angle, averaged per run into `DirectBeamSegment.tthd`.
TTHD_LOG = "tthd"


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


def _new_name(identifier: str, token: str, role: str) -> MantidWorkspaceName:
    """An unused workspace name, ``<identifier>__<token>_<role>[_<n>]``.

    Stands in for the workspace store's naming until the store exists.
    """
    base = f"{identifier}__{token}_{role}"
    name, counter = base, 2
    while workspace_exists(name):
        name, counter = f"{base}_{counter}", counter + 1
    return name


def _release(*names: MantidWorkspaceName) -> None:
    for name in names:
        if workspace_exists(name):
            api.DeleteWorkspace(Workspace=name, EnableLogging=False)


def _mean_and_spread(values: list[float]) -> tuple[float, float]:
    """The mean and population standard deviation of `values`."""
    return float(np.mean(values)), float(np.std(values))


def _run_number(run: RunData) -> ID:
    return run.run_numbers[0]


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
        self._token = token_hex(4)
        #: Every workspace name the operation issued, and those of them that are results.
        self._created: list[MantidWorkspaceName] = []
        self._kept: set[MantidWorkspaceName] = set()

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
        """Correct every run, merge the segments into the composite, and record what was done.

        The runs are processed in ascending Cd thickness, each on its own working copy cut to the
        direct beam's box, so the loaded workspaces are never modified. Only the composite, and the
        merge weights in diagnostic mode, are kept; `cleanup` deletes everything else the
        operation created, and everything when a step fails.

        Raises
        ------
        NotImplementedError
            A step whose algorithm does not exist yet was reached.
        """
        direct_beam = self.config.direct_beams[self.direct_beam_name]
        params = self.config.effective(direct_beam)
        roi = direct_beam_roi(params)
        peak_fit = params.peak_fit or PeakFitConfig()
        ordered = cd_order(self.data, direct_beam.flip_atten)

        identifier = f"REF_L_{_run_number(ordered[0][0])}"
        composite = self._new_name(identifier, "db")
        weights = self._new_name(identifier, "db_weights") if app_config["diagnostic_plots"] else None

        working_copies = []
        per_run = []
        for index, (run, cd_thickness) in enumerate(ordered):
            work = self._new_name(f"REF_L_{_run_number(run)}", "work")
            working_copies.append(work)
            lowest_wavelength = index == len(ordered) - 1
            per_run.append(self._correct_run(run, work, cd_thickness, lowest_wavelength, roi, peak_fit, direct_beam))

        merge_outputs = {"OutputWeightsWorkspace": weights} if weights is not None else {}
        merge = api.MergeDirectBeamSegments(InputWorkspaces=working_copies, OutputWorkspace=composite, **merge_outputs)
        bands = zip(merge.WavelengthMin, merge.WavelengthMax, strict=True)
        segments = tuple(
            DirectBeamSegment(**fields, wavelength_band=(float(low), float(high)))
            for fields, (low, high) in zip(per_run, bands, strict=True)
        )
        gaps = tuple((float(low), float(high)) for low, high in zip(merge.GapMin, merge.GapMax, strict=True))
        self._write_metadata(composite, segments)

        result = CompositeDirectBeam(
            workspace=composite,
            run_numbers=tuple(segment.run_number for segment in segments),
            segments=segments,
            gaps=gaps,
            weights_workspace=weights,
        )
        self._kept = {name for name in (composite, weights) if name is not None}
        return result

    def _new_name(self, identifier: str, role: str) -> MantidWorkspaceName:
        """A new workspace name for this operation, recorded so `cleanup` can delete it."""
        name = _new_name(identifier, self._token, role)
        self._created.append(name)
        return name

    def _correct_run(
        self,
        run: RunData,
        work: MantidWorkspaceName,
        cd_thickness: float,
        lowest_wavelength: bool,
        roi: tuple[tuple[int, int], tuple[int, int]],
        peak_fit: PeakFitConfig,
        direct_beam: DirectBeamConfig,
    ) -> dict:
        """Run the per-run steps on a working copy of `run` named `work`.

        Returns the run's `DirectBeamSegment` fields, all but the wavelength band, which the merge
        reports.
        """
        curve = self._new_name(f"REF_L_{_run_number(run)}", "dead_time_curve")
        dead_time = self.config.corrections.dead_time
        tof_trim = direct_beam.tof_trim
        emission_time = self.config.corrections.emission_time

        # 1. The dead-time curve, from the full detector. The trim reads it even when the events
        #    are not corrected.
        if dead_time.apply or tof_trim.apply:
            error_events = {}
            if run.error_events_workspace is not None:
                error_events["InputErrorEventsWorkspace"] = run.error_events_workspace
            api.SingleReadoutDeadTimeCorrection(
                InputWorkspace=run.workspace,
                DeadTime=dead_time.value,
                TOFStep=dead_time.tof_step,
                Paralyzable=dead_time.paralyzable,
                UseDeadTimeThreshold=False,
                OutputWorkspace=curve,
                **error_events,
            )

        # 2. The direct beam's box; this creates the working copy.
        x_range, y_range = roi
        api.ExtractDetectorROI(
            InputWorkspace=run.workspace,
            XRange=list(x_range),
            YRange=list(y_range),
            OutputWorkspace=work,
        )
        # 3. The instrument settings of the sequence.
        self._apply_instrument_settings(work)
        # 4. Proton-charge normalization, so the trim's intensity floor is per µA·h.
        api.NormaliseByCurrent(InputWorkspace=work, OutputWorkspace=work)
        # 5. Dead time, in TOF.
        if dead_time.apply:
            api.ApplyDeadTimeCorrection(InputWorkspace=work, DeadTimeCorrectionWorkspace=curve, OutputWorkspace=work)
        # 6. TOF trim, in TOF; the thickest-Cd run covers the shortest wavelengths.
        trim_fields = {
            "dead_time_threshold": None,
            "tof_range": None,
            "chopper_cut": None,
            "dead_time_cut": None,
            "intensity_cut": None,
        }
        if tof_trim.apply:
            threshold = (
                tof_trim.dead_time_threshold_lowest_wavelength if lowest_wavelength else tof_trim.dead_time_threshold
            )
            trim = api.TrimDirectBeamTof(
                InputWorkspace=work,
                DeadTimeCorrectionWorkspace=curve,
                DeadTimeThreshold=threshold,
                IntensityThreshold=tof_trim.intensity_threshold,
                IntensityTOFStep=tof_trim.intensity_tof_step,
                ChopperCutCoefficients=list(tof_trim.chopper_cut.coefficients),
                ChopperCutOffset=tof_trim.chopper_cut.offset_ms,
                OutputWorkspace=work,
            )
            trim_fields = {
                "dead_time_threshold": threshold,
                "tof_range": (trim.TOFMin, trim.TOFMax),
                "chopper_cut": trim.ChopperCut,
                "dead_time_cut": trim.DeadTimeCut,
                "intensity_cut": trim.IntensityCut,
            }
        # 7. TOF -> wavelength; T0 and T1 come from the run's own chopper logs unless configured.
        coefficients = {}
        if emission_time.coefficients is not None:
            coefficients = {"T0": emission_time.coefficients[0], "T1": emission_time.coefficients[1]}
        api.ConvertTofToWavelength(
            InputWorkspace=work,
            OutputWorkspace=work,
            ApplyEmissionTimeCorrection=emission_time.apply,
            **coefficients,
        )
        # 8. Cd transmission, in wavelength.
        api.CdAttenuationCorrection(InputWorkspace=work, FlipAtten=direct_beam.flip_atten, OutputWorkspace=work)
        # 9. Slit-scan rescale.
        api.SlitScanCorrection(InputWorkspace=work, OutputWorkspace=work)
        # 10. The beam-peak centroid.
        db_pixel, db_pixel_err = self._fit_beam_peak(work, peak_fit)
        # 11. The detector arm angle.
        logs = SampleLogs(work)
        return {
            "run_number": _run_number(run),
            "cd_thickness_um": cd_thickness * 1e4,
            "slit_scan_factor": float(logs[FACTOR_LOG]),
            **trim_fields,
            "db_pixel": db_pixel,
            "db_pixel_err": db_pixel_err,
            "tthd": logs.time_average(TTHD_LOG),
        }

    def _apply_instrument_settings(self, work: MantidWorkspaceName) -> None:
        """Step 3: put the sequence's instrument settings in force on `work`."""
        raise NotImplementedError(
            "ApplyInstrumentSettings does not exist yet, so the composite direct beam cannot be built"
        )

    def _fit_beam_peak(self, work: MantidWorkspaceName, peak_fit: PeakFitConfig) -> tuple[float, float]:
        """Step 10: the fitted beam-peak centroid of `work` and its uncertainty, in pixels."""
        raise NotImplementedError("FitBeamPeak does not exist yet, so the composite direct beam cannot be built")

    @staticmethod
    def _write_metadata(composite: MantidWorkspaceName, segments: tuple[DirectBeamSegment, ...]) -> None:
        """Record the beam position, the arm angle, and the composite's constituents as sample logs."""
        logs = SampleLogs(composite)
        db_pixel, db_pixel_err = _mean_and_spread([segment.db_pixel for segment in segments])
        tthd, tthd_err = _mean_and_spread([segment.tthd for segment in segments])
        logs.insert("db_pixel", db_pixel)
        logs.insert("db_pixel_err", db_pixel_err)
        logs.insert("tthd", tthd)
        logs.insert("tthd_err", tthd_err)
        logs.insert("Cd_attenuator", [segment.cd_thickness_um for segment in segments], unit="um")
        logs.insert("db_runs", [segment.run_number for segment in segments])

    def cleanup(self) -> None:
        """Delete every workspace the operation created except the results it returned.

        After a successful `process` that is the working copies and the dead-time curves; after a
        failure, it is everything, the composite included. Safe to call more than once.
        """
        _release(*(name for name in self._created if name not in self._kept))
        self._created = [name for name in self._created if name in self._kept]
