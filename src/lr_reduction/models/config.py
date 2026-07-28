"""Pydantic schema for the reduction workflow configuration.

Sequence-wide parameters and per-run parameters are modeled separately.
A single `RunParameters` bundle is used both as the sequence-wide `defaults` block and,
with every field optional, as per-run overrides on `DirectBeamConfig`/`ReflectedRunConfig`:
`None` means "inherit from defaults".

Composite direct beams live in their own named key space and are referenced by name from reflected runs.
"""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from lr_reduction.types import ID

###################
### Value types ###
###################


class PixelRange(BaseModel):
    """Inclusive detector-pixel range [min, max]."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    min: int = Field(ge=0)
    max: int = Field(ge=0)

    @model_validator(mode="after")
    def _ordered(self) -> PixelRange:
        if self.max < self.min:
            raise ValueError(f"pixel range max ({self.max}) < min ({self.min})")
        return self


class AcceptanceWindow(BaseModel):
    """Wavelength and/or TOF acceptance window.

    TODO: §2.3.1 (Q2):
        the requirements allow either wavelength or TOF acceptance, or a single
        standardized representation; both are modeled here until instrument scientists settle on one.
    """

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    wavelength: tuple[float, float] | None = None  # (lambda_min, lambda_max), Angstrom
    tof: tuple[float, float] | None = None  # (tof_min, tof_max), microseconds

    @model_validator(mode="after")
    def _ordered(self) -> AcceptanceWindow:
        if self.wavelength is not None:
            lo, hi = self.wavelength
            if lo <= 0:
                raise ValueError(f"wavelength min ({lo}) must be positive")
            if hi < lo:
                raise ValueError(f"wavelength max ({hi}) < min ({lo})")
        if self.tof is not None:
            lo, hi = self.tof
            if lo < 0:
                raise ValueError(f"tof min ({lo}) must be non-negative")
            if hi < lo:
                raise ValueError(f"tof max ({hi}) < min ({lo})")
        return self


class QBinning(BaseModel):
    """Final Q grid and Q-construction method (§3.3.6)."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    q_min: float = Field(0.001, gt=0)
    q_max: float = Field(0.5, gt=0)
    q_step: float = Field(0.005, gt=0, description="Logarithmic step size")
    method: Literal["constantTOF", "meanTheta", "constantQ"] = Field("constantQ")
    coverage_threshold: float = Field(1.0, ge=0, le=1, description="Q-line coverage threshold")

    @model_validator(mode="after")
    def _ordered(self) -> QBinning:
        if self.q_max <= self.q_min:
            raise ValueError(f"q_max ({self.q_max}) must exceed q_min ({self.q_min})")
        return self


class BackgroundConfig(BaseModel):
    """Background subtraction settings.

    TODO: §2.3.4 (Q3):
        keep both side-region and functional background as options; legacy also
        supports a second side-background window (`two_backgrounds`), carried here as `region2`.
    """

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    apply: bool = Field(True, description="Enable background subtraction")
    region: PixelRange | None = Field(None, description="Side-region background pixel range")
    region2: PixelRange | None = Field(None, description="Secondary side-region background pixel range")
    mode: Literal["side", "functional"] = Field("side", description="Background subtraction method")


class AngleConfig(BaseModel):
    """Scattering-angle determination.

    TODO: §2.3.3:
        the priority science-review item — old vs. new theta-determination methods
        differ and need to be confirmed with instrument scientists.
    """

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    use_calculated_theta: bool = Field(False, description="False: use motor angle; True: fit/calculate theta")
    theta_offset: float = Field(0.0, description="Offset angle in degrees")
    theta_offset_error: float | None = Field(None, description="Error in offset angle in degrees")


class DeadTimeConfig(BaseModel):
    """Dead-time correction parameters."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    apply: bool = Field(True, description="Enable dead-time correction")
    value: float = Field(4.2, description="Dead-time value")
    tof_step: float = Field(50, gt=0, description="Time-of-Flight step size")
    paralyzable: bool = Field(False)
    threshold: float | None = Field(None)


class EmissionTimeConfig(BaseModel):
    """Emission-time (moderator) correction parameters.

    TODO: §2.1.3 (Q1):
        whether emission-time coefficients should be caller-configurable at all,
        or always derived instrument-side, is still open.
    """

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    apply: bool = Field(True, description="Enable emission-time correction")
    coefficients: list[float] | None = Field(None)

    @model_validator(mode="after")
    def _coefficients_non_empty(self) -> EmissionTimeConfig:
        if self.coefficients is not None and not self.coefficients:
            raise ValueError("`coefficients` must not be empty when provided")
        return self


class Corrections(BaseModel):
    """Correction enable flags and parameters (§3.3.6), applied globally to every run.

    Background subtraction is not modeled here: unlike these corrections, it varies per run and
    is carried by `RunParameters.background` instead.
    """

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    dead_time: DeadTimeConfig = Field(default_factory=DeadTimeConfig)
    emission_time: EmissionTimeConfig = Field(default_factory=EmissionTimeConfig)
    # TODO: legacy `GravityDirection` supports UP/DOWN/OFF; §3.3.6 only asks for an enable flag.
    #       Confirm whether direction reversal is still operationally needed.
    gravity: bool = Field(False, description="Enable gravity correction")
    background: BackgroundConfig = Field(default_factory=BackgroundConfig)


class PeakFitConfig(BaseModel):
    """Specular-peak fitting parameters."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    type: Literal["gauss", "supergauss"] = Field("supergauss")
    pad: int = Field(1, ge=0)


class ResolutionConfig(BaseModel):
    """Detector resolution function."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    function: Literal["rectangular", "gaussian"] = Field("rectangular")
    sigma: float = Field(0.8, gt=0)


class GeometryOverride(BaseModel):
    """Instrument-geometry overrides; unset fields default to the IDF value."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    l1: float | None = Field(None, description="source-to-sample distance, meters")
    l2: float | None = Field(None, description="sample-to-detector distance, meters")


class RunFilter(BaseModel):
    """Time and/or log-value filter applied when loading a source run."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    start_time: float | None = Field(None, description="Seconds, relative to run start")
    stop_time: float | None = Field(None, description="Seconds, relative to run start")
    log_name: str | None = Field(None, description="Log name")
    log_min: float | None = Field(None, description="Minimum log value")
    log_max: float | None = Field(None, description="Maximum log value")

    @model_validator(mode="after")
    def _has_criterion(self) -> RunFilter:
        has_time = self.start_time is not None or self.stop_time is not None
        has_log = self.log_name is not None
        if not (has_time or has_log):
            raise ValueError("RunFilter requires at least one time or log-value criterion")
        return self

    @model_validator(mode="after")
    def _ordered(self) -> RunFilter:
        if self.start_time is not None and self.stop_time is not None and self.stop_time <= self.start_time:
            raise ValueError(f"stop_time ({self.stop_time}) must exceed start_time ({self.start_time})")
        if self.log_min is not None and self.log_max is not None and self.log_max <= self.log_min:
            raise ValueError(f"log_max ({self.log_max}) must exceed log_min ({self.log_min})")
        return self


##########################
### Per-run parameters ###
##########################


class RunParameters(BaseModel):
    """Per-run-capable scientific parameters.

    Used both as part of the global `defaults` block and, with every field optional, as per-run
    overrides on `DirectBeamConfig` (`None` means "inherit from defaults").
    """

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    peak: PixelRange | None = None
    background: BackgroundConfig | None = None
    low_res: PixelRange | None = None  # low-resolution-direction pixel range


class ReflectedRunParameters(RunParameters):
    """The per-run-capable scientific parameters that only make sense for a reflected run.

    Extends `RunParameters` with the fields a direct beam has no use for (acceptance window,
    Q binning, angle, stitch scale, curve trimming, …). Used both as the global `defaults` block
    and, with every field optional, as per-run overrides on `ReflectedRunConfig` (`None` means
    "inherit from defaults").
    """

    acceptance: AcceptanceWindow | None = Field(None)
    q_binning: QBinning | None = Field(None)
    angle: AngleConfig | None = Field(None)
    peak_fit: PeakFitConfig | None = Field(None)
    resolution: ResolutionConfig | None = Field(None)
    geometry: GeometryOverride | None = Field(None)
    scale_factor: float | None = Field(None, description="Manual stitch scale override")
    trim: tuple[int, int] | None = Field(None, description="(pre_cut, post_cut) points trimmed from curve ends")

    @model_validator(mode="after")
    def _trim_non_negative(self) -> ReflectedRunParameters:
        if self.trim is not None and (self.trim[0] < 0 or self.trim[1] < 0):
            raise ValueError(f"trim values must be non-negative, got {self.trim}")
        return self


#################
### Run types ###
#################


class DirectBeamConfig(RunParameters):
    """A named composite direct beam: the individual DB runs it is built from.

    Inherits `RunParameters` (rather than carrying only `db_runs`) because the direct beam has
    its own peak/background/low-res characterization window, distinct from any reflected run
    that references it, and should inherit from the same global `defaults`.
    """

    db_runs: list[ID] = Field(min_length=1)
    # TODO:
    #   Not explicitly specified for DB runs; included for symmetry with reflected runs
    #   Determine whether this is actually needed.
    filter: RunFilter | None = Field(None)


class ReflectedRunConfig(ReflectedRunParameters):
    """One reflected run keyed by sequence_number."""

    sequence_number: ID = Field(..., ge=1)
    direct_beam: str = Field(..., description="named DB reference (§3.3.4)")

    # Only one of these two fields is set: either a single run number, or a list of source runs to sum.
    run_number: ID | None = Field(None)  # TODO: could instead be resolved from NeXus logs by sequence_number
    source_runs: list[ID] | None = Field(None, min_length=1, description="List of reflected runs to sum")
    # NOTE: hard-fail on mismatched conditions is enforced at the summing operation, not at schema load

    filter: RunFilter | None = Field(None, description="§3.3.6.2")

    @model_validator(mode="after")
    def _one_source(self) -> ReflectedRunConfig:
        has_run_number = self.run_number is not None
        has_source_runs = self.source_runs is not None
        if has_run_number == has_source_runs:
            raise ValueError(
                f"sequence_number={self.sequence_number}: exactly one of `run_number` or `source_runs` must be set"
            )
        if has_source_runs and not self.source_runs:
            raise ValueError(f"sequence_number={self.sequence_number}: `source_runs` must not be empty")
        return self

    @property
    def resolved_source_runs(self) -> list[ID]:
        """The run(s) to load and (if more than one) sum for this reflected run."""
        if self.source_runs:
            return list(self.source_runs)
        assert self.run_number is not None
        return [self.run_number]


################################
### Sequence-wide containers ###
################################


class OutputConfig(BaseModel):
    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    stem: str = Field(
        "reduction_output",
        description="Output file stem; directory and extension are added automatically",
    )
    subtitle: str | None = Field(None, description="Optional subtitle for plots")
    eight_column: bool = Field(False, description="Use legacy 8-column output format (T, L, dT, dL)")


class DiagnosticsConfig(BaseModel):
    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    plots: bool = Field(True, description="Generate diagnostic plots")
    rq4: bool = Field(False, description="R*Q^4 view")


class AssemblyConfig(BaseModel):
    """Sequence-assembly parameters (§2.4) — apply to the combined curve, not one run.

    TODO:
        `autoscale` and `stitching_type` overlap conceptually (legacy `AutoScale` flag vs.
        `StitchingType.AUTOMATIC_AVERAGE`); reconcile with instrument scientists.
    """

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    critical_edge_normalization: bool = Field(False, description="Normalize at the critical edge")
    q_norm: float = Field(0.015, gt=0)
    autoscale: bool = Field(False, description="Automatically scale the data")
    stitching_type: Literal["manual", "automatic_average", "absolute_normalization"] = Field(
        "manual", description="Type of stitching to use"
    )
    scale_factor_q_range: tuple[float, float] | None = Field(None, description="Range of Q values for scaling factor")
    normalize_first_angle: bool = Field(False, description="Normalize the first angle")

    @model_validator(mode="after")
    def _ordered(self) -> AssemblyConfig:
        if self.scale_factor_q_range is not None:
            q_min, q_max = self.scale_factor_q_range
            if q_min <= 0:
                raise ValueError(f"scale_factor_q_range min ({q_min}) must be positive")
            if q_max <= q_min:
                raise ValueError(f"scale_factor_q_range max ({q_max}) must exceed min ({q_min})")
        return self


##################################
### Top level reduction config ###
##################################


class ReductionConfig(BaseModel):
    """Global configuration (§3.3.1)."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    # global settings
    instrument: str = Field("BL4B")
    output: OutputConfig = Field(default_factory=OutputConfig)
    diagnostics: DiagnosticsConfig = Field(default_factory=DiagnosticsConfig)
    assembly: AssemblyConfig = Field(default_factory=AssemblyConfig)
    corrections: Corrections = Field(default_factory=Corrections)

    # per-run defaults inherited by every run
    defaults: ReflectedRunParameters = Field(default_factory=ReflectedRunParameters)

    # the two per-run key spaces (§3.3.2/.3/.4)
    direct_beams: dict[str, DirectBeamConfig] = Field(default_factory=dict)
    runs: list[ReflectedRunConfig] = Field(default_factory=list)

    @model_validator(mode="after")
    def _referential_integrity(self) -> ReductionConfig:
        """Enforce §3.3.5: every reflected run's `direct_beam` reference must resolve."""
        names = set(self.direct_beams)
        for run in self.runs:
            if run.direct_beam not in names:
                raise ValueError(
                    f"run sequence_number={run.sequence_number} references unknown "
                    f"direct_beam '{run.direct_beam}'; defined: {sorted(names)}"
                )
        seqs = [r.sequence_number for r in self.runs]
        if len(seqs) != len(set(seqs)):
            raise ValueError(f"duplicate sequence_number values: {seqs}")
        return self

    def for_run(self, sequence_number: int) -> ReductionConfig:
        """Extract a valid one-run config (§3.3.8): same global settings, only the
        selected run and the direct beam it references. Re-validated like any other config,
        so no separate single-run format exists."""
        try:
            run = next(r for r in self.runs if r.sequence_number == sequence_number)
        except StopIteration:
            raise ValueError(f"no run with sequence_number={sequence_number}") from None
        data = self.model_dump()
        data["runs"] = [run.model_dump()]
        data["direct_beams"] = {run.direct_beam: self.direct_beams[run.direct_beam].model_dump()}
        return ReductionConfig(**data)

    def effective(self, run: ReflectedRunConfig | DirectBeamConfig) -> RunParameters:
        """The effective per-run parameters: `defaults` merged with `run`'s own overrides.

        A `DirectBeamConfig` only draws on the shared `RunParameters` fields; a
        `ReflectedRunConfig` draws on the full `ReflectedRunParameters` set.
        """
        param_cls = ReflectedRunParameters if isinstance(run, ReflectedRunConfig) else RunParameters
        fields = set(param_cls.model_fields)
        defaults = {k: v for k, v in self.defaults.model_dump(exclude_none=True).items() if k in fields}
        overrides = run.model_dump(include=fields, exclude_none=True)
        return param_cls(**_deep_merge(defaults, overrides))


def apply_overrides(config: ReductionConfig, overrides: dict) -> ReductionConfig:
    """Apply caller-supplied overrides on top of `config` for a single invocation (§3.3.7).

    `config` itself is never mutated; overrides take precedence over every value in `config`
    (file < overrides), applied as a deep merge over the full configuration dict, then the whole
    tree is re-validated — including cross-field checks like referential integrity (§3.3.5).

    For a single override, setting an attribute directly on a loaded config (or a nested value
    within it, e.g. `config.assembly.q_norm = 0.02`) is simpler and is still validated: every
    model here sets `validate_assignment=True`, so an assignment is checked against that field's
    own constraints (type, range, ordering) immediately. The one thing direct assignment does
    *not* do is re-run a *different* model's validators — e.g. setting
    `config.runs[0].direct_beam = "..."` re-checks `ReflectedRunConfig`'s own fields but does not
    re-run `ReductionConfig._referential_integrity`, since that validator lives on the parent
    model. Use `apply_overrides` (or `ReductionConfig.model_validate(config.model_dump())` after
    manual edits) when the override could affect cross-model checks.
    """
    merged = _deep_merge(config.model_dump(), overrides)
    return ReductionConfig(**merged)


def _deep_merge(base: dict, overrides: dict) -> dict:
    merged = dict(base)
    for key, value in overrides.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged
