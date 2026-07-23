"""Pydantic schema for the reduction workflow configuration (§3.3).

Global parameters and per-run parameters are modeled separately (§3.3.1). `RunParameters` holds
the ROI-selection fields (`peak`, `background`, `low_res`) shared by both direct beam and
reflected runs; `ReflectedRunParameters` extends it with the fields that only make sense for a
reflected run (acceptance window, Q binning, etc.); global-only settings such as `Corrections`
live directly on `ReductionConfig`. With every field optional, each
is used both as the global `defaults` block and as per-run overrides on
`DirectBeamConfig`/`ReflectedRunConfig` (§3.3.6) — `None` means "inherit from defaults". Composite
direct beams live in their own named key space (§3.3.3) and are referenced by name from reflected
runs (§3.3.4); referential integrity is enforced at load (§3.3.5).
"""

from __future__ import annotations

from typing import Literal, Optional, Union

from pydantic import BaseModel, ConfigDict, Field, model_validator

# --------------------------------------------------------------------------- value types ----


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

    OPEN §2.3.1 (Q2): the requirements allow either wavelength or TOF acceptance, or a single
    standardized representation; both are modeled here until instrument scientists settle on one.
    """

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    wavelength: Optional[tuple[float, float]] = None  # (lambda_min, lambda_max), Angstrom
    tof: Optional[tuple[float, float]] = None  # (tof_min, tof_max), microseconds

    @model_validator(mode="after")
    def _ordered(self) -> AcceptanceWindow:
        if self.wavelength is not None and self.wavelength[1] < self.wavelength[0]:
            raise ValueError(f"wavelength max ({self.wavelength[1]}) < min ({self.wavelength[0]})")
        if self.tof is not None and self.tof[1] < self.tof[0]:
            raise ValueError(f"tof max ({self.tof[1]}) < min ({self.tof[0]})")
        return self


class QBinning(BaseModel):
    """Final Q grid and Q-construction method (§3.3.6)."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    q_min: float = Field(0.001, gt=0)
    q_max: float = Field(0.5, gt=0)
    q_step: float = Field(0.005, gt=0)  # logarithmic step
    method: Literal["constantTOF", "meanTheta", "constantQ"] = "constantQ"
    coverage_threshold: float = Field(1.0, ge=0, le=1)  # Q-line coverage threshold

    @model_validator(mode="after")
    def _ordered(self) -> QBinning:
        if self.q_max <= self.q_min:
            raise ValueError(f"q_max ({self.q_max}) must exceed q_min ({self.q_min})")
        return self


class BackgroundConfig(BaseModel):
    """Background subtraction settings.

    OPEN §2.3.4 (Q3): keep both side-region and functional background as options; legacy also
    supports a second side-background window (`two_backgrounds`), carried here as `region2`.
    """

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    apply: bool = True
    region: Optional[PixelRange] = None
    region2: Optional[PixelRange] = None  # legacy `two_backgrounds` secondary side window
    mode: Literal["side", "functional"] = "side"


class AngleConfig(BaseModel):
    """Scattering-angle determination.

    OPEN §2.3.3: the priority science-review item — old vs. new theta-determination methods
    differ and need to be confirmed with instrument scientists.
    """

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    use_calculated_theta: bool = False  # False: use motor angle; True: fit/calculate theta
    theta_offset: float = 0.0  # degrees
    theta_offset_error: Optional[float] = None  # degrees


class DeadTimeConfig(BaseModel):
    """Dead-time correction parameters."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    apply: bool = True
    value: float = 4.2
    tof_step: float = Field(50, gt=0)
    paralyzable: bool = True
    threshold: Optional[float] = None


class EmissionTimeConfig(BaseModel):
    """Emission-time (moderator) correction parameters.

    OPEN §2.1.3 (Q1): whether emission-time coefficients should be caller-configurable at all,
    or always derived instrument-side, is still open.
    """

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    apply: bool = True
    coefficients: Optional[list[float]] = None


class Corrections(BaseModel):
    """Correction enable flags and parameters (§3.3.6), applied globally to every run.

    Background subtraction is not modeled here: unlike these corrections, it varies per run and
    is carried by `RunParameters.background` instead.
    """

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    dead_time: DeadTimeConfig = Field(default_factory=DeadTimeConfig)
    emission_time: EmissionTimeConfig = Field(default_factory=EmissionTimeConfig)
    # OPEN: legacy `GravityDirection` supports UP/DOWN/OFF; §3.3.6 only asks for an enable flag.
    # Confirm whether direction reversal is still operationally needed.
    gravity: bool = True


class PeakFitConfig(BaseModel):
    """Specular-peak fitting parameters."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    type: Literal["gauss", "supergauss"] = "supergauss"
    pad: int = Field(1, ge=0)


class ResolutionConfig(BaseModel):
    """Detector resolution function."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    function: Literal["rectangular", "gaussian"] = "rectangular"
    sigma: float = Field(0.8, gt=0)


class GeometryOverride(BaseModel):
    """Instrument-geometry overrides; unset fields default to the IDF value (§2.1.4.3)."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    l1: Optional[float] = None  # source-to-sample distance, meters
    l2: Optional[float] = None  # sample-to-detector distance, meters


class RunFilter(BaseModel):
    """Time and/or log-value filter applied when loading a source run (§3.1.4)."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    start_time: Optional[float] = None  # seconds, relative to run start
    stop_time: Optional[float] = None
    log_name: Optional[str] = None
    log_min: Optional[float] = None
    log_max: Optional[float] = None

    @model_validator(mode="after")
    def _has_criterion(self) -> RunFilter:
        has_time = self.start_time is not None or self.stop_time is not None
        has_log = self.log_name is not None
        if not (has_time or has_log):
            raise ValueError("RunFilter requires at least one time or log-value criterion")
        return self


# ------------------------------------------------------------------- per-run parameter bundle ----


class RunParameters(BaseModel):
    """ROI-selection parameters shared by both direct beam and reflected runs (§3.3.6).

    Used both as part of the global `defaults` block and, with every field optional, as per-run
    overrides on `DirectBeamConfig` (`None` means "inherit from defaults").
    """

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    peak: Optional[PixelRange] = None
    background: Optional[BackgroundConfig] = None
    low_res: Optional[PixelRange] = None  # low-resolution-direction pixel range


class ReflectedRunParameters(RunParameters):
    """The §3.3.6 per-run-capable scientific parameters that only make sense for a reflected run.

    Extends `RunParameters` with the fields a direct beam has no use for (acceptance window,
    Q binning, angle, stitch scale, curve trimming, …). Used both as the global `defaults` block
    and, with every field optional, as per-run overrides on `ReflectedRunConfig` (`None` means
    "inherit from defaults").
    """

    acceptance: Optional[AcceptanceWindow] = None
    q_binning: Optional[QBinning] = None
    angle: Optional[AngleConfig] = None
    peak_fit: Optional[PeakFitConfig] = None
    resolution: Optional[ResolutionConfig] = None
    geometry: Optional[GeometryOverride] = None
    scale_factor: Optional[float] = None  # manual stitch scale override
    trim: Optional[tuple[int, int]] = None  # (pre_cut, post_cut) points trimmed from curve ends


# --- the two run kinds ------------------------------------------------------


class DirectBeamConfig(RunParameters):
    """A named composite direct beam (§3.3.3): the individual DB runs it is built from.

    Inherits `RunParameters` (rather than carrying only `db_runs`) because the direct beam has
    its own peak/background/low-res characterization window, distinct from any reflected run
    that references it, and should inherit from the same global `defaults`.
    """

    db_runs: list[int] = Field(min_length=1)
    # OPEN: not explicitly specified for DB runs; included for symmetry with reflected runs (§3.1.4)
    filter: Optional[RunFilter] = None


class ReflectedRunConfig(ReflectedRunParameters):
    """One reflected run keyed by sequence_number (§3.3.4).

    Inherits the optional `ReflectedRunParameters` fields as per-run overrides atop the global
    `defaults`, and carries the run's source data (§3.3.6.1) and optional filter (§3.3.6.2).
    """

    sequence_number: int = Field(ge=1)
    direct_beam: str  # named DB reference (§3.3.4)

    # Exactly one of these identifies the reflected run's source data (§3.3.6.1).
    run_number: Optional[int] = None  # OPEN: could instead be resolved from NeXus logs by sequence_number
    source_runs: list[int] | None = None  # multiple runs to sum, reflected only (hard-fail
    # on mismatched conditions is enforced at the summing operation, not at schema load)

    filter: Optional[RunFilter] = None  # §3.3.6.2

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
    def resolved_source_runs(self) -> list[int]:
        """The run(s) to load and (if more than one) sum for this reflected run."""
        return list(self.source_runs) if self.source_runs else [self.run_number]


# --------------------------------------------------------------------------- global containers ----


class OutputConfig(BaseModel):
    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    stem: str = "reduction_output"
    subtitle: Optional[str] = None
    eight_column: bool = False  # legacy 8-column output format (T, L, dT, dL)


class DiagnosticsConfig(BaseModel):
    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    plots: bool = True
    rq4: bool = False  # R*Q^4 view


class AssemblyConfig(BaseModel):
    """Sequence-assembly parameters (§2.4) — apply to the combined curve, not one run.

    OPEN: `autoscale` and `stitching_type` overlap conceptually (legacy `AutoScale` flag vs.
    `StitchingType.AUTOMATIC_AVERAGE`); reconcile with instrument scientists.
    """

    model_config = ConfigDict(extra="forbid", validate_assignment=True)
    critical_edge_normalization: bool = False
    q_norm: float = Field(0.015, gt=0)
    autoscale: bool = False
    stitching_type: Literal["manual", "automatic_average", "absolute_normalization"] = "manual"
    scale_factor_q_range: Optional[tuple[float, float]] = None
    normalize_first_angle: bool = False


# ------------------------------------------------------------------------------- top level ----


class ReductionConfig(BaseModel):
    """Global configuration (§3.3.1)."""

    model_config = ConfigDict(extra="forbid", validate_assignment=True)

    schema_version: Literal["1.0"] = "1.0"

    # global settings
    instrument: str = "BL4B"
    output: OutputConfig = Field(default_factory=OutputConfig)
    diagnostics: DiagnosticsConfig = Field(default_factory=DiagnosticsConfig)
    assembly: AssemblyConfig = Field(default_factory=AssemblyConfig)
    corrections: Corrections = Field(default_factory=Corrections)

    # per-run defaults inherited by every run (§3.3.1 global vs. per-run split); the superset
    # type since a reflected run's defaults may cover fields a direct beam has no use for
    defaults: ReflectedRunParameters = Field(default_factory=ReflectedRunParameters)

    # the two per-run key spaces (§3.3.2/.3/.4)
    direct_beams: dict[str, DirectBeamConfig]
    runs: list[ReflectedRunConfig]

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

    def effective(self, run: Union[ReflectedRunConfig, DirectBeamConfig]) -> RunParameters:
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
