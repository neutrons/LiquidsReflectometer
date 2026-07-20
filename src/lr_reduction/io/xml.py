"""Legacy XML I/O module.

Reads legacy RefRed XML templates and converts them to a `ReductionConfig`, kept here for
backward compatibility during the coexistence period (this branch is `@deprecated`, see
`io/config_loader.py`). Per §3.3.10 the new workflow itself shall not read XML — the sanctioned
path from a legacy template into the new workflow is a standalone conversion utility (§3.3.10.2:
each legacy single direct-beam run is wrapped as a one-run composite direct beam), which reuses
the field-mapping helpers below so the ~40-field mapping logic lives in exactly one place.
"""

from lr_reduction.legacy.reduction_template_reader import ReductionParameters, from_xml
from lr_reduction.models.config import (
    AcceptanceWindow,
    AngleConfig,
    AssemblyConfig,
    BackgroundConfig,
    Corrections,
    DeadTimeConfig,
    DirectBeamConfig,
    EmissionTimeConfig,
    GeometryOverride,
    PixelRange,
    QBinning,
    ReductionConfig,
    ReflectedRunConfig,
)
from lr_reduction.utils import get_logger, get_sequence_id_from_path

logger = get_logger(__name__)

_STITCHING_TYPE_MAP = {
    "None": "manual",
    "AutomaticAverage": "automatic_average",
    "AbsoluteNormalization": "absolute_normalization",
}


def read_config(filepath: str) -> ReductionConfig:
    """Read a legacy XML template and convert it to a ReductionConfig."""
    logger.info(f"Reading legacy XML configuration from {filepath}")
    with open(filepath) as f:
        data_sets = from_xml(f.read())

    direct_beams: dict[str, DirectBeamConfig] = {}
    runs = []
    for sequence_number, params in enumerate(data_sets, start=1):
        run, db_name, direct_beam = _run_config_from_legacy(params, sequence_number)
        runs.append(run)
        direct_beams.setdefault(db_name, direct_beam)

    assembly = _assembly_config_from_legacy(data_sets[0]) if data_sets else AssemblyConfig()

    return ReductionConfig(
        sequence_id=int(get_sequence_id_from_path(filepath)),
        direct_beams=direct_beams,
        runs=runs,
        assembly=assembly,
    )


def _run_config_from_legacy(
    params: ReductionParameters, sequence_number: int
) -> tuple[ReflectedRunConfig, str, DirectBeamConfig]:
    """Map one legacy `ReductionParameters` instance to (reflected run, DB name, DB config).

    §3.3.10.2: the single legacy direct-beam run (`norm_file`) is wrapped as a one-run composite
    direct beam. §3.3.10.1: any new-workflow field with no legacy source simply keeps its schema
    default, since the models below are built via keyword construction rather than a raw dict.
    """
    db_name = f"db_{params.norm_file}"
    direct_beam = DirectBeamConfig(
        db_runs=[params.norm_file],
        peak=PixelRange(min=params.norm_peak_range[0], max=params.norm_peak_range[1]),
        background=BackgroundConfig(
            apply=params.subtract_norm_background,
            region=PixelRange(min=params.norm_background_roi[0], max=params.norm_background_roi[1]),
        ),
        low_res=(
            PixelRange(min=params.norm_x_range[0], max=params.norm_x_range[1]) if params.norm_x_range_flag else None
        ),
    )

    background = BackgroundConfig(
        apply=params.subtract_background,
        region=PixelRange(min=params.background_roi[0], max=params.background_roi[1]),
    )
    if params.two_backgrounds and params.background_roi[2:] != [0, 0]:
        background.region2 = PixelRange(min=params.background_roi[2], max=params.background_roi[3])

    geometry = None
    if params.apply_instrument_settings:
        # OPEN: legacy `source_detector_distance` is moderator-to-detector, not source-to-sample;
        # l1 is derived here rather than copied directly. Confirm with instrument scientists.
        geometry = GeometryOverride(
            l1=params.source_detector_distance - params.sample_detector_distance,
            l2=params.sample_detector_distance,
        )

    run = ReflectedRunConfig(
        sequence_number=sequence_number,
        direct_beam=db_name,
        run_number=params.data_files[0] if len(params.data_files) == 1 else None,
        source_runs=list(params.data_files) if len(params.data_files) > 1 else None,
        peak=PixelRange(min=params.data_peak_range[0], max=params.data_peak_range[1]),
        background=background,
        low_res=(
            PixelRange(min=params.data_x_range[0], max=params.data_x_range[1]) if params.data_x_range_flag else None
        ),
        acceptance=(AcceptanceWindow(tof=tuple(params.tof_range)) if params.select_tof_range else None),
        q_binning=QBinning(
            q_min=params.q_min,
            q_step=params.q_step,
            method="constantQ" if params.const_q else "constantTOF",
        ),
        angle=AngleConfig(theta_offset=params.angle_offset, theta_offset_error=params.angle_offset_error),
        corrections=Corrections(
            dead_time=DeadTimeConfig(
                apply=params.dead_time,
                value=params.dead_time_value,
                tof_step=params.dead_time_tof_step,
                paralyzable=params.paralyzable,
                threshold=(params.dead_time_threshold if params.use_dead_time_threshold else None),
            ),
            emission_time=EmissionTimeConfig(apply=params.use_emission_time),
            # OPEN: legacy `GravityDirection` also carries UP/DOWN; only the enable flag survives here.
            gravity=params.gravity_direction is not None,
        ),
        geometry=geometry,
        scale_factor=params.stitching_reflectivity_scale_factor,
        trim=(params.pre_cut, params.post_cut),
    )
    return run, db_name, direct_beam


def _assembly_config_from_legacy(params: ReductionParameters) -> AssemblyConfig:
    """Map the sequence-wide stitching configuration.

    Legacy repeats the same stitching configuration on every per-run block (a parallel-array
    artifact); it is genuinely sequence-wide, so only the first block's value is used.
    """
    stitching = params.stitching_configuration
    return AssemblyConfig(
        stitching_type=_STITCHING_TYPE_MAP[stitching.type.value],
        scale_factor_q_range=(stitching.scale_factor_qmin, stitching.scale_factor_qmax),
        normalize_first_angle=stitching.normalize_first_angle,
    )
