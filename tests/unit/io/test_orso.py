from datetime import datetime

import numpy as np
from orsopy.fileio import Software

from lr_reduction.io.orso import _run_label, write_orso
from lr_reduction.models import CombinedReductionResult, ReductionConfig, ReductionResult, ReflectivityCurve
from lr_reduction.models.config import DirectBeamConfig, ReflectedRunConfig


def test_read_orso(template_dir): ...


def test_read_partials(template_dir): ...


# TODO: read output file and check contents
def test_write_orso(tmp_path):
    """Test writing an ORSO file from a ReductionResult."""
    config = ReductionConfig()
    reduction_result = ReductionResult(
        curve=ReflectivityCurve(
            q=np.array([0.1, 0.2, 0.3]),
            r=np.array([0.1, 0.2, 0.3]),
            dr=np.array([0.01, 0.02, 0.03]),
            dq=np.array([0.01, 0.02, 0.03]),
        ),
        run_numbers=[1, 2, 3],
        sequence_id="1",
        sequence_number=1,
        reduction_config=config,
        lr_reduction_info=Software(name="lr_reduction", version="1.0.0"),
        mantid_info=Software(name="mantid", version="1.0.0"),
        reduction_timestamp=datetime.fromisoformat("2024-01-01T12:00:00"),
    )
    result = write_orso(results=reduction_result, output_dir=tmp_path, title="Test ORSO Output")
    assert result.exists()


def _software(name: str) -> Software:
    return Software(name=name, version="1.0.0")


def _partial(run_numbers: list[int], config: ReductionConfig | None = None) -> ReductionResult:
    return ReductionResult(
        curve=ReflectivityCurve.empty(),
        run_numbers=run_numbers,
        sequence_id=1,
        sequence_number=1,
        reduction_config=config or ReductionConfig(),
        lr_reduction_info=_software("lr_reduction"),
        mantid_info=_software("mantid"),
        reduction_timestamp=datetime.fromisoformat("2024-01-01T12:00:00"),
    )


def _combined(config: ReductionConfig, partials: list[ReductionResult]) -> CombinedReductionResult:
    return CombinedReductionResult(
        curve=ReflectivityCurve.empty(),
        reduction_config=config,
        lr_reduction_info=_software("lr_reduction"),
        mantid_info=_software("mantid"),
        reduction_timestamp=datetime.fromisoformat("2024-01-01T12:00:00"),
        partials=partials,
    )


def _sequence_config() -> ReductionConfig:
    return ReductionConfig(
        direct_beams={"db": DirectBeamConfig(run_numbers=[11111])},
        runs={
            1: ReflectedRunConfig(sequence_number=1, direct_beam="db", run_number=12360),
            2: ReflectedRunConfig(sequence_number=2, direct_beam="db", source_runs=[12361, 12362]),
        },
    )


def test_run_label_names_a_single_run():
    assert _run_label(_partial([12360])) == "run 12360"


def test_run_label_joins_summed_source_runs():
    assert _run_label(_partial([12361, 12362])) == "run 12361+12362"


def test_run_label_names_a_combined_result_by_its_partials():
    combined = _combined(_sequence_config(), [_partial([12360]), _partial([12361, 12362])])

    assert _run_label(combined) == "combined runs 12360, 12361+12362"


def test_run_label_falls_back_to_the_configured_runs_without_partials():
    assert _run_label(_combined(_sequence_config(), [])) == "combined runs 12360, 12361+12362"


def test_run_label_without_any_runs():
    assert _run_label(_combined(ReductionConfig(), [])) == "combined runs <unknown>"
