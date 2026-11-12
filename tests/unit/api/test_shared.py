from lr_reduction.api._shared import get_direct_beam_name, get_direct_beams
from lr_reduction.models.config import DirectBeamConfig, ReductionConfig, ReflectedRunConfig


def _config() -> ReductionConfig:
    return ReductionConfig(
        direct_beams={
            "db_8mm": DirectBeamConfig(direct_beam_run_numbers=[11111, 11112]),
            "db_20mm": DirectBeamConfig(direct_beam_run_numbers=[22222]),
        },
        runs={
            1: ReflectedRunConfig(sequence_number=1, direct_beam="db_8mm", run_number=12345),
            2: ReflectedRunConfig(sequence_number=2, direct_beam="db_20mm", run_number=12346),
        },
    )


def test_get_direct_beam_name_follows_the_runs_reference():
    config = _config()
    assert get_direct_beam_name(1, config) == "db_8mm"
    assert get_direct_beam_name(2, config) == "db_20mm"


def test_get_direct_beams_loads_the_named_composites_runs_in_order():
    class _RecordingLoader:
        def load(self, run_number):
            return run_number

    assert get_direct_beams(_RecordingLoader(), _config(), "db_8mm") == [11111, 11112]
