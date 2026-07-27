import pytest
from pydantic import ValidationError

from lr_reduction.models.config import PixelRange, QBinning, ReductionConfig, apply_overrides


def _minimal_config(**overrides):
    data = dict(
        direct_beams={"db_8mm": {"db_runs": [12345, 12346, 12347]}},
        runs=[
            {"sequence_number": 1, "direct_beam": "db_8mm", "run_number": 111},
            {"sequence_number": 2, "direct_beam": "db_8mm", "run_number": 112},
        ],
    )
    data.update(overrides)
    return ReductionConfig(**data)


def test_minimal_config_is_valid():
    config = _minimal_config()
    assert len(config.runs) == 2


def test_referential_integrity_rejects_unknown_direct_beam():
    with pytest.raises(ValidationError, match="unknown direct_beam"):
        _minimal_config(runs=[{"sequence_number": 1, "direct_beam": "does_not_exist", "run_number": 111}])


def test_rejects_duplicate_sequence_number():
    with pytest.raises(ValidationError, match="duplicate sequence_number"):
        _minimal_config(
            runs=[
                {"sequence_number": 1, "direct_beam": "db_8mm", "run_number": 111},
                {"sequence_number": 1, "direct_beam": "db_8mm", "run_number": 112},
            ]
        )


def test_reflected_run_requires_exactly_one_source():
    # neither run_number nor source_runs
    with pytest.raises(ValidationError, match="exactly one of"):
        _minimal_config(runs=[{"sequence_number": 1, "direct_beam": "db_8mm"}])
    # both run_number and source_runs
    with pytest.raises(ValidationError, match="exactly one of"):
        _minimal_config(
            runs=[{"sequence_number": 1, "direct_beam": "db_8mm", "run_number": 111, "source_runs": [111, 112]}]
        )


def test_reflected_run_source_runs_summing():
    config = _minimal_config(runs=[{"sequence_number": 1, "direct_beam": "db_8mm", "source_runs": [111, 112, 113]}])
    assert config.runs[0].resolved_source_runs == [111, 112, 113]


def test_pixel_range_rejects_inverted_bounds():
    with pytest.raises(ValidationError, match="max .* < min"):
        PixelRange(min=150, max=100)


def test_q_binning_rejects_non_increasing_range():
    with pytest.raises(ValidationError, match="q_max .* must exceed q_min"):
        QBinning(q_min=0.5, q_max=0.1)


def test_for_run_extracts_single_run_and_referenced_direct_beam():
    config = _minimal_config(
        direct_beams={
            "db_8mm": {"db_runs": [12345]},
            "db_20mm": {"db_runs": [12350]},
        },
        runs=[
            {"sequence_number": 1, "direct_beam": "db_8mm", "run_number": 111},
            {"sequence_number": 2, "direct_beam": "db_20mm", "run_number": 112},
        ],
    )
    single = config.for_run(2)
    assert [r.sequence_number for r in single.runs] == [2]
    assert set(single.direct_beams) == {"db_20mm"}
    # still passes the same validation as any other ReductionConfig
    assert isinstance(single, ReductionConfig)


def test_for_run_raises_for_unknown_sequence_number():
    config = _minimal_config()
    with pytest.raises(ValueError, match="no run with sequence_number"):
        config.for_run(999)


def test_effective_merges_defaults_with_run_overrides():
    config = _minimal_config(
        defaults={"peak": {"min": 140, "max": 150}, "scale_factor": 1.0},
        runs=[
            {
                "sequence_number": 1,
                "direct_beam": "db_8mm",
                "run_number": 111,
                "peak": {"min": 100, "max": 110},
            }
        ],
    )
    effective = config.effective(config.runs[0])
    # run-level override wins over the global default
    assert effective.peak.min == 100
    assert effective.peak.max == 110
    # fields not overridden by the run are inherited from defaults
    assert effective.scale_factor == 1.0


def test_apply_overrides_takes_precedence_over_file():
    config = _minimal_config(assembly={"q_norm": 0.015})
    overridden = apply_overrides(config, {"assembly": {"q_norm": 0.05}})
    assert overridden.assembly.q_norm == 0.05
    # original config is untouched
    assert config.assembly.q_norm == 0.015


def test_apply_overrides_does_not_mutate_original():
    config = _minimal_config()
    apply_overrides(config, {"instrument": "SOMETHING_ELSE"})
    assert config.instrument == "BL4B"


def test_direct_attribute_assignment_is_validated():
    config = _minimal_config()

    # a valid direct override is simpler than building an overrides dict
    config.assembly.q_norm = 0.02
    assert config.assembly.q_norm == 0.02

    # out-of-range values are still rejected immediately, per that field's own constraints
    with pytest.raises(ValidationError):
        config.assembly.q_norm = -1


def test_direct_attribute_assignment_coerces_and_validates_nested_models():
    config = _minimal_config()

    # assigning a plain dict to a nested model field constructs and validates it, same as at load
    config.runs[0].peak = {"min": 100, "max": 110}
    assert config.runs[0].peak == PixelRange(min=100, max=110)

    with pytest.raises(ValidationError, match="max .* < min"):
        config.runs[0].peak = {"min": 200, "max": 100}


def test_direct_attribute_assignment_does_not_recheck_other_models():
    # documented limitation: mutating a nested model only re-runs *that* model's own
    # validators, not a parent model's cross-field checks like referential integrity;
    # apply_overrides (or re-validating the whole tree) is required for that.
    config = _minimal_config()
    config.runs[0].direct_beam = "does_not_exist"
    assert config.runs[0].direct_beam == "does_not_exist"

    with pytest.raises(ValidationError, match="unknown direct_beam"):
        apply_overrides(config, {})
