import pytest
from pydantic import ValidationError

from lr_reduction.models.config import (
    BackgroundConfig,
    DirectBeamTrimConfig,
    EmissionTimeConfig,
    PeakFitConfig,
    PixelRange,
    QBinning,
    ReductionConfig,
    apply_overrides,
)


def _minimal_config(**overrides):
    data = dict(
        direct_beams={"db_8mm": {"direct_beam_run_numbers": [12345, 12346, 12347]}},
        runs={
            1: {"sequence_number": 1, "direct_beam": "db_8mm", "run_number": 111},
            2: {"sequence_number": 2, "direct_beam": "db_8mm", "run_number": 112},
        },
    )
    data.update(overrides)
    return ReductionConfig(**data)


def test_minimal_config_is_valid():
    config = _minimal_config()
    assert len(config.runs) == 2


def test_referential_integrity_rejects_unknown_direct_beam():
    with pytest.raises(ValidationError, match="unknown direct_beam"):
        _minimal_config(runs={1: {"sequence_number": 1, "direct_beam": "does_not_exist", "run_number": 111}})


def test_rejects_runs_key_sequence_number_mismatch():
    with pytest.raises(ValidationError, match="does not match"):
        _minimal_config(runs={1: {"sequence_number": 2, "direct_beam": "db_8mm", "run_number": 111}})


def test_reflected_run_requires_exactly_one_source():
    # neither run_number nor source_runs
    with pytest.raises(ValidationError, match="exactly one of"):
        _minimal_config(runs={1: {"sequence_number": 1, "direct_beam": "db_8mm"}})
    # both run_number and source_runs
    with pytest.raises(ValidationError, match="exactly one of"):
        _minimal_config(
            runs={1: {"sequence_number": 1, "direct_beam": "db_8mm", "run_number": 111, "source_runs": [111, 112]}}
        )


def test_reflected_run_source_runs_summing():
    config = _minimal_config(runs={1: {"sequence_number": 1, "direct_beam": "db_8mm", "source_runs": [111, 112, 113]}})
    assert config.runs[1].resolved_source_runs == [111, 112, 113]


def test_pixel_range_rejects_inverted_bounds():
    with pytest.raises(ValidationError, match="max .* < min"):
        PixelRange(min=150, max=100)


def test_q_binning_rejects_non_increasing_range():
    with pytest.raises(ValidationError, match="q_max .* must exceed q_min"):
        QBinning(q_min=0.5, q_max=0.1)


def test_emission_time_coefficients_accept_t0_t1():
    assert EmissionTimeConfig(coefficients=[10.0, -2.5]).coefficients == [10.0, -2.5]


@pytest.mark.parametrize("coefficients", [[], [1.0], [1.0, 2.0, 3.0]])
def test_emission_time_coefficients_must_be_exactly_two(coefficients):
    with pytest.raises(ValidationError, match=r"exactly \[t0, t1\]"):
        EmissionTimeConfig(coefficients=coefficients)


@pytest.mark.parametrize("bad", [float("nan"), float("inf")])
def test_emission_time_coefficients_must_be_finite(bad):
    with pytest.raises(ValidationError, match="finite"):
        EmissionTimeConfig(coefficients=[1.0, bad])


def test_for_run_extracts_single_run_and_referenced_direct_beam():
    config = _minimal_config(
        direct_beams={
            "db_8mm": {"direct_beam_run_numbers": [12345]},
            "db_20mm": {"direct_beam_run_numbers": [12350]},
        },
        runs={
            1: {"sequence_number": 1, "direct_beam": "db_8mm", "run_number": 111},
            2: {"sequence_number": 2, "direct_beam": "db_20mm", "run_number": 112},
        },
    )
    single = config.for_run(2)
    assert [r.sequence_number for r in single.runs.values()] == [2]
    assert set(single.direct_beams) == {"db_20mm"}
    # still passes the same validation as any other ReductionConfig
    assert isinstance(single, ReductionConfig)


def test_for_run_raises_for_unknown_sequence_number():
    config = _minimal_config()
    with pytest.raises(ValueError, match="No run with sequence_number"):
        config.for_run(999)


def test_effective_merges_defaults_with_run_overrides():
    config = _minimal_config(
        defaults={"peak": {"min": 140, "max": 150}, "scale_factor": 1.0},
        runs={
            1: {
                "sequence_number": 1,
                "direct_beam": "db_8mm",
                "run_number": 111,
                "peak": {"min": 100, "max": 110},
            }
        },
    )
    effective = config.effective(config.runs[1])
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
    config.runs[1].peak = {"min": 100, "max": 110}
    assert config.runs[1].peak == PixelRange(min=100, max=110)

    with pytest.raises(ValidationError, match="max .* < min"):
        config.runs[1].peak = {"min": 200, "max": 100}


def test_direct_attribute_assignment_does_not_recheck_other_models():
    # documented limitation: mutating a nested model only re-runs *that* model's own
    # validators, not a parent model's cross-field checks like referential integrity;
    # apply_overrides (or re-validating the whole tree) is required for that.
    config = _minimal_config()
    config.runs[1].direct_beam = "does_not_exist"
    assert config.runs[1].direct_beam == "does_not_exist"

    with pytest.raises(ValidationError, match="unknown direct_beam"):
        apply_overrides(config, {})


def test_background_requires_strips_when_applied():
    with pytest.raises(ValidationError, match="sequence_number=1: background subtraction is on but no strips"):
        _minimal_config(defaults={"background": {"apply": True}})
    _minimal_config(defaults={"background": {"apply": False}})


def test_background_partial_override_inherits_from_defaults():
    defaults = {
        "peak": {"min": 140, "max": 150},
        "background": {"model": "constant", "strips": [{"min": 120, "max": 130}, {"min": 160, "max": 170}]},
    }
    runs = {
        1: {"sequence_number": 1, "direct_beam": "db_8mm", "run_number": 111, "background": {"model": "linear"}},
        2: {
            "sequence_number": 2,
            "direct_beam": "db_8mm",
            "run_number": 112,
            "background": {"strips": [{"min": 120, "max": 130}]},
        },
    }
    config = _minimal_config(defaults=defaults, runs=runs)
    first, second = (config.effective(config.runs[n]).background for n in (1, 2))
    assert (first.model, len(first.strips)) == ("linear", 2)
    assert (second.model, len(second.strips)) == ("constant", 1)


def test_background_accepts_at_most_two_strips():
    strip = {"min": 1, "max": 2}
    with pytest.raises(ValidationError):
        BackgroundConfig(strips=[strip, strip, strip])


def test_background_strips_must_leave_rows_outside_peak():
    defaults = {"peak": {"min": 140, "max": 150}, "background": {"strips": [{"min": 128, "max": 137}]}}
    assert _minimal_config(defaults=defaults).defaults.background.model == "linear"

    # one row left outside the peak: enough for `constant`, too few for `linear`
    defaults["background"] = {"model": "constant", "strips": [{"min": 139, "max": 150}]}
    _minimal_config(defaults=defaults)
    defaults["background"]["model"] = "linear"
    with pytest.raises(ValidationError, match="sequence_number=1.*linear model needs at least 2"):
        _minimal_config(defaults=defaults)

    # a per-run override is checked against the run's effective peak
    defaults["background"] = {"strips": [{"min": 128, "max": 137}]}
    runs = {1: {"sequence_number": 1, "direct_beam": "db_8mm", "run_number": 111, "peak": {"min": 120, "max": 140}}}
    with pytest.raises(ValidationError, match="sequence_number=1"):
        _minimal_config(defaults=defaults, runs=runs)


def test_background_partial_override_survives_for_run_and_apply_overrides():
    defaults = {
        "peak": {"min": 140, "max": 150},
        "background": {"model": "constant", "strips": [{"min": 120, "max": 130}, {"min": 160, "max": 170}]},
    }
    runs = {
        1: {
            "sequence_number": 1,
            "direct_beam": "db_8mm",
            "run_number": 111,
            "background": {"strips": [{"min": 120, "max": 130}]},
        },
    }
    config = _minimal_config(defaults=defaults, runs=runs)
    for derived in (config.for_run(1), apply_overrides(config, {})):
        assert derived.effective(derived.runs[1]).background.model == "constant"


def test_direct_beam_tof_trim_defaults_match_algorithm_defaults():
    from lr_reduction.operations.algorithms.trim_direct_beam_tof import TrimDirectBeamTof

    algorithm = TrimDirectBeamTof()
    algorithm.initialize()
    trim = DirectBeamTrimConfig()
    assert list(trim.chopper_cut.coefficients) == list(algorithm.getProperty("ChopperCutCoefficients").value)
    assert trim.chopper_cut.offset_ms == algorithm.getProperty("ChopperCutOffset").value
    assert trim.intensity_tof_step == algorithm.getProperty("IntensityTOFStep").value


def test_direct_beam_tof_trim_from_yaml_shaped_dict():
    config = _minimal_config(
        direct_beams={
            "db_8mm": {
                "direct_beam_run_numbers": [12345],
                "flip_atten": True,
                "tof_trim": {"apply": False, "chopper_cut": {"coefficients": [2.0, 16000.0], "offset_ms": 0.5}},
            }
        }
    )
    db = config.direct_beams["db_8mm"]
    assert db.flip_atten is True
    assert db.tof_trim.apply is False
    assert db.tof_trim.chopper_cut.coefficients == (2.0, 16000.0)
    assert db.tof_trim.chopper_cut.offset_ms == 0.5
    # unspecified fields keep their defaults
    assert db.tof_trim.intensity_tof_step == 50.0


@pytest.mark.parametrize(
    "tof_trim",
    [
        {"dead_time_threshold": 0.9},
        {"dead_time_threshold_lowest_wavelength": 0.9},
        {"intensity_threshold": -1.0},
        {"intensity_tof_step": 0.0},
        {"chopper_cut": {"coefficients": [1.0]}},
    ],
)
def test_direct_beam_tof_trim_rejects_invalid_values(tof_trim):
    with pytest.raises(ValidationError):
        DirectBeamTrimConfig(**tof_trim)


def test_peak_fit_defaults():
    peak_fit = PeakFitConfig()
    assert (peak_fit.type, peak_fit.background, peak_fit.pad) == ("supergauss", "linear", 1)
    with pytest.raises(ValidationError):
        PeakFitConfig(background="quadratic")


def test_direct_beam_inherits_peak_fit_from_defaults():
    config = _minimal_config(defaults={"peak_fit": {"type": "gauss", "pad": 3}})
    assert config.effective(config.direct_beams["db_8mm"]).peak_fit == PeakFitConfig(type="gauss", pad=3)


def test_edits_to_default_created_blocks_survive_for_run_and_apply_overrides():
    config = _minimal_config()
    config.assembly.q_norm = 0.02
    config.corrections.dead_time.apply = not config.corrections.dead_time.apply
    config.direct_beams["db_8mm"].tof_trim.chopper_cut.offset_ms = 2.5
    config.defaults.peak = PixelRange(min=140, max=150)
    for derived in (config.for_run(1), apply_overrides(config, {})):
        assert derived.assembly.q_norm == 0.02
        assert derived.corrections.dead_time.apply == config.corrections.dead_time.apply
        assert derived.direct_beams["db_8mm"].tof_trim.chopper_cut.offset_ms == 2.5
        assert derived.defaults.peak == PixelRange(min=140, max=150)


def test_background_min_rows_match_algorithm():
    from lr_reduction.operations.algorithms.subtract_strip_background import MIN_BACKGROUND_ROWS

    # the operation passes `model.capitalize()` as BackgroundModel
    assert {model.capitalize(): rows for model, rows in BackgroundConfig.MIN_ROWS.items()} == MIN_BACKGROUND_ROWS
    assert set(BackgroundConfig.MIN_ROWS) == set(BackgroundConfig.model_fields["model"].annotation.__args__)


def test_invalid_merge_with_defaults_names_the_run():
    runs = {3: {"sequence_number": 3, "direct_beam": "db_8mm", "run_number": 113, "q_binning": {"q_min": 0.3}}}
    with pytest.raises(ValidationError, match=r"(?s)sequence_number=3, merged with `defaults`:.*q_max \(0.2\)"):
        _minimal_config(defaults={"q_binning": {"q_min": 0.01, "q_max": 0.2}}, runs=runs)


def _partial_override_config():
    """Partial per-run blocks over `defaults`, in a background block and in two others."""
    return _minimal_config(
        defaults={
            "peak": {"min": 140, "max": 150},
            "peak_fit": {"type": "gauss"},
            "q_binning": {"method": "constantTOF"},
            "background": {"model": "constant", "strips": [{"min": 120, "max": 130}]},
        },
        direct_beams={"db_8mm": {"direct_beam_run_numbers": [12345], "peak_fit": {"pad": 3}}},
        runs={
            1: {
                "sequence_number": 1,
                "direct_beam": "db_8mm",
                "run_number": 111,
                "q_binning": {"q_step": 0.01},
                "background": {"model": "linear"},
            }
        },
    )


def _resolved(config):
    run, beam = config.runs[1], config.direct_beams["db_8mm"]
    return config.effective(run), config.effective(beam)


def test_partial_blocks_inherit_unwritten_keys():
    run, beam = _resolved(_partial_override_config())
    assert (beam.peak_fit.type, beam.peak_fit.pad) == ("gauss", 3)
    assert (run.q_binning.method, run.q_binning.q_step) == ("constantTOF", 0.01)
    assert (run.background.model, run.background.strips) == ("linear", [PixelRange(min=120, max=130)])


@pytest.mark.parametrize("mode", ["python", "json"])
def test_model_dump_round_trip_keeps_effective_parameters(mode):
    config = _partial_override_config()
    rebuilt = ReductionConfig.model_validate(config.model_dump(mode=mode))
    assert _resolved(rebuilt) == _resolved(config)


def test_in_place_edit_of_a_run_block_is_kept():
    config = _partial_override_config()
    config.runs[1].background.strips.append(PixelRange(min=160, max=170))
    expected = [PixelRange(min=160, max=170)]
    for derived in (config, config.for_run(1), apply_overrides(config, {})):
        assert derived.effective(derived.runs[1]).background.strips == expected


def test_serialization_schema_keeps_run_fields():
    definitions = ReductionConfig.model_json_schema(mode="serialization")["$defs"]
    for model in ("ReflectedRunConfig", "DirectBeamConfig", "ReflectedRunParameters"):
        assert {"peak", "background"} <= set(definitions[model]["properties"])
