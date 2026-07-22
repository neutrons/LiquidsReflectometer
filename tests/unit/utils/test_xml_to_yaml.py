from pathlib import Path

from lr_reduction.io.yaml import read_config as read_yaml_config
from lr_reduction.legacy.reduction_template_reader import ReductionParameters
from lr_reduction.utils.xml_to_yaml import _run_config_from_legacy, convert, read_config


def test_read_config_maps_runs_and_dedups_direct_beams(template_dir):
    config = read_config(str(Path(template_dir) / "template_201282.xml"))

    assert [r.sequence_number for r in config.runs] == list(range(1, len(config.runs) + 1))
    # composite direct beams are named after the legacy single DB run and shared across runs
    assert len(config.direct_beams) < len(config.runs)
    for db in config.direct_beams.values():
        assert len(db.db_runs) == 1


def test_read_config_maps_stitching_and_q_binning(template_dir):
    config = read_config(str(Path(template_dir) / "template_with_const_q_true.xml"))
    assert config.runs[0].q_binning.method == "constantQ"

    config = read_config(str(Path(template_dir) / "template_stitching_automatic_average.xml"))
    assert config.assembly.stitching_type == "automatic_average"


def test_read_config_maps_instrument_geometry_override(template_dir):
    config = read_config(str(Path(template_dir) / "template_with_instrument_settings.xml"))
    assert config.runs[0].geometry is not None
    assert config.runs[0].geometry.l2 == config.runs[0].geometry.l2  # sample_detector_distance, sanity check


def test_run_config_from_legacy_maps_second_background_region():
    params = ReductionParameters()
    params.two_backgrounds = True
    params.background_roi = [130, 140, 150, 160]
    params.data_files = [111]
    params.norm_file = 222

    run, db_name, direct_beam = _run_config_from_legacy(params, sequence_number=1)

    assert run.background.region2 is not None
    assert (run.background.region2.min, run.background.region2.max) == (150, 160)
    assert db_name == "db_222"
    assert direct_beam.db_runs == [222]


def test_run_config_from_legacy_source_runs_summing():
    params = ReductionParameters()
    params.data_files = [111, 112, 113]
    params.norm_file = 222

    run, _, _ = _run_config_from_legacy(params, sequence_number=1)

    assert run.run_number is None
    assert run.source_runs == [111, 112, 113]
    assert run.resolved_source_runs == [111, 112, 113]


def test_convert_round_trips_through_native_yaml_loader(tmp_path, template_dir):
    xml_path = Path(template_dir) / "template_201282.xml"
    yaml_path = tmp_path / "converted.yaml"

    convert(str(xml_path), str(yaml_path))

    direct = read_config(str(xml_path))
    via_yaml = read_yaml_config(str(yaml_path))

    assert set(via_yaml.direct_beams) == set(direct.direct_beams)
    assert [r.sequence_number for r in via_yaml.runs] == [r.sequence_number for r in direct.runs]
    assert [r.run_number for r in via_yaml.runs] == [r.run_number for r in direct.runs]
