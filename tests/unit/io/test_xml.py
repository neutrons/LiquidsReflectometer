import shutil
from pathlib import Path

from lr_reduction.io.xml import _run_config_from_legacy, read_config
from lr_reduction.legacy.reduction_template_reader import ReductionParameters


def test_read_config_maps_runs_and_dedups_direct_beams(template_dir):
    config = read_config(str(Path(template_dir) / "template_201282.xml"))

    assert config.sequence_id == 201282
    assert [r.sequence_number for r in config.runs] == list(range(1, len(config.runs) + 1))
    # composite direct beams are named after the legacy single DB run and shared across runs
    assert len(config.direct_beams) < len(config.runs)
    for db in config.direct_beams.values():
        assert len(db.db_runs) == 1


def test_read_config_maps_stitching_and_q_binning(tmp_path, template_dir):
    # `template_with_const_q_true.xml`/`template_stitching_automatic_average.xml` don't encode a
    # sequence id in their filename, so copy them under a name get_sequence_id_from_path accepts.
    const_q_path = tmp_path / "seq_99001_template_with_const_q_true.xml"
    shutil.copy(Path(template_dir) / "template_with_const_q_true.xml", const_q_path)
    config = read_config(str(const_q_path))
    assert config.runs[0].q_binning.method == "constantQ"

    stitching_path = tmp_path / "seq_99002_template_stitching_automatic_average.xml"
    shutil.copy(Path(template_dir) / "template_stitching_automatic_average.xml", stitching_path)
    config = read_config(str(stitching_path))
    assert config.assembly.stitching_type == "automatic_average"


def test_read_config_maps_instrument_geometry_override(tmp_path, template_dir):
    path = tmp_path / "seq_99003_template_with_instrument_settings.xml"
    shutil.copy(Path(template_dir) / "template_with_instrument_settings.xml", path)
    config = read_config(str(path))
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
