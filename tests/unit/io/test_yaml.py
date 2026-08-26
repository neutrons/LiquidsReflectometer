from pathlib import Path

from lr_reduction.io.yaml import read_config


def test_read_config(template_dir):
    config = read_config(str(Path(template_dir) / "config_example.yaml"))
    assert config.instrument == "BL4B"
    assert set(config.direct_beams) == {"db_8mm", "db_20mm"}
    assert [r.sequence_number for r in config.runs.values()] == [1, 2, 3]
    assert config.runs[1].direct_beam == "db_8mm"
    assert config.runs[1].peak.min == 140
    # a run may inherit global defaults instead of overriding them
    assert config.runs[3].peak is None
    assert config.effective(config.runs[3]).q_binning.q_max == 0.5


def test_read_config_preserves_zero_padded_run_number(tmp_path):
    # an unquoted, zero-padded run number must not be misread as octal (e.g. 012345 -> 5349)
    config_file = tmp_path / "config.yaml"
    config_file.write_text(
        "direct_beams:\n"
        "  db_8mm:\n"
        "    db_run_numbers: [012345]\n"
        "runs:\n"
        "  - sequence_number: 1\n"
        "    direct_beam: db_8mm\n"
        "    run_number: 012345\n"
    )
    config = read_config(str(config_file))
    assert config.direct_beams["db_8mm"].db_run_numbers == [12345]
    assert config.runs[1].run_number == 12345
