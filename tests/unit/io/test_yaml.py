from pathlib import Path

from lr_reduction.io.yaml import read_config


def test_read_config(template_dir):
    config = read_config(str(Path(template_dir) / "config_example.yaml"))
    assert config.instrument == "BL4B"
    assert set(config.direct_beams) == {"db_8mm", "db_20mm"}
    assert [r.sequence_number for r in config.runs] == [1, 2, 3]
    assert config.runs[0].direct_beam == "db_8mm"
    assert config.runs[0].peak.min == 140
    # a run may inherit global defaults instead of overriding them
    assert config.runs[2].peak is None
    assert config.effective(config.runs[2]).q_binning.q_max == 0.5
