from pathlib import Path

from lr_reduction.io.xml import read_config as read_xml_config
from lr_reduction.io.xml_to_yaml import convert
from lr_reduction.io.yaml import read_config as read_yaml_config


def test_convert_round_trips_through_native_yaml_loader(tmp_path, template_dir):
    xml_path = Path(template_dir) / "template_201282.xml"
    yaml_path = tmp_path / "converted.yaml"

    convert(str(xml_path), str(yaml_path))

    direct = read_xml_config(str(xml_path))
    via_yaml = read_yaml_config(str(yaml_path))

    assert via_yaml.sequence_id == direct.sequence_id
    assert set(via_yaml.direct_beams) == set(direct.direct_beams)
    assert [r.sequence_number for r in via_yaml.runs] == [r.sequence_number for r in direct.runs]
    assert [r.run_number for r in via_yaml.runs] == [r.run_number for r in direct.runs]
