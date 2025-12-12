import pytest

from lr_reduction.template import get_default_template_file


def test_get_default_template_file_up_down_and_default(tmp_path):
    out = tmp_path

    # No template -> ValueError
    with pytest.raises(ValueError):
        get_default_template_file(str(out), 0)

    # Case: fallback default template
    default_file = out / "template.xml"
    default_file.touch()
    # tthd doesn't matter when only default exists
    assert get_default_template_file(str(out), 0) == str(default_file)

    # Case: up geometry
    up_file = out / "template_up.xml"
    up_file.touch()
    assert get_default_template_file(str(out), 1) == str(up_file)

    # Case: down geometry
    down_file = out / "template_down.xml"
    down_file.touch()
    assert get_default_template_file(str(out), -1) == str(down_file)
