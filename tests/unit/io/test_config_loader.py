from pathlib import Path

import pytest
from pydantic import BaseModel, ValidationError

from lr_reduction.io.config_loader import (
    ConfigFileTypeError,
    ConfigLoader,
    ConfigNotFoundError,
    ConfigParseError,
    ConfigValidationError,
)
from lr_reduction.models.config import ReductionConfig


class _Dummy(BaseModel):
    """A model unrelated to ReductionConfig, used only to produce a generic ValidationError."""

    x: int


def _generic_validation_error() -> ValidationError:
    try:
        _Dummy(x="not-an-int")
    except ValidationError as exc:
        return exc
    raise AssertionError("expected ValidationError")  # pragma: no cover


def test_load_dispatches_valid_yaml(template_dir):
    config = ConfigLoader().load(str(Path(template_dir) / "config_example.yaml"))
    assert isinstance(config, ReductionConfig)
    assert config.instrument == "BL4B"


def test_load_raises_not_found_for_missing_file(tmp_path):
    missing = tmp_path / "does_not_exist.yaml"
    with pytest.raises(ConfigNotFoundError, match="not found") as exc_info:
        ConfigLoader().load(str(missing))
    assert exc_info.value.filepath == str(missing)


def test_load_raises_file_type_error_for_unsupported_extension(tmp_path):
    unsupported = tmp_path / "config.json"
    unsupported.touch()
    with pytest.raises(ConfigFileTypeError, match="Unsupported"):
        ConfigLoader().load(str(unsupported))


def test_load_raises_parse_error_for_invalid_yaml_syntax(tmp_path):
    bad_yaml = tmp_path / "broken.yaml"
    bad_yaml.write_text("instrument: [unclosed\n")
    with pytest.raises(ConfigParseError) as exc_info:
        ConfigLoader().load(str(bad_yaml))
    assert exc_info.value.filepath == str(bad_yaml)


def test_load_raises_parse_error_when_yaml_is_not_a_mapping(tmp_path):
    not_a_mapping = tmp_path / "list.yaml"
    not_a_mapping.write_text("- 1\n- 2\n")
    with pytest.raises(ConfigParseError, match="mapping"):
        ConfigLoader().load(str(not_a_mapping))


def test_load_raises_validation_error_on_pydantic_failure(tmp_path, monkeypatch):
    # The specific Pydantic model doesn't matter here; only that ConfigLoader correctly
    # translates any ValidationError raised by the underlying loader into a ConfigValidationError.
    config_file = tmp_path / "config.yaml"
    config_file.write_text("instrument: BL4B\n")
    validation_error = _generic_validation_error()
    monkeypatch.setattr(
        "lr_reduction.io.config_loader.read_yaml_config",
        lambda _path: (_ for _ in ()).throw(validation_error),
    )

    with pytest.raises(ConfigValidationError, match="x") as exc_info:
        ConfigLoader().load(str(config_file))

    assert exc_info.value.filepath == str(config_file)
    assert exc_info.value.errors == validation_error.errors()


def test_load_raises_not_found_when_underlying_loader_raises_file_not_found(tmp_path, monkeypatch):
    config_file = tmp_path / "config.yaml"
    config_file.write_text("instrument: BL4B\n")
    monkeypatch.setattr(
        "lr_reduction.io.config_loader.read_yaml_config",
        lambda path: (_ for _ in ()).throw(FileNotFoundError(path)),
    )

    with pytest.raises(ConfigNotFoundError):
        ConfigLoader().load(str(config_file))
