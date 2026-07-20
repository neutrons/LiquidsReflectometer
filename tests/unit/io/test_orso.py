import numpy as np
import pytest
from orsopy import fileio
from orsopy.fileio.base import Column, Person, Value
from orsopy.fileio.data_source import DataSource, Experiment, Measurement, Sample
from orsopy.fileio.data_source import InstrumentSettings as OrsoInstrumentSettings
from orsopy.fileio.reduction import Reduction

from lr_reduction.io.orso import read_config
from lr_reduction.models.config import ReductionConfig


def _write_orso_file(path, configuration_yaml=None):
    instrument_settings = OrsoInstrumentSettings(incident_angle=Value(0.5, "deg"), wavelength=Value(4.0, "angstrom"))
    data_source = DataSource(
        owner=Person(name="tester", affiliation="ORNL"),
        experiment=Experiment(
            title="test", instrument="BL4B", start_date="2026-01-01T00:00:00", probe="neutron", facility="SNS"
        ),
        sample=Sample(name="sample"),
        measurement=Measurement(instrument_settings=instrument_settings, data_files=[]),
    )
    reduction = Reduction(software=fileio.Software(name="lr_reduction"))
    extra = {"configuration_yaml": configuration_yaml} if configuration_yaml is not None else {}
    header = fileio.orso.Orso(
        data_source=data_source, reduction=reduction, columns=[Column("Qz"), Column("R")], **extra
    )
    data = np.array([[0.01, 1.0], [0.02, 0.9]])
    fileio.save_orso([fileio.orso.OrsoDataset(info=header, data=data)], str(path))


def _example_config_dict():
    return ReductionConfig(
        sequence_id=12340,
        direct_beams={"db_8mm": {"db_runs": [12345]}},
        runs=[{"sequence_number": 1, "direct_beam": "db_8mm", "run_number": 111}],
    ).model_dump()


def test_read_config(tmp_path):
    path = tmp_path / "output.ort"
    _write_orso_file(path, configuration_yaml=_example_config_dict())

    config = read_config(str(path))

    assert config.sequence_id == 12340
    assert set(config.direct_beams) == {"db_8mm"}
    assert config.runs[0].run_number == 111


def test_read_config_raises_without_embedded_config(tmp_path):
    path = tmp_path / "no_config.ort"
    _write_orso_file(path)

    with pytest.raises(ValueError, match="no embedded"):
        read_config(str(path))


def test_read_partials(template_dir): ...


# TODO: read output file and check contents
def test_write_orso(template_dir): ...
