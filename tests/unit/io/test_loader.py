from lr_reduction.io.run_loader import RunLoader
from lr_reduction.models.run_data import RunData


def test_load_returns_run_data():
    assert isinstance(RunLoader().load(12345), RunData)


def test_load_from_path_returns_run_data(tmp_path):
    nexus_file = tmp_path / "REF_L_12345.nxs.h5"
    nexus_file.touch()

    assert isinstance(RunLoader().load_from_path(nexus_file), RunData)
