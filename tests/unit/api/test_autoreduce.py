from lr_reduction.api.autoreduce import AutoreduceSingleRun, FromDiskSequence, reduce_auto
from lr_reduction.io import RunData
from lr_reduction.models.config import DirectBeamConfig, ReductionConfig, ReflectedRunConfig
from lr_reduction.models.results import CombinedReductionResult, ReductionResult


def _config() -> ReductionConfig:
    return ReductionConfig(
        direct_beams={"db": DirectBeamConfig(db_run_numbers=[11111])},
        runs={1: ReflectedRunConfig(sequence_number=1, direct_beam="db", run_number=178200)},
    )


def test_reduce_auto_executes_end_to_end(tmp_path, monkeypatch):
    """Runs the single-run leaf then FromDiskSequence, proving the incremental on-disk path (§6.4.4) works today."""
    monkeypatch.setattr("lr_reduction.io.config_loader.ConfigLoader.load", lambda _self, _path: _config())
    nexus_file = tmp_path / "REF_L_178200.nxs.h5"
    nexus_file.touch()

    result = reduce_auto(str(nexus_file), str(tmp_path))

    assert isinstance(result, CombinedReductionResult)


def test_reduce_auto_loads_a_relocated_nexus_file(tmp_path, monkeypatch):
    """The nexus path is loaded directly (no run-number rediscovery), so a file that doesn't
    follow the REF_L_<run_number> naming convention or live at the conventional location
    still reduces successfully."""
    monkeypatch.setattr("lr_reduction.io.config_loader.ConfigLoader.load", lambda _self, _path: _config())
    nexus_file = tmp_path / "reprocessed_data.nxs.h5"
    nexus_file.touch()

    result = reduce_auto(str(nexus_file), str(tmp_path))

    assert isinstance(result, CombinedReductionResult)


def test_autoreduce_single_run_loads_data_from_the_given_path(tmp_path, monkeypatch):
    monkeypatch.setattr("lr_reduction.io.config_loader.ConfigLoader.load", lambda _self, _path: _config())
    captured_paths = []

    def _capture_path(_self, nexus_file_path):
        captured_paths.append(nexus_file_path)
        return RunData()

    monkeypatch.setattr("lr_reduction.io.run_loader.RunLoader.load_from_path", _capture_path)
    nexus_file = tmp_path / "reprocessed_data.nxs.h5"
    nexus_file.touch()

    result = AutoreduceSingleRun(nexus_file, tmp_path).execute()

    assert isinstance(result, ReductionResult)
    assert captured_paths == [nexus_file]


def test_from_disk_sequence_reuses_read_partials(tmp_path, monkeypatch):
    monkeypatch.setattr("lr_reduction.io.config_loader.ConfigLoader.load", lambda _self, _path: _config())

    result = FromDiskSequence(tmp_path).execute()

    assert isinstance(result, CombinedReductionResult)


def test_from_disk_sequence_reads_partials_for_the_output_directory(tmp_path, monkeypatch):
    monkeypatch.setattr("lr_reduction.io.config_loader.ConfigLoader.load", lambda _self, _path: _config())
    captured_sequence_ids = []

    def _capture_sequence_id(_partial_dir, sequence_id):
        captured_sequence_ids.append(sequence_id)
        return []

    monkeypatch.setattr("lr_reduction.api.autoreduce.read_partials", _capture_sequence_id)

    FromDiskSequence(tmp_path).execute()

    assert captured_sequence_ids == [0]
