from lr_reduction.api.autoreduce import FromDiskSequence, reduce_autoreduce
from lr_reduction.models.config import DirectBeamConfig, ReductionConfig, ReflectedRunConfig
from lr_reduction.models.results import CombinedReductionResult


def _config() -> ReductionConfig:
    return ReductionConfig(
        direct_beams={"db": DirectBeamConfig(db_runs=[11111])},
        runs=[ReflectedRunConfig(sequence_number=1, direct_beam="db", run_number=178200)],
    )


def test_reduce_autoreduce_executes_end_to_end(tmp_path, monkeypatch):
    """Runs the single-run leaf then FromDiskSequence, proving the incremental on-disk path (§6.4.4) works today."""
    monkeypatch.setattr("lr_reduction.io.config_loader.ConfigLoader.load", lambda _self, _path: _config())
    nexus_file = tmp_path / "REF_L_178200.nxs.h5"
    nexus_file.touch()

    result = reduce_autoreduce(str(nexus_file), str(tmp_path))

    assert isinstance(result, CombinedReductionResult)


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
