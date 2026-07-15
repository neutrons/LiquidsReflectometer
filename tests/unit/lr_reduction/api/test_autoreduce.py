from lr_reduction.api.autoreduce import FromDiskSequence, reduce_autoreduce
from lr_reduction.types import SequenceResult


def test_reduce_autoreduce_executes_end_to_end(tmp_path):
    """Runs the single-run leaf then FromDiskSequence, proving the incremental on-disk path (§6.4.4) works today."""
    nexus_file = tmp_path / "REF_L_178200.nxs.h5"
    nexus_file.touch()

    result = reduce_autoreduce(str(nexus_file), str(tmp_path))

    assert isinstance(result, SequenceResult)


def test_from_disk_sequence_reuses_read_partials(tmp_path):
    result = FromDiskSequence(tmp_path).execute()
    assert isinstance(result, SequenceResult)


def test_from_disk_sequence_converts_sequence_id_to_int(tmp_path, monkeypatch):
    captured_sequence_ids = []

    def _capture_sequence_id(_partial_dir, sequence_id):
        captured_sequence_ids.append(sequence_id)
        return []

    monkeypatch.setattr("lr_reduction.api.autoreduce.read_partials", _capture_sequence_id)

    FromDiskSequence(tmp_path).execute()

    assert captured_sequence_ids == [0]
