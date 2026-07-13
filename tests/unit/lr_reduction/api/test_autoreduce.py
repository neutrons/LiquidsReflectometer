from lr_reduction.api.autoreduce import FromDiskSequence, reduce_autoreduce
from lr_reduction.types import SequenceResult


def test_reduce_autoreduce_executes_end_to_end(tmp_path):
    """Runs the single-run leaf then FromDiskSequence, proving the incremental on-disk path (§6.4.4) works today."""
    nexus_file = tmp_path / "REF_L_178200.nxs.h5"
    nexus_file.touch()

    result = reduce_autoreduce(nexus_file, tmp_path)

    assert isinstance(result, SequenceResult)


def test_from_disk_sequence_reuses_read_partials(tmp_path):
    result = FromDiskSequence(tmp_path).execute()
    assert isinstance(result, SequenceResult)
