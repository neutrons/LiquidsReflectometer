from lr_reduction.api.batch import reduce_batch


def test_reduce_batch_skips_missing_runs(tmp_path):
    """Non-existent nexus files are skipped rather than raising (mirrors scripts/shared/batch_reduce.py)."""
    reduce_batch("IPTS-00000", 1, 3, tmp_path)
