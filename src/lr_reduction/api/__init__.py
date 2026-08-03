"""Public entrypoint API (§6.4, §8)."""

from lr_reduction.api.autoreduce import reduce_auto
from lr_reduction.api.batch import reduce_batch
from lr_reduction.api.live import reduce_live
from lr_reduction.api.manual import reduce_run, reduce_run_sequence

__all__ = [
    "reduce_auto",
    "reduce_batch",
    "reduce_live",
    "reduce_run",
    "reduce_run_sequence",
]
