"""Batch driver entrypoint (§6.4.9)

Loops the autoreduction path over a contiguous run-number range for one IPTS experiment
as prototyped in scripts/shared/batch_reduce.py."""

from __future__ import annotations

import argparse
from pathlib import Path

from lr_reduction.api.autoreduce import reduce_autoreduce
from lr_reduction.utils import get_logger

logger = get_logger(__name__)


def reduce_batch(experiment_id: str, first_run: int, last_run: int, output_dir: str | Path, **overrides) -> None:
    """Batch driver (§6.4.9, §11.6.6).

    Re-invokes the single-run reduction (via the autoreduction path) for each existing
    run in [first_run, last_run], relying on on-disk incremental assembly (§2.4.1.a) to
    complete each sequence.
    """
    output_dir = Path(output_dir)
    for run_number in range(first_run, last_run + 1):
        nexus_file_path = Path("/SNS", "REF_L", experiment_id, "nexus", f"REF_L_{run_number}.nxs.h5")
        if not nexus_file_path.is_file():
            logger.warning(f"File does not exist: {nexus_file_path}")
            continue
        logger.info(f"Processing {nexus_file_path}")
        reduce_autoreduce(nexus_file_path, output_dir, **overrides)


def main(argv: list[str] | None = None) -> None:
    """CLI adapter for the batch driver (tier 3, §8)."""
    parser = argparse.ArgumentParser(
        description="Batch driver: autoreduce a contiguous run range for one IPTS (§6.4.9)."
    )
    parser.add_argument("experiment_id", type=str, help='Experiment id, e.g. "IPTS-20406".')
    parser.add_argument("first_run", type=int)
    parser.add_argument("last_run", type=int)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args(argv)
    reduce_batch(args.experiment_id, args.first_run, args.last_run, args.output_dir)


if __name__ == "__main__":
    main()
