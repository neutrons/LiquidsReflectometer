"""Autoreduction entrypoint (§6.4.1)

A thin composite in the autoreduce workflow.
Re-grows the sequence's on-disk assembly (§6.4.4, §2.4.1.a).
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from lr_reduction.api._shared import (
    locate_configuration_relative_to,
    parse_numeric_identifier,
    placeholder_sequence_result,
)
from lr_reduction.api.interfaces import Entrypoint
from lr_reduction.api.manual import reduce_run
from lr_reduction.config.model import ReductionConfig
from lr_reduction.io import ConfigLoader
from lr_reduction.io.orso import read_partials
from lr_reduction.io.orso import write as write_orso
from lr_reduction.io.report import html_report
from lr_reduction.types import SequenceResult
from lr_reduction.utils import get_logger

logger = get_logger(__name__)

_NEXUS_RUN_NUMBER_RE = re.compile(r"REF_L_(\d+)")


class FromDiskSequence(Entrypoint[SequenceResult]):
    """Assembles a sequence from its on-disk partial ORSO files (§2.4.1.a)."""

    def __init__(self, output_directory: str | Path, **overrides):
        self.output_directory = Path(output_directory)
        self.overrides = overrides
        self._config_loader = ConfigLoader()

    def load_configuration(self) -> ReductionConfig:
        config_path = locate_configuration_relative_to(self.output_directory)
        return self._config_loader.load(str(config_path))

    def load_data(self, config: ReductionConfig):
        sequence_id = parse_numeric_identifier(config.sequence_id, field_name="sequence_id")
        return read_partials(str(self.output_directory), sequence_id)

    def call_operations(self, config: ReductionConfig, _data) -> SequenceResult:
        # Placeholder until Op 3 (assemble-only) exists.
        reflected_runs = [reflected_run.run_id for reflected_run in config.reflected_runs]
        return placeholder_sequence_result(config, reflected_runs)

    def save_output(self, result: SequenceResult) -> None:
        result.html_report = html_report()
        write_orso(result, output_dir=str(self.output_directory))


def _run_number_from_nexus_path(nexus_file_path: Path) -> int:
    match = _NEXUS_RUN_NUMBER_RE.search(nexus_file_path.name)
    if not match:
        raise ValueError(f"Could not extract run number from nexus file path: {nexus_file_path}")
    return int(match.group(1))


def reduce_autoreduce(nexus_file_path: str | Path, output_directory_path: str | Path, **overrides) -> SequenceResult:
    """Autoreduction (on-disk assembly) (§6.4.1, §11.6.1).

    A thin wrapper (§6.4.3): reduces the newly-arrived run, writes its partial, then
    re-assembles the sequence from all partials on disk (§6.4.4).
    """
    nexus_file_path = Path(nexus_file_path)
    output_directory_path = Path(output_directory_path)
    run_number = _run_number_from_nexus_path(nexus_file_path)
    config_path = locate_configuration_relative_to(output_directory_path)

    logger.info(f"Autoreducing run {run_number} into {output_directory_path}")
    reduce_run(run_number, config_path, output_dir=str(output_directory_path), **overrides)
    return FromDiskSequence(output_directory_path, **overrides).execute()


def main(argv: list[str] | None = None) -> None:
    """CLI adapter for the autoreduction entrypoint (tier 3, §8); invoked by the BL4B autoreduction service."""
    parser = argparse.ArgumentParser(description="Autoreduction entrypoint for REF_L (§6.4.1).")
    parser.add_argument("nexus_file_path", type=Path)
    parser.add_argument("output_directory_path", type=Path)
    args = parser.parse_args(argv)
    reduce_autoreduce(args.nexus_file_path, args.output_directory_path)


if __name__ == "__main__":
    main()
