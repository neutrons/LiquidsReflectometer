"""Autoreduction entrypoint (§6.4.1)

A thin composite in the autoreduce workflow.
Re-grows the sequence's on-disk assembly (§6.4.4, §2.4.1.a).
"""

from __future__ import annotations

import argparse
from pathlib import Path

from lr_reduction.api._shared import locate_configuration_relative_to
from lr_reduction.api._single_run import SingleRunReduction
from lr_reduction.api.interfaces import Entrypoint
from lr_reduction.io import ConfigLoader, RunData, RunLoader
from lr_reduction.io.orso import read_partials, write_orso
from lr_reduction.io.report import html_report
from lr_reduction.models.config import ReductionConfig
from lr_reduction.models.results import CombinedReductionResult, ReductionResult
from lr_reduction.operations import CombineResultsOperation
from lr_reduction.utils import get_logger

logger = get_logger(__name__)


class AutoreduceSingleRun(SingleRunReduction):
    """Single run from a NeXus path, configuration relative to the output directory (§6.4.1)."""

    def __init__(self, nexus_file_path: str | Path, output_directory: str | Path, **overrides):
        super().__init__(**overrides)
        self.nexus_file_path = Path(nexus_file_path)
        self._output_directory = Path(output_directory)
        self._run_loader = RunLoader()

    def load_configuration(self) -> ReductionConfig:
        config_path = locate_configuration_relative_to(self._output_directory)
        return self._config_loader.load(str(config_path))

    def load_data(self, _config: ReductionConfig) -> RunData:
        return self._run_loader.load_from_path(self.nexus_file_path)

    @property
    def output_directory(self) -> Path:
        return self._output_directory


class FromDiskSequence(Entrypoint[CombinedReductionResult]):
    """Assembles a sequence from its on-disk partial ORSO files (§2.4.1.a)."""

    def __init__(self, output_directory: str | Path, **overrides):
        self.output_directory = Path(output_directory)
        self.overrides = overrides
        self._config_loader = ConfigLoader()

    def load_configuration(self) -> ReductionConfig:
        config_path = locate_configuration_relative_to(self.output_directory)
        return self._config_loader.load(str(config_path))

    def load_data(self, config: ReductionConfig):  # noqa: ARG002
        # Placeholder: ReductionConfig does not yet carry a sequence_id (§3.3.1); a real
        # identifier will come from the acquired run's metadata (§3.1.2.1) once available.
        sequence_id = 0
        return read_partials(str(self.output_directory), sequence_id)

    def call_operations(self, data: list[ReductionResult], config: ReductionConfig) -> CombinedReductionResult:
        op = CombineResultsOperation(data, config)
        op.validate_input()
        result = op.process()
        op.cleanup()
        return result

    def save_output(self, result: CombinedReductionResult) -> None:
        result.html_report = html_report()
        write_orso(result, output_dir=str(self.output_directory))


def reduce_auto(nexus_file_path: str | Path, output_directory_path: str | Path, **overrides) -> CombinedReductionResult:
    """Autoreduction (on-disk assembly) (§6.4.1, §11.6.1).

    A thin wrapper (§6.4.3): reduces the newly-arrived run, writes its partial, then
    re-assembles the sequence from all partials on disk (§6.4.4).
    """
    nexus_file_path = Path(nexus_file_path)
    output_directory_path = Path(output_directory_path)
    logger.info(f"Autoreducing {nexus_file_path.name} into {output_directory_path}")
    AutoreduceSingleRun(nexus_file_path, output_directory_path, **overrides).execute()
    return FromDiskSequence(output_directory_path, **overrides).execute()


def main(argv: list[str] | None = None) -> None:
    """CLI adapter for the autoreduction entrypoint (tier 3, §8); invoked by the BL4B autoreduction service."""
    parser = argparse.ArgumentParser(description="Autoreduction entrypoint for REF_L (§6.4.1).")
    parser.add_argument("nexus_file_path", type=Path)
    parser.add_argument("output_directory_path", type=Path)
    args = parser.parse_args(argv)
    reduce_auto(args.nexus_file_path, args.output_directory_path)


if __name__ == "__main__":
    main()
