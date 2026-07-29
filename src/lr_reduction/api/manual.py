"""Manual entrypoints: single-run (§6.4.7) and full-sequence (§6.4.6) reduction."""

from __future__ import annotations

import argparse
from pathlib import Path

from lr_reduction.api._shared import placeholder_sequence_result, placeholder_single_run_result
from lr_reduction.api.interfaces import Entrypoint
from lr_reduction.io import ConfigLoader, RunData, RunLoader
from lr_reduction.io.orso import write as write_orso
from lr_reduction.io.report import html_report
from lr_reduction.models.config import ReductionConfig
from lr_reduction.models.results import CombinedReductionResult, ReductionResult


class ManualSingleRun(Entrypoint[ReductionResult]):
    """Manual reduction of a single run, by run number (§6.4.7, §11.6.5)."""

    def __init__(self, run_number: int, configuration: str | Path, **overrides):
        self.run_number = run_number
        self.configuration = Path(configuration)
        self.overrides = overrides
        self._config_loader = ConfigLoader()
        self._run_loader = RunLoader()

    def load_configuration(self) -> ReductionConfig:
        return self._config_loader.load(str(self.configuration))

    def load_data(self, _config: ReductionConfig) -> RunData:
        return self._run_loader.load(self.run_number)

    def call_operations(self, config: ReductionConfig, _data: RunData) -> ReductionResult:
        # Placeholder until Op 1 (build composite direct beam) / Op 2 (single-run reduce) exist.
        return placeholder_single_run_result(config)

    def save_output(self, result: ReductionResult) -> None:
        write_orso(result, output_dir=self.overrides.get("output_dir", "."))


class ManualRunSequence(Entrypoint[CombinedReductionResult]):
    """Manual reduction of a full run sequence, assembled in-memory (§6.4.6, §11.6.4)."""

    def __init__(self, sequence_id: int, configuration: str | Path, **overrides):
        self.sequence_id = sequence_id
        self.configuration = Path(configuration)
        self.overrides = overrides
        self._config_loader = ConfigLoader()
        self._run_loader = RunLoader()

    def load_configuration(self) -> ReductionConfig:
        return self._config_loader.load(str(self.configuration))

    def load_data(self, config: ReductionConfig) -> list[RunData]:
        return [
            self._run_loader.load(run_number)
            for reflected_run in config.runs
            for run_number in reflected_run.resolved_source_runs
        ]

    def call_operations(self, config: ReductionConfig, _data: list[RunData]) -> CombinedReductionResult:
        # Placeholder until Op 1 / Op 2 x N / Op 3 (in-memory assembly, §2.4.1.b) exist.
        return placeholder_sequence_result(config)

    def save_output(self, result: CombinedReductionResult) -> None:
        result.html_report = html_report()
        write_orso(result, output_dir=self.overrides.get("output_dir", "."))


def reduce_run(run_number: int, configuration: str | Path, **overrides) -> ReductionResult:
    """Manual single-run reduction (§6.4.7, §11.6.5)."""
    return ManualSingleRun(run_number, configuration, **overrides).execute()


def reduce_run_sequence(sequence_id: int, configuration: str | Path, **overrides) -> CombinedReductionResult:
    """Manual run-sequence reduction, in-memory assembly (§6.4.6, §11.6.4)."""
    return ManualRunSequence(sequence_id, configuration, **overrides).execute()


def main(argv: list[str] | None = None) -> None:
    """CLI adapter for the manual entrypoints (tier 3, §8)."""
    parser = argparse.ArgumentParser(description="Manual reduction: a single run, or a full run sequence.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Reduce a single run (§6.4.7).")
    run_parser.add_argument("run_number", type=int)
    run_parser.add_argument("configuration", type=Path)

    sequence_parser = subparsers.add_parser("sequence", help="Reduce a full run sequence (§6.4.6).")
    sequence_parser.add_argument("sequence_id", type=int)
    sequence_parser.add_argument("configuration", type=Path)

    args = parser.parse_args(argv)
    if args.command == "run":
        reduce_run(args.run_number, args.configuration)
    else:
        reduce_run_sequence(args.sequence_id, args.configuration)


if __name__ == "__main__":
    main()
