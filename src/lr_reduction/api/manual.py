"""Manual entrypoints: single-run (§6.4.7) and full-sequence (§6.4.6) reduction."""

from __future__ import annotations

import argparse
from pathlib import Path

from lr_reduction.api._single_run import SingleRunReduction
from lr_reduction.api.interfaces import Entrypoint
from lr_reduction.io import ConfigLoader, RunData, RunLoader
from lr_reduction.io.orso import write_orso
from lr_reduction.io.report import html_report
from lr_reduction.models.config import ReductionConfig
from lr_reduction.models.results import CombinedReductionResult, ReductionResult
from lr_reduction.operations import CombineResultsOperation, DirectBeamCompositionOperation, SingleRunReductionOperation
from lr_reduction.types import CompositeDirectBeam


class ManualSingleRun(SingleRunReduction):
    """Manual reduction of a single run, by run number."""

    def __init__(self, run_number: int, configuration: str | Path, **overrides):
        super().__init__(**overrides)
        self.run_number = run_number
        self.configuration = Path(configuration)
        self._run_loader = RunLoader()

    def load_configuration(self) -> ReductionConfig:
        return self._config_loader.load(str(self.configuration))

    def load_data(self, _config: ReductionConfig) -> RunData:
        return self._run_loader.load(self.run_number)


class ManualRunSequence(Entrypoint[list[RunData], CombinedReductionResult]):
    """Manual reduction of a full run sequence, assembled in-memory."""

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

    def call_operations(self, data: list[RunData], config: ReductionConfig) -> CombinedReductionResult:
        _results: list[ReductionResult] = []
        _comp_dbs: dict[str, CompositeDirectBeam] = {}
        for idx, d in enumerate(data):
            # Compose the direct beam for this run
            db_name = config.runs[idx].direct_beam
            db_op = DirectBeamCompositionOperation(data=[d], config=config.direct_beams[db_name])
            comp_db = db_op.execute()
            _comp_dbs[db_name] = comp_db

            # Perform the single-run reduction for this run
            single_run_config = config.for_run(config.runs[idx].sequence_number)
            op = SingleRunReductionOperation(d, single_run_config, comp_db)
            op.validate_input()
            result = op.process()
            op.cleanup()
            _results.append(result)

        # Combine the results
        combine_op = CombineResultsOperation(_results, config)
        combine_op.validate_input()
        combined_result = combine_op.process()
        combine_op.cleanup()
        return combined_result

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
