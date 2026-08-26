"""Manual entrypoints: single-run (§6.4.7) and full-sequence (§6.4.6) reduction."""

from __future__ import annotations

import argparse
from pathlib import Path

from lr_reduction.api._shared import get_direct_beam_config
from lr_reduction.api._single_run import SingleRunReduction, reduce_one
from lr_reduction.api.interfaces import Entrypoint
from lr_reduction.exceptions import LrValidationError
from lr_reduction.io import ConfigLoader, RunLoader
from lr_reduction.io.orso import write_orso
from lr_reduction.io.report import html_report
from lr_reduction.models.config import ReductionConfig
from lr_reduction.models.results import CombinedReductionResult, ReductionResult
from lr_reduction.operations import CombineResultsOperation
from lr_reduction.types import ID, SingleReductionInput


class ManualSingleRun(SingleRunReduction):
    """Manual reduction of a single run, by run number."""

    def __init__(self, run_number: ID, configuration: str | Path, sequence_number: ID | None = None, **overrides):
        super().__init__(**overrides)
        self.configuration = Path(configuration)
        self.run_number = run_number
        self.sequence_number = sequence_number
        self._run_loader = RunLoader()

    def load_configuration(self) -> ReductionConfig:
        return self._config_loader.load(str(self.configuration))

    def load_data(self, config: ReductionConfig) -> SingleReductionInput:
        if self.sequence_number is None:
            # TODO: if sequence_number is None, use sequence_number from workspace sample logs
            #       and assign to self.sequence_number for use in parent call_operations
            self.sequence_number = 1  # Placeholder: sequence_number is not yet available
        run = self._run_loader.load(self.run_number)
        db_config = get_direct_beam_config(self.sequence_number, config)
        direct_beams = [self._run_loader.load(run_number) for run_number in db_config.direct_beam_run_numbers]
        return SingleReductionInput(run_data=run, direct_beams=direct_beams, direct_beam_config=db_config)


class ManualRunSequence(Entrypoint[list[SingleReductionInput], CombinedReductionResult]):
    """Manual reduction of a full run sequence, assembled in-memory.

    NOTE: `sequence_numbers` as an override is currently not implemented
    """

    def __init__(
        self, run_numbers: list[ID], configuration: str | Path, sequence_numbers: list[ID] | None = None, **overrides
    ):
        self.run_numbers = run_numbers
        self.configuration = Path(configuration)
        self.sequence_numbers = sequence_numbers
        self.overrides = overrides
        self._config_loader = ConfigLoader()
        self._run_loader = RunLoader()

        # Validate_sequence_numbers
        if self.sequence_numbers is not None:
            if len(self.sequence_numbers) != len(self.run_numbers):
                raise LrValidationError(
                    f"Override sequences numbers {self.sequence_numbers}"
                    + f" do not match the number of runs in the configuration ({len(self.run_numbers)})."
                )

    def load_configuration(self) -> ReductionConfig:
        return self._config_loader.load(str(self.configuration))

    def load_data(self, config: ReductionConfig) -> list[SingleReductionInput]:
        run_data = []
        for run_number in self.run_numbers:
            run = self._run_loader.load(run_number)
            db_config = get_direct_beam_config(run.sequence_number, config)
            run_data.append(
                SingleReductionInput(
                    run_data=run,
                    direct_beams=[self._run_loader.load(db) for db in db_config.direct_beam_run_numbers],
                    direct_beam_config=db_config,
                )
            )
        return run_data

    # TODO: figure out how to address override sequence numbers
    def call_operations(self, data: list[SingleReductionInput], config: ReductionConfig) -> CombinedReductionResult:
        """Perform the run-sequence reduction operations."""
        # for each run, compose the direct beam and reduce the run
        reduced_runs = []
        for run in data:
            reduced_run = reduce_one(run, config, run.run_data.sequence_number)
            reduced_runs.append(reduced_run)
        # combine the reduced runs into a single result
        combine_op = CombineResultsOperation(reduced_runs, config)
        combined_result = combine_op.execute()
        return combined_result

    def save_output(self, result: CombinedReductionResult) -> None:
        result.html_report = html_report()
        write_orso(result, output_dir=self.overrides.get("output_dir", "."))


def reduce_run(
    run_number: int, configuration: str | Path, sequence_number: ID | None = None, **overrides
) -> ReductionResult:
    """Manual single-run reduction (§6.4.7, §11.6.5)."""
    return ManualSingleRun(run_number, configuration, sequence_number, **overrides).execute()


def reduce_and_combine_runs(
    run_numbers: list[ID], configuration: str | Path, sequence_numbers: list[ID] | None = None, **overrides
) -> CombinedReductionResult:
    """Manual run-sequence reduction, in-memory assembly (§6.4.6, §11.6.4)."""
    return ManualRunSequence(run_numbers, configuration, sequence_numbers, **overrides).execute()


def main(argv: list[str] | None = None) -> None:
    """CLI adapter for the manual entrypoints (tier 3, §8)."""
    parser = argparse.ArgumentParser(description="Manual reduction: a single run, or a full run sequence.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Reduce a single run (§6.4.7).")
    run_parser.add_argument("run_number", type=int)
    run_parser.add_argument("--configuration", type=Path, required=True)

    sequence_parser = subparsers.add_parser("sequence", help="Reduce a full run sequence.")
    sequence_parser.add_argument("--run-numbers", type=int, nargs="+", required=True)
    sequence_parser.add_argument("--configuration", type=Path, required=True)
    sequence_parser.add_argument("--sequence-numbers", type=int, nargs="*", default=None)

    args = parser.parse_args(argv)
    if args.command == "run":
        reduce_run(args.run_number, args.configuration)
    else:
        reduce_and_combine_runs(args.run_numbers, args.configuration, sequence_numbers=args.sequence_numbers)


if __name__ == "__main__":
    main()
