"""Live-reduction entrypoint (§6.4.5): drives one in-memory reflected-run workspace,
re-reducing as more events arrive."""

from __future__ import annotations

import argparse

from mantid.dataobjects import EventWorkspace

from lr_reduction.api._shared import locate_standard_configuration
from lr_reduction.api._single_run import SingleRunReduction
from lr_reduction.api.autoreduce import FromDiskSequence
from lr_reduction.io.plotting import diagnostic_plot
from lr_reduction.models.config import ReductionConfig
from lr_reduction.models.results import CombinedReductionResult, ReductionResult


class LiveEntrypoint(SingleRunReduction):
    """Live reduction of one in-memory reflected-run workspace; the single-run leaf (§6.4.5).

    Sequence assembly is driven by `reduce_live`, exactly as `reduce_auto` drives
    `FromDiskSequence` around `AutoreduceSingleRun` (§6.4.4).
    """

    def __init__(self, reflected_run: EventWorkspace, **overrides):
        super().__init__(**overrides)
        self.reflected_run = reflected_run

    def load_configuration(self) -> ReductionConfig:
        # Configuration is discovered from the run number carried by the workspace (§6.4.5.1),
        # mapping to the same configuration the autoreduction entrypoint would apply to that run.
        run_number = self.reflected_run.getRunNumber()
        config_path = locate_standard_configuration(run_number)
        return self._config_loader.load(str(config_path))

    def load_data(self, _config: ReductionConfig) -> EventWorkspace:
        return self.reflected_run

    def publish(self, _result: ReductionResult) -> None:
        diagnostic_plot()


def reduce_live(reflected_run: EventWorkspace, **overrides) -> CombinedReductionResult:
    """Live reduction (§6.4.5, §11.6.3).

    Mirrors autoreduction (§6.4.4): reduces the run in-memory, writes its partial, then
    re-assembles the sequence from all partials on disk (§2.4.1.a).
    """
    entrypoint = LiveEntrypoint(reflected_run, **overrides)
    entrypoint.execute()
    return FromDiskSequence(entrypoint.output_directory, **overrides).execute()


def main(argv: list[str] | None = None) -> None:
    """CLI adapter placeholder (tier 3, §8).

    Live reduction is driven in-process by the live-data service (§1.3) with an
    in-memory event workspace; it has no standalone command-line invocation.
    """
    parser = argparse.ArgumentParser(
        description="Live reduction is driven by the live-data service and cannot be invoked from the command line."
    )
    parser.parse_args(argv)
    raise SystemExit("reduce_live requires an in-memory EventWorkspace supplied by the live-data service.")


if __name__ == "__main__":
    main()
