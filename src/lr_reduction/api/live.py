"""Live-reduction entrypoint (§6.4.5): drives one in-memory reflected-run workspace,
re-reducing as more events arrive."""

from __future__ import annotations

import argparse

from mantid.dataobjects import EventWorkspace

from lr_reduction.api._shared import locate_standard_configuration, placeholder_single_run_result
from lr_reduction.api.interfaces import Entrypoint
from lr_reduction.io import ConfigLoader
from lr_reduction.io.orso import write as write_orso
from lr_reduction.io.plotting import diagnostic_plot
from lr_reduction.models.config import ReductionConfig
from lr_reduction.models.results import ReductionResult


class LiveEntrypoint(Entrypoint[ReductionResult]):
    """Live reduction of one in-memory reflected-run workspace (§6.4.5)."""

    def __init__(self, reflected_run: EventWorkspace, **overrides):
        self.reflected_run = reflected_run
        self.overrides = overrides
        self._config_loader = ConfigLoader()

    def load_configuration(self) -> ReductionConfig:
        # Configuration is discovered from the run number carried by the workspace (§6.4.5.1),
        # mapping to the same configuration the autoreduction entrypoint would apply to that run.
        run_number = self.reflected_run.getRunNumber()
        config_path = locate_standard_configuration(run_number)
        return self._config_loader.load(str(config_path))

    def load_data(self, _config: ReductionConfig) -> EventWorkspace:
        return self.reflected_run

    def call_operations(self, config: ReductionConfig, _data: EventWorkspace) -> ReductionResult:
        # Placeholder until Op 1 (build composite direct beam) / Op 2 (single-run reduce) exist.
        return placeholder_single_run_result(config)

    def save_output(self, result: ReductionResult) -> None:
        write_orso(result, output_dir=self.overrides.get("output_dir", "."))
        diagnostic_plot()


def reduce_live(reflected_run: EventWorkspace, **overrides) -> ReductionResult:
    """Live reduction (in-memory single run) (§6.4.5, §11.6.3)."""
    return LiveEntrypoint(reflected_run, **overrides).execute()


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
