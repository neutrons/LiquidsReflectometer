"""Shared machinery for the single-run entrypoint leaves (§2.3).

Manual, live, and autoreduce single-run reduction differ only in where their
configuration and data come from (load_configuration / load_data); this base
carries the identical call_operations / save_output steps.
"""

from __future__ import annotations

from pathlib import Path

from lr_reduction.api._shared import placeholder_single_run_result
from lr_reduction.api.interfaces import Entrypoint
from lr_reduction.io import ConfigLoader
from lr_reduction.io.orso import write as write_orso
from lr_reduction.models.config import ReductionConfig
from lr_reduction.models.results import ReductionResult


class SingleRunReduction(Entrypoint[ReductionResult]):
    """Shared machinery for the single-run surfaces (§2.3).

    Leaves differ only in load_configuration / load_data.
    """

    def __init__(self, **overrides):
        self.overrides = overrides
        self._config_loader = ConfigLoader()

    def call_operations(self, config: ReductionConfig, _data) -> ReductionResult:
        # Placeholder until Op 1 (composite direct beam) / Op 2 (single-run reduce) exist (§1.2).
        return placeholder_single_run_result(config)

    def save_output(self, result: ReductionResult) -> None:
        write_orso(result, output_dir=str(self.output_directory))
        self.publish(result)

    def publish(self, _result: ReductionResult) -> None:
        """Emission beyond the ORSO partial. Default: nothing; live overrides it."""

    @property
    def output_directory(self) -> Path:
        return Path(self.overrides.get("output_dir", "."))

    # load_configuration / load_data remain abstract (inherited from Entrypoint)
