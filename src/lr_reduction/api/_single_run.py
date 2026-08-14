"""Shared machinery for the single-run entrypoint leaves (§2.3).

Manual, live, and autoreduce single-run reduction differ only in where their
configuration and data come from (load_configuration / load_data); this base
carries the identical call_operations / save_output steps.
"""

from __future__ import annotations

from pathlib import Path

from lr_reduction.api.interfaces import Entrypoint
from lr_reduction.io import ConfigLoader
from lr_reduction.io.orso import write_orso
from lr_reduction.models import ReductionConfig, ReductionResult, RunData
from lr_reduction.operations import DirectBeamCompositionOperation, SingleRunReductionOperation


class SingleRunReduction(Entrypoint[RunData, ReductionResult]):
    """Shared machinery for the single-run surfaces (§2.3).

    Leaves differ only in load_configuration / load_data.
    """

    def __init__(self, **overrides):
        self.overrides = overrides
        self._config_loader = ConfigLoader()

    def call_operations(self, data: RunData, config: ReductionConfig) -> ReductionResult:
        """Perform the single-run reduction operations."""
        db_name= config.runs[0].direct_beam
        db_op = DirectBeamCompositionOperation(
            data=[data],
            config = config.direct_beams[db_name],
        )
        comp_db = db_op.execute()
        single_run_config = config.for_run(config.runs[0].sequence_number)
        op = SingleRunReductionOperation(data, single_run_config, comp_db)
        result = op.execute()
        return result

    def save_output(self, result: ReductionResult) -> None:
        write_orso(result, output_dir=str(self.output_directory))
        self.publish(result)

    def publish(self, _result: ReductionResult) -> None:
        """Emission beyond the ORSO partial. Default: nothing; live overrides it."""

    @property
    def output_directory(self) -> Path:
        return Path(self.overrides.get("output_dir", "."))

    # load_configuration / load_data remain abstract (inherited from Entrypoint)
