"""Shared machinery for the single-run entrypoint leaves (§2.3).

Manual, live, and autoreduce single-run reduction differ only in where their
configuration and data come from (load_configuration / load_data); this base
carries the identical call_operations / save_output steps.
"""

from __future__ import annotations

from pathlib import Path

from lr_reduction.api.interfaces import Entrypoint
from lr_reduction.exceptions import LrReductionError
from lr_reduction.io import ConfigLoader
from lr_reduction.io.orso import write_orso
from lr_reduction.models import ReductionConfig, ReductionResult
from lr_reduction.operations import DirectBeamCompositionOperation, SingleRunReductionOperation
from lr_reduction.types import ID, SingleReductionInput


def reduce_one(data: SingleReductionInput, config: ReductionConfig, sequence_number: ID) -> ReductionResult:
    """Perform a single-run reduction."""
    db_op = DirectBeamCompositionOperation(data=data.direct_beams, config=data.direct_beam_config)
    comp_db = db_op.execute()
    op = SingleRunReductionOperation(data.run_data, config, comp_db, sequence_number)
    result = op.execute()
    return result


class SingleRunReduction(Entrypoint[SingleReductionInput, ReductionResult]):
    """Shared machinery for the single-run surfaces (§2.3).

    Leaves differ only in load_configuration / load_data.
    """

    sequence_number: ID | None

    def __init__(self, **overrides):
        """Initialize the single-run reduction entrypoint."""
        self.overrides = overrides
        self._config_loader = ConfigLoader()
        # self.sequence_number = self.overrides.get("sequence_number", None)  # Maybe?

    def call_operations(self, data: SingleReductionInput, config: ReductionConfig) -> ReductionResult:
        """Perform the single-run reduction operations."""
        if self.sequence_number is None:
            raise LrReductionError("sequence_number must be set before calling operations")
        return reduce_one(data, config, self.sequence_number)

    def save_output(self, result: ReductionResult) -> None:
        """Save the reduction output."""
        write_orso(result, output_dir=str(self.output_directory))
        self.publish(result)

    def publish(self, _result: ReductionResult) -> None:
        """Emission beyond the ORSO partial. Default: nothing; live overrides it."""

    @property
    def output_directory(self) -> Path:
        return Path(self.overrides.get("output_dir", "."))

    # load_configuration / load_data remain abstract (inherited from Entrypoint)
