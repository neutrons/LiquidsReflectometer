from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple

from mantid.api import Workspace

if TYPE_CHECKING:
    from lr_reduction.models.run_data import RunData

type MantidWorkspace = str | Workspace


"""Alias for the name of a workspace registered in the analysis data service (§11.1.6)."""
type MantidWorkspaceName = str


"""Alias for ID numbers (run numbers, sequence id's, etc.)"""
type ID = int


class SingleReductionInput(NamedTuple):
    """A named tuple for the reduction input.

    `direct_beam_name` names the entry of `ReductionConfig.direct_beams` that `direct_beams` were
    loaded for, so `call_operations` can build that composite without re-resolving the reference.
    """

    run_data: RunData
    direct_beams: list[RunData]
    direct_beam_name: str
