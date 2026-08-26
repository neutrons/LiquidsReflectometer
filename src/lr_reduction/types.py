from __future__ import annotations

from typing import TYPE_CHECKING, NamedTuple, TypeAlias, Union

from mantid.api import Workspace
from mantid.dataobjects import EventWorkspace

if TYPE_CHECKING:
    from lr_reduction.models.config import DirectBeamConfig
    from lr_reduction.models.run_data import RunData

MantidWorkspace = Union[str, Workspace]


"""The assembled direct beam to be used in single run reduction."""
CompositeDirectBeam: TypeAlias = EventWorkspace


"""Alias for ID numbers (run numbers, sequence id's, etc.)"""
ID: TypeAlias = int


class SingleReductionInput(NamedTuple):
    """A named tuple for the reduction input.

    `direct_beam_config` is the configuration `direct_beams` was resolved against (§3.3.2), carried
    alongside so `call_operations` can reuse it directly instead of re-deriving it from `ReductionConfig`.
    """

    run_data: RunData | EventWorkspace
    direct_beams: list[RunData]
    direct_beam_config: DirectBeamConfig
