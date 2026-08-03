from __future__ import annotations

from dataclasses import dataclass

from lr_reduction.types import MantidWorkspace


@dataclass
class RunData:
    """
    Raw neutron event data for a single run.

    Returned by a RunLoaderInterface instance.
    """

    data: MantidWorkspace
