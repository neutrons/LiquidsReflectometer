from __future__ import annotations

from dataclasses import dataclass

from lr_reduction.types import ID


@dataclass
class RunData:
    """
    Raw neutron event data for a single run.

    Returned by a RunLoaderInterface instance.
    """

    run_number: ID = 123
