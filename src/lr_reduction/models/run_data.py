from __future__ import annotations

from dataclasses import dataclass

from lr_reduction.types import ID, MantidWorkspace


@dataclass
class RunData:
    """
    Raw neutron event data for a single run.

    Returned by a RunLoaderInterface instance.
    """

    run_number: ID = 123
    sequence_number: ID = 1
    workspace: MantidWorkspace = ""  # Placeholder: actual workspace will be loaded at runtime
