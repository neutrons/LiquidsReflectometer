from typing import TypeAlias, Union

from mantid.api import Workspace
from mantid.dataobjects import EventWorkspace

MantidWorkspace = Union[str, Workspace]

"""
The assembled direct beam produced by Op 1 (§2.2)
Consumed by single-run reduction (§3.2.1, §4.1.1)
"""
CompositeDirectBeam: TypeAlias = EventWorkspace
