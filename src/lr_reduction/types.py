from typing import TypeAlias, Union

from mantid.api import Workspace
from mantid.dataobjects import EventWorkspace

MantidWorkspace = Union[str, Workspace]

"""The assembled direct beam to be used in single run reduction."""
CompositeDirectBeam: TypeAlias = EventWorkspace


"""Alias for ID numbers (run numbers, sequence id's, etc.)"""
ID: TypeAlias = int
