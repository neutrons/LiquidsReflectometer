from pydantic import BaseModel

from lr_reduction.types import ID


class AssemblyConfig(BaseModel):
    """Configuration for combining reduced data from multiple runs."""

    ...


class DirectBeamConfig(BaseModel):
    """
    One composite direct-beam definition.
    """

    name: str
    db_runs: list[ID]


class ReflectedRunConfig(BaseModel):
    """
    One reflected run referencing a named direct beam.
    """

    run_id: ID
    direct_beam: str


class ReductionConfig(BaseModel):
    """
    Sequence-wide configuration.
    """

    sequence_id: ID
    direct_beams: list[DirectBeamConfig]
    reflected_runs: list[ReflectedRunConfig]
