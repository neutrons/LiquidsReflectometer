from pydantic import BaseModel


class DirectBeamConfig(BaseModel):
    """
    One composite direct-beam definition.
    """

    name: str
    db_runs: list[str]


class ReflectedRunConfig(BaseModel):
    """
    One reflected run referencing a named direct beam.
    """

    run_id: str
    direct_beam: str


class ReductionConfig(BaseModel):
    """
    Sequence-wide configuration.
    """

    sequence_id: str
    direct_beams: list[DirectBeamConfig]
    reflected_runs: list[ReflectedRunConfig]
