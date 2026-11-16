"""Mantid algorithms of the reduction, registered with the AlgorithmFactory on import."""

from lr_reduction.operations.algorithms import (
    fit_beam_peak,  # noqa: F401
    merge_direct_beam_segments,  # noqa: F401
    trim_direct_beam_tof,  # noqa: F401
)
