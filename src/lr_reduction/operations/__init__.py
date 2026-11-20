# Importing the algorithms package registers every reduction algorithm with `mantid.simpleapi`,
# so each operation can call them through it.
import lr_reduction.operations.algorithms  # noqa: F401
from lr_reduction.operations.combine import CombineResultsOperation
from lr_reduction.operations.direct_beam import DirectBeamCompositionOperation
from lr_reduction.operations.single_run import SingleRunReductionOperation

__all__ = [
    "CombineResultsOperation",
    "DirectBeamCompositionOperation",
    "SingleRunReductionOperation",
]
