"""Registration of this package's Mantid algorithms.

`AlgorithmFactory.subscribe` alone makes an algorithm available to `AlgorithmManager`, but not
as a `mantid.simpleapi` function: `simpleapi` builds its functions from the algorithms known when
it is first imported, and every algorithm module imports `simpleapi` before it subscribes. Mantid
fills that gap for algorithms loaded from a configured plugin directory; for algorithms shipped in
a package, its documented alternative is to create the function at runtime with
`mantid.simpleapi._create_algorithm_function`, as `register` does
(https://developer.mantidproject.org/PythonAlgorithmsInExternalProjects.html).
"""

import mantid.simpleapi as simpleapi
from mantid.api import AlgorithmFactory, PythonAlgorithm


def register(algorithm_class: type[PythonAlgorithm]) -> None:
    """Subscribe `algorithm_class` and give it a `mantid.simpleapi` function, as Mantid's own algorithms have.

    Called once per algorithm, at the bottom of the module that defines it, so importing the
    module is all a caller needs to do before calling the algorithm through `mantid.simpleapi`.

    Parameters
    ----------
    algorithm_class
        The algorithm to register.
    """
    AlgorithmFactory.subscribe(algorithm_class)
    algorithm = algorithm_class()
    algorithm.initialize()
    simpleapi._create_algorithm_function(algorithm.name(), algorithm.version(), algorithm)
