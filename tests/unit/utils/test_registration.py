"""Unit tests for `lr_reduction.utils.registration`."""

import subprocess
import sys
import textwrap
from pathlib import Path

import mantid.simpleapi as simpleapi
import pytest
from mantid.api import AlgorithmFactory, PythonAlgorithm
from mantid.kernel import Direction

import lr_reduction.operations.algorithms as algorithms_package
from lr_reduction.utils.registration import register

ALGORITHMS_DIR = Path(algorithms_package.__file__).parent


class LrRegistrationTestAlgorithm(PythonAlgorithm):
    """Doubles `Value`; refuses a negative one from PyExec, where only a rethrowing call raises."""

    def category(self):
        return "Testing"

    def PyInit(self):  # noqa: N802
        self.declareProperty("Value", 0.0)
        self.declareProperty("Doubled", 0.0, direction=Direction.Output)

    def PyExec(self):  # noqa: N802
        value = self.getProperty("Value").value
        if value < 0:
            raise ValueError("Value must not be negative")
        self.setProperty("Doubled", 2 * value)


@pytest.fixture(scope="module")
def registered():
    register(LrRegistrationTestAlgorithm)
    yield "LrRegistrationTestAlgorithm"
    AlgorithmFactory.unsubscribe("LrRegistrationTestAlgorithm", 1)


def test_register_subscribes_the_algorithm(registered):
    assert AlgorithmFactory.exists(registered)


def test_register_creates_a_simpleapi_function(registered):
    assert getattr(simpleapi, registered)(Value=1.5) == 3.0


def test_simpleapi_call_raises_on_a_failure_inside_pyexec(registered):
    with pytest.raises(RuntimeError, match="must not be negative"):
        getattr(simpleapi, registered)(Value=-1.0)


def test_every_algorithm_module_registers_with_simpleapi():
    """Import each module in a fresh interpreter without plugin directories, as on a machine whose
    Mantid configuration does not point at this package, then check `simpleapi` for its algorithms."""
    modules = sorted(path.stem for path in ALGORITHMS_DIR.glob("*.py") if path.stem != "__init__")
    script = textwrap.dedent(
        f"""
        from mantid.kernel import config
        config["user.python.plugins.directories"] = ""
        import importlib
        import inspect
        import mantid.simpleapi as simpleapi
        from mantid.api import PythonAlgorithm

        missing = []
        for module_name in {modules!r}:
            module = importlib.import_module(f"lr_reduction.operations.algorithms.{{module_name}}")
            for name, member in inspect.getmembers(module, inspect.isclass):
                if issubclass(member, PythonAlgorithm) and member.__module__ == module.__name__:
                    if not hasattr(simpleapi, member().name()):
                        missing.append(f"{{module_name}}.{{name}}")
        print("missing:", missing)
        """
    )
    completed = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, timeout=300)
    assert completed.returncode == 0, completed.stderr
    assert "missing: []" in completed.stdout, completed.stdout
