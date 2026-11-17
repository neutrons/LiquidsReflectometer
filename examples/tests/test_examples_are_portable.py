"""
What keeps an example runnable anywhere, checked without data.

An example names no folder of its own: inputs come from example_paths.data_folder(), results go to
example_paths.output_folder(). Importing an example runs nothing, and no example function is named test_*,
which pytest would collect and run against whatever folder it named.
"""

import ast
import os
import subprocess
import sys
from pathlib import Path

import pytest

import example_paths

EXAMPLES = Path(example_paths.__file__).parent
FILES = sorted(EXAMPLES.glob("example_*.py"))


@pytest.mark.parametrize("path", FILES, ids=[path.name for path in FILES])
def test_an_example_names_no_folder_of_its_own(path):
    absolute = [
        (node.lineno, node.value)
        for node in ast.walk(ast.parse(path.read_text()))
        if isinstance(node, ast.Constant) and isinstance(node.value, str) and node.value.lstrip().startswith(("/", "~"))
        and "\n" not in node.value
    ]
    assert absolute == []


@pytest.mark.parametrize("path", FILES, ids=[path.name for path in FILES])
def test_no_example_function_is_named_as_a_test(path):
    names = [node.name for node in ast.walk(ast.parse(path.read_text())) if isinstance(node, ast.FunctionDef)]
    assert [name for name in names if name.startswith("test")] == []


@pytest.mark.parametrize("path", FILES, ids=[path.name for path in FILES])
def test_importing_an_example_runs_nothing(path, tmp_path):
    """Imported where there is no data folder and no way to write: an example that ran on import would fail."""
    environment = {**os.environ, example_paths.DATA_VARIABLE: str(tmp_path / "no data here"),
                   example_paths.OUTPUT_VARIABLE: str(tmp_path / "output"), "MPLBACKEND": "Agg"}
    done = subprocess.run([sys.executable, "-c", f"import {path.stem}"], cwd=EXAMPLES, env=environment,
                          capture_output=True, text=True, check=False)
    assert done.returncode == 0, done.stderr
    assert not (tmp_path / "output").exists()


def test_the_data_folder_is_the_one_chosen(tmp_path, monkeypatch):
    monkeypatch.setenv(example_paths.DATA_VARIABLE, str(tmp_path))
    assert example_paths.data_folder() == tmp_path


def test_a_data_folder_that_is_not_there_says_how_to_get_one(tmp_path, monkeypatch):
    monkeypatch.setenv(example_paths.DATA_VARIABLE, str(tmp_path / "nowhere"))
    with pytest.raises(FileNotFoundError, match=example_paths.DATA_VARIABLE):
        example_paths.data_folder()


def test_the_output_folder_is_created(tmp_path, monkeypatch):
    monkeypatch.setenv(example_paths.OUTPUT_VARIABLE, str(tmp_path / "a" / "b"))
    assert example_paths.output_folder() == tmp_path / "a" / "b"
    assert (tmp_path / "a" / "b").is_dir()


def test_the_output_folder_defaults_to_one_below_the_working_folder(tmp_path, monkeypatch):
    monkeypatch.delenv(example_paths.OUTPUT_VARIABLE, raising=False)
    monkeypatch.chdir(tmp_path)
    assert example_paths.output_folder() == tmp_path / example_paths.DEFAULT_OUTPUT_NAME


def test_named_examples_run_in_the_order_named_and_all_run_when_none_is_named():
    ran = []

    def example_a():
        ran.append("a")

    def example_b():
        ran.append("b")

    example_paths.run_examples([example_a, example_b], arguments=["example_b"])
    example_paths.run_examples([example_a, example_b], arguments=[])
    assert ran == ["b", "a", "b"]


def test_an_unknown_example_name_lists_the_ones_there_are():
    def example_a():
        pass

    with pytest.raises(SystemExit, match="example_a"):
        example_paths.run_examples([example_a], arguments=["example_z"])
