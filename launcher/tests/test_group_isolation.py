"""The launcher suite runs each process-risk group in a pytest process of its own (conftest.py).

Three properties, each through real pytest processes:
- the census: the three groups hold every test of the suite once, none twice, none lost;
- a native crash in one group fails that group, with the test and the signal in the output, while the other groups
  run and report;
- a hang in one group is stopped by its timeout, with the stacks dumped, while the other groups run and report.

The two failure cases run on a scratch tree laid out as the repository is (its pyproject.toml and this conftest.py at
launcher/tests/), so the orchestration is on there exactly as in the repository's own session.
"""

import ast
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
GROUP_VARIABLE = "LR_LAUNCHER_TEST_GROUP"
RECORD_VARIABLE = "LR_LAUNCHER_TEST_RECORD"


def _environment(**extra):
    """This process's environment without the group's variables (a test that runs in a group's process inherits them)."""
    environment = {key: value for key, value in os.environ.items() if key not in (GROUP_VARIABLE, RECORD_VARIABLE)}
    environment.update(extra)
    return environment


def _collected(group=None):
    environment = _environment(**({GROUP_VARIABLE: group} if group else {}))
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "--collect-only", "-q", "-p", "no:cacheprovider", "launcher/tests/"],
        cwd=ROOT, env=environment, capture_output=True, text=True, timeout=600, check=False,
    )
    assert result.returncode == 0, result.stdout[-2000:] + result.stderr[-2000:]
    return [line for line in result.stdout.splitlines() if "::" in line]


@pytest.mark.timeout(900)
def test_the_groups_hold_every_test_once():
    """T2: the union of the three groups' tests is the whole suite's, with no test in two groups and none in none."""
    everything = _collected()
    groups = {group: _collected(group) for group in ("model", "data", "ui")}
    together = [test for tests in groups.values() for test in tests]
    assert sorted(together) == sorted(everything)
    assert len(together) == len(set(together))
    assert all(groups.values()), {group: len(tests) for group, tests in groups.items()}


def test_every_real_data_test_in_the_model_module_is_marked():
    """A test of test_roi_estimate.py that reads the real event files (the `nexus_dir` fixture) carries `realdata`, so
    it runs in the data group, not with the Qt-free model tests."""
    tree = ast.parse((HERE / "test_roi_estimate.py").read_text(encoding="utf-8"))
    unmarked = [
        node.name
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name.startswith("test_")
        and "nexus_dir" in [argument.arg for argument in node.args.args]
        and "pytest.mark.realdata" not in [ast.unparse(decorator) for decorator in node.decorator_list]
    ]
    assert unmarked == []


_MODEL = "def test_model():\n    assert True\n"
_DATA = (
    "import pytest\n\n\n@pytest.mark.realdata\ndef test_data():\n    assert True\n"
)
_UI = "def test_ui():\n    assert True\n"
_CRASH = (
    "import os\nimport signal\n\nimport pytest\n\n\n@pytest.mark.realdata\ndef test_data_crashes():\n"
    "    os.kill(os.getpid(), signal.SIGSEGV)\n"
)
_HANG = "import time\n\n\ndef test_ui_hangs():\n    time.sleep(200)\n"
_PATCHES = """\
import builtins
import json
import os

import pytest


def _refuse(*_args, **_kwargs):
    raise OSError(28, "No space left on device")


@pytest.mark.realdata
def test_data_patches_what_a_writer_uses(monkeypatch):
    for module, name in ((json, "dumps"), (builtins, "open"), (os, "open"), (os, "replace")):
        monkeypatch.setattr(module, name, _refuse)
"""


_ABORT = "import os\n\n\ndef test_model_aborts():\n    os.abort()\n"


def _scratch_suite(tmp_path, data, ui, model=_MODEL):
    """A tree laid out as the repository: pyproject.toml at its root, conftest.py at launcher/tests/, one module per
    group (the model module under one of the model group's names)."""
    tests = tmp_path / "launcher" / "tests"
    tests.mkdir(parents=True)
    shutil.copy(ROOT / "pyproject.toml", tmp_path / "pyproject.toml")
    shutil.copy(HERE / "conftest.py", tests / "conftest.py")
    (tmp_path / "launcher" / "__init__.py").write_text("")
    (tests / "__init__.py").write_text("")
    (tests / "test_settings_document.py").write_text(model)
    (tests / "test_data_probe.py").write_text(data)
    (tests / "test_ui_probe.py").write_text(ui)
    return tmp_path


def _group_line(group, outcome, output):
    """The group table has the line for `group`, one test, with `outcome` (and what it names)."""
    return re.search(rf"^{group}\s+1 tests\s+{re.escape(outcome)}$", output, re.MULTILINE) is not None


def _run(root, *arguments, environment=None):
    return subprocess.run(
        [sys.executable, "-m", "pytest", "-vv", "-p", "no:cacheprovider", *arguments, "launcher/tests/"],
        cwd=root, env=environment or _environment(), capture_output=True, text=True, timeout=300, check=False,
    )


@pytest.mark.timeout(300)
def test_a_crash_fails_its_group_and_the_others_report(tmp_path):
    """T3: a segmentation fault in the data group fails the session, names the test and the signal, and the model and
    ui groups still run and pass."""
    result = _run(_scratch_suite(tmp_path, _CRASH, _UI), "--timeout=60", "--timeout-method=thread")
    output = result.stdout + result.stderr
    assert result.returncode == 1, output[-3000:]
    assert "Fatal Python error: Segmentation fault" in output
    assert "test_data_crashes" in output[output.index("Fatal Python error"):]
    assert "test_settings_document.py::test_model PASSED" in output
    assert "test_ui_probe.py::test_ui PASSED" in output
    assert _group_line("data", "died: SIGSEGV (in launcher/tests/test_data_probe.py::test_data_crashes)", output)
    assert _group_line("model", "passed", output) and _group_line("ui", "passed", output)
    assert "FAILED launcher/tests/test_data_probe.py::test_data_crashes" in output
    assert re.search(r"=+ 1 failed, 2 passed in ", output)


@pytest.mark.timeout(300)
def test_a_hang_fails_its_group_and_the_others_report(tmp_path):
    """T4: a test that hangs in the ui group is stopped by the timeout with the stacks dumped, the session fails, and
    the model and data groups still run and pass."""
    result = _run(_scratch_suite(tmp_path, _DATA, _HANG), "--timeout=5", "--timeout-method=thread")
    output = result.stdout + result.stderr
    assert result.returncode == 1, output[-3000:]
    assert "Timeout" in output and "test_ui_hangs" in output
    assert "test_settings_document.py::test_model PASSED" in output
    assert "test_data_probe.py::test_data PASSED" in output
    assert _group_line("ui", "failed: exit status 1 (in launcher/tests/test_ui_probe.py::test_ui_hangs)", output)
    assert _group_line("model", "passed", output) and _group_line("data", "passed", output)
    assert "FAILED launcher/tests/test_ui_probe.py::test_ui_hangs" in output
    assert re.search(r"=+ 1 failed, 2 passed in ", output)


@pytest.mark.timeout(300)
def test_a_test_that_patches_json_open_and_os_is_still_recorded(tmp_path):
    """The record is written while the test's own monkeypatches are in force: a test that patches json.dumps, open,
    os.open and os.replace (as test_settings_document.py's write-failure tests do) is recorded like any other."""
    result = _run(_scratch_suite(tmp_path, _PATCHES, _UI), "--timeout=60", "--timeout-method=thread")
    output = result.stdout + result.stderr
    assert result.returncode == 0, output[-3000:]
    assert "INTERNALERROR" not in output
    assert re.search(r"=+ 3 passed in ", output)


@pytest.mark.timeout(300)
def test_an_abort_in_the_first_group_still_runs_the_other_two(tmp_path):
    """T7: the groups run in order, and all of them run: an os.abort() in the first group (model) fails that group with
    SIGABRT, named, and the data and ui groups still run and pass."""
    result = _run(_scratch_suite(tmp_path, _DATA, _UI, model=_ABORT), "--timeout=60", "--timeout-method=thread")
    output = result.stdout + result.stderr
    assert result.returncode == 1, output[-3000:]
    assert output.index("launcher group model") < output.index("launcher group data") < output.index("launcher group ui")
    assert _group_line("model", "died: SIGABRT (in launcher/tests/test_settings_document.py::test_model_aborts)", output)
    assert _group_line("data", "passed", output) and _group_line("ui", "passed", output)
    assert re.search(r"=+ 1 failed, 2 passed in ", output)


@pytest.mark.timeout(300)
def test_a_session_given_a_group_runs_it_itself(tmp_path):
    """T7: with the group variable set (a group's own process), the session spawns nothing: it keeps that group's
    tests, deselects the others and runs in-process, printing no group table."""
    result = _run(_scratch_suite(tmp_path, _DATA, _UI), environment=_environment(**{GROUP_VARIABLE: "ui"}))
    output = result.stdout + result.stderr
    assert result.returncode == 0, output[-3000:]
    assert "launcher group" not in output and "launcher groups" not in output
    assert "test_ui_probe.py::test_ui PASSED" in output
    assert re.search(r"=+ 1 passed, 2 deselected in ", output)


@pytest.mark.timeout(300)
def test_one_group_runs_in_the_session_itself(tmp_path):
    """With one group only (here, a session on the model module alone), there is nothing to isolate: the session runs
    its tests itself and prints no group table."""
    root = _scratch_suite(tmp_path, _DATA, _UI)
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-vv", "-p", "no:cacheprovider", "launcher/tests/test_settings_document.py"],
        cwd=root, env=_environment(), capture_output=True, text=True, timeout=300, check=False,
    )
    assert result.returncode == 0, result.stdout[-3000:]
    assert "1 passed" in result.stdout
    assert "launcher groups" not in result.stdout


@pytest.mark.timeout(300)
def test_every_group_passing_reads_as_one_session(tmp_path):
    """When every group passes, the session's own summary counts every group's tests and its exit status is 0."""
    result = _run(_scratch_suite(tmp_path, _DATA, _UI), "--timeout=60", "--timeout-method=thread")
    output = result.stdout + result.stderr
    assert result.returncode == 0, output[-3000:]
    assert re.search(r"=+ 3 passed in ", output)
    assert all(_group_line(group, "passed", output) for group in ("model", "data", "ui"))
