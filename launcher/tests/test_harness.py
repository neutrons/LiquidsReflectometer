# launcher/tests/test_harness.py
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from qtpy import QtCore, QtWidgets

# Constructed at import, before any fixture can run: this is what caches Qt's
# settings root for the process. If the conftest's redirect were fixture-scoped
# (or format-partial), test_qsettings_file_lands_in_tmp below would fail
# regardless of test order — which is the point of poisoning it here rather
# than relying on which test happens to run first.
_ROOT_POISON = QtCore.QSettings()


def test_no_qmessagebox_fixture_neutralizes_modals(isolated_qapp, no_qmessagebox):
    # Without the fixture this call blocks forever under offscreen Qt
    # (the 10-hour-orphan class); with it, it returns immediately.
    assert QtWidgets.QMessageBox.warning(None, "t", "must not block") is None


def test_isolated_qsettings_roundtrip(isolated_qapp):
    settings = QtCore.QSettings()
    settings.setValue("harness/probe", "x")
    settings.sync()
    assert settings.value("harness/probe") == "x"
    assert "test-org" in QtCore.QCoreApplication.organizationName()


def test_qsettings_file_lands_in_tmp(isolated_qapp, tmp_path):
    """The roundtrip test passes even with isolation broken; this one does not.

    Order-independent: the module-level QSettings above has already cached the
    root by the time this runs.
    """
    settings = QtCore.QSettings()
    settings.setValue("harness/probe", "x")
    settings.sync()
    assert QtCore.QSettings().fileName().startswith(str(tmp_path))


def test_dialog_exec_neutralized(isolated_qapp, no_qmessagebox):
    """The surface that actually blocks in launcher/: a plain QDialog.

    exec_ is defined on QDialog and only inherited by QMessageBox, so a patch
    aimed at the subclass leaves every production site live.
    """
    dialog = QtWidgets.QDialog()
    assert dialog.exec_() == QtWidgets.QDialog.Accepted


def test_dialog_exec_returns_accepted_not_ok(isolated_qapp, no_qmessagebox):
    """Guards the value, not just the non-blocking.

    Call sites compare against 1 / QDialog.Accepted. QMessageBox.Ok is 1024, so
    returning it from the base-class patch would silently convert every accept
    into a cancel.
    """
    assert QtWidgets.QDialog().exec_() == 1
    assert QtWidgets.QDialog().exec_() != QtWidgets.QMessageBox.Ok


def test_file_dialogs_neutralized_by_default(isolated_qapp):
    """no_qfiledialog is autouse — no test opts in, and all 26 sites block."""
    assert QtWidgets.QFileDialog.getOpenFileName(None, "t") == ("", "")
    assert QtWidgets.QFileDialog.getSaveFileName(None, "t") == ("", "")
    assert QtWidgets.QFileDialog.getExistingDirectory(None, "t") == ""
    first = QtWidgets.QFileDialog.getOpenFileNames(None, "t")
    assert first == ([], "")
    first[0].append("mutated")
    assert QtWidgets.QFileDialog.getOpenFileNames(None, "t") == ([], "")


def test_collection_hook_arms_only_launcher_items():
    """Scope and precedence of the timeout hook, checked directly.

    A subprocess cannot show this cheaply, and the property that matters most
    is a negative one: a combined `pytest tests launcher/tests` invocation must
    not time-limit the reduction suite.
    """
    from launcher.tests.conftest import pytest_collection_modifyitems

    class _Config:
        def __init__(self, cli_timeout=None):
            self._cli_timeout = cli_timeout

        def getoption(self, name, default=None):
            return self._cli_timeout if name == "timeout" else default

    class _Item:
        def __init__(self, path, config, marker=None):
            self.path = path
            self.config = config
            self._marker = marker
            self.added = []

        def get_closest_marker(self, name):
            return self._marker if name == "timeout" else None

        def add_marker(self, marker):
            self.added.append(marker)

    here = Path(__file__).parent
    reduction = here.parent.parent / "tests" / "test_something.py"
    config = _Config()

    launcher_item = _Item(here / "test_x.py", config)
    reduction_item = _Item(reduction, config)
    marked_item = _Item(here / "test_y.py", config, marker=object())
    cli_item = _Item(here / "test_z.py", _Config(cli_timeout=5))

    pytest_collection_modifyitems([launcher_item, reduction_item, marked_item, cli_item])

    assert len(launcher_item.added) == 1
    assert launcher_item.added[0].args == (120,)
    assert reduction_item.added == [], "the reduction suite must not be time-limited"
    assert marked_item.added == [], "an explicit marker must win"
    assert cli_item.added == [], "an explicit --timeout must win"


# The inner run hangs inside Qt's C++ event loop, which is the case that
# discriminates the timeout method: a pure-Python `while True` loop is killed by
# the signal method too, so a busy-loop self-test would pass even if the method
# silently regressed. Budgeted at 10 s via a marker in the generated file, so
# this costs seconds rather than the conftest's 120 s default.
@pytest.mark.timeout(120)
def test_timeout_backstop_fires(tmp_path):
    """Self-test for the shipped configuration: if pytest-timeout goes missing
    or timeout_method regresses to signal, this goes red instead of a future
    slug hanging for hours."""
    shutil.copy(Path(__file__).parent / "conftest.py", tmp_path / "conftest.py")
    (tmp_path / "pytest.ini").write_text("[pytest]\ntimeout = 10\ntimeout_method = thread\n")
    (tmp_path / "test_hangs.py").write_text(
        "import pytest\n"
        "from qtpy import QtWidgets\n"
        "\n"
        "@pytest.mark.timeout(10)\n"
        "def test_hangs(isolated_qapp):\n"
        "    QtWidgets.QMessageBox.warning(None, 't', 'blocks in the C++ event loop')\n"
    )

    proc = subprocess.run(
        [sys.executable, "-m", "pytest", "-p", "no:cacheprovider", str(tmp_path / "test_hangs.py")],
        cwd=str(tmp_path),
        capture_output=True,
        timeout=100,
        check=False,
    )
    combined = proc.stdout + proc.stderr
    tail = combined.decode(errors="replace")[-1500:]
    assert proc.returncode != 0, tail
    assert b"timeout method: thread" in combined, tail
    assert b"Timeout" in combined, tail
    assert b"QMessageBox" in combined or b"test_hangs" in combined, tail
