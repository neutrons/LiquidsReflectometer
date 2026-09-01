# launcher/tests/test_harness.py
import os
import shutil
import signal
import subprocess
import sys
from pathlib import Path

import pytest
from qtpy import QtCore, QtWidgets


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
    """The roundtrip test passes even with isolation broken; this one does not."""
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
# silently regressed.
@pytest.mark.timeout(30)
def test_timeout_backstop_fires(tmp_path):
    """Self-test for the shipped configuration: if pytest-timeout goes missing
    or timeout_method regresses to signal, this goes red instead of a future
    slug hanging for hours.

    Two mechanics matter as much as the assertions:

    * The inner run reads the repo's real pyproject.toml, copied in, rather
      than a hand-written pytest.ini. An ini of our own making would detach the
      test from the shipped config — deleting `timeout_method = "thread"` from
      pyproject would leave this green, which is exactly the regression it
      exists to catch.
    * Output goes to files and the process is killed in a finally. With pipes,
      a SIGKILL of the outer pytest leaves the inner hang alive forever:
      pytest-timeout's own `finally` raises BrokenPipeError on the dead pipe,
      which pre-empts its os._exit. One such orphan ran 52 minutes at 99.9% CPU
      on this host before it was found and killed.
    """
    shutil.copy(Path(__file__).parent / "conftest.py", tmp_path / "conftest.py")
    shutil.copy(Path(__file__).resolve().parents[2] / "pyproject.toml", tmp_path / "pyproject.toml")
    (tmp_path / "test_hangs.py").write_text(
        "import pytest\n"
        "from qtpy import QtWidgets\n"
        "\n"
        "@pytest.mark.timeout(3)\n"
        "def test_hangs(isolated_qapp):\n"
        "    QtWidgets.QMessageBox.warning(None, 't', 'blocks in the C++ event loop')\n"
    )
    out_path = tmp_path / "inner-stdout.txt"
    err_path = tmp_path / "inner-stderr.txt"

    with open(out_path, "wb") as out, open(err_path, "wb") as err:
        proc = subprocess.Popen(
            [sys.executable, "-m", "pytest", "-p", "no:cacheprovider", str(tmp_path / "test_hangs.py")],
            cwd=str(tmp_path),
            stdout=out,
            stderr=err,
            start_new_session=True,
        )
        try:
            returncode = proc.wait(timeout=25)
        finally:
            if proc.poll() is None:
                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                proc.wait(timeout=10)

    combined = out_path.read_bytes() + err_path.read_bytes()
    tail = combined.decode(errors="replace")[-1500:]
    # No header assertion: pytest-timeout prints "timeout: …/method: …" only for
    # a session-level timeout, and the budget here comes from a marker. The kill
    # itself is the method evidence — the hang is inside Qt's C++ event loop,
    # where the signal method provably never fires (measured: signal ran past a
    # 45 s external kill; thread ends it on schedule).
    assert returncode != 0, tail
    assert b"Timeout" in combined, tail
    assert b"test_hangs" in combined, tail


def test_native_format_also_lands_in_tmp(isolated_qapp, tmp_path):
    """The redirect covers both formats. A NativeFormat store constructed
    explicitly would escape a single-format redirect — which is the shape of the
    bug v1 shipped."""
    settings = QtCore.QSettings(QtCore.QSettings.NativeFormat, QtCore.QSettings.UserScope, "test-org", "test-app")
    settings.setValue("harness/probe", "x")
    settings.sync()
    assert settings.fileName().startswith(str(tmp_path))


def test_import_scope_redirect_holds_without_any_fixture(tmp_path):
    """Pins the import block: a QSettings built at import time, before any
    fixture runs, must not touch the real config.

    Runs in a subprocess with a throwaway HOME because the property is about
    what happens *before* a fixture can intervene — an in-process test has
    already imported the conftest and cannot observe it.
    """
    home = tmp_path / "home"
    (home / ".config").mkdir(parents=True)
    probe = tmp_path / "probe.py"
    probe.write_text(
        "import sys\n"
        "sys.path.insert(0, %r)\n" % str(Path(__file__).parent)
        + "import conftest  # installs the redirect at import\n"
        "from qtpy import QtCore\n"
        "QtCore.QCoreApplication.setOrganizationName('probe-org')\n"
        "QtCore.QCoreApplication.setApplicationName('probe-app')\n"
        "paths = [QtCore.QSettings().fileName()]\n"
        "paths.append(QtCore.QSettings(QtCore.QSettings.NativeFormat, QtCore.QSettings.UserScope,"
        " 'probe-org', 'probe-app').fileName())\n"
        "print('PATHS' + repr(paths))\n"
    )
    env = dict(os.environ, HOME=str(home), XDG_CONFIG_HOME=str(home / ".config"), QT_QPA_PLATFORM="offscreen")
    proc = subprocess.run(
        [sys.executable, str(probe)],
        cwd=str(Path(__file__).resolve().parents[2]),
        env=env,
        capture_output=True,
        timeout=120,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr.decode(errors="replace")[-2000:]
    paths = eval(proc.stdout.decode().split("PATHS", 1)[1].strip())  # noqa: S307 — our own literal
    for path in paths:
        assert not path.startswith(str(home)), f"settings escaped to the real config root: {path}"
        assert "launcher-tests-settings-" in path, f"not inside the scratch root: {path}"
