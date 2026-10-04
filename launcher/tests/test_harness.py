# launcher/tests/test_harness.py
import concurrent.futures
import os
import shutil
import signal
import subprocess
import sys
from pathlib import Path

import pytest
from qtpy import QtCore, QtWidgets, sip


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


def _run_inner_pytest(directory, module_name, source, wait):
    """Run pytest in a subprocess on one generated test module, under the shipped conftest and pyproject.

    Returns (exit status, stdout + stderr). The mechanics are the ones test_timeout_backstop_fires explains:
    the real conftest.py and pyproject.toml are copied in, output goes to files, and the process group is
    killed in a finally. ``--basetemp`` inside `directory` keeps runs made side by side from cleaning up one
    another's temporary trees.
    """
    directory.mkdir(parents=True, exist_ok=True)
    shutil.copy(Path(__file__).parent / "conftest.py", directory / "conftest.py")
    shutil.copy(Path(__file__).resolve().parents[2] / "pyproject.toml", directory / "pyproject.toml")
    (directory / module_name).write_text(source)
    out_path = directory / "inner-stdout.txt"
    err_path = directory / "inner-stderr.txt"

    with open(out_path, "wb") as out, open(err_path, "wb") as err:
        proc = subprocess.Popen(
            [sys.executable, "-m", "pytest", "-p", "no:cacheprovider", f"--basetemp={directory / 'basetemp'}",
             str(directory / module_name)],
            cwd=str(directory),
            stdout=out,
            stderr=err,
            start_new_session=True,
        )
        try:
            returncode = proc.wait(timeout=wait)
        finally:
            if proc.poll() is None:
                os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
                proc.wait(timeout=10)
    return returncode, out_path.read_bytes() + err_path.read_bytes()


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
    returncode, combined = _run_inner_pytest(
        tmp_path,
        "test_hangs.py",
        "import pytest\n"
        "from qtpy import QtWidgets\n"
        "\n"
        "@pytest.mark.timeout(3)\n"
        "def test_hangs(isolated_qapp):\n"
        "    QtWidgets.QMessageBox.warning(None, 't', 'blocks in the C++ event loop')\n",
        wait=25,
    )
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


# --------------------------------------------------------------------------
# launcher-test-teardown: the teardown frees only the windows nothing else owns
# --------------------------------------------------------------------------

# The old teardown crashed in about half the runs: 54 of 96 here, about 0.48 for the test reviewer, and 4 of 12 in
# its lowest invocation. At those rates, 20 passing runs under it would be an event of 7e-8 to 3e-4.
_TEARDOWN_RUNS = 20

_FREED_TYPES = ("Window", "Dialog", "Tool", "Sheet", "Drawer")
_OWNED_TYPES = ("Popup", "ToolTip", "SplashScreen", "SubWindow", "ForeignWindow", "CoverWindow")

_COMPLETER_POPUP_MODULE = """\
from qtpy import QtWidgets

# The test returns with its window open, as a test with a reference cycle does: the teardown frees it.
_LEFT_OPEN = []


def test_opens_a_completer_popup(isolated_qapp):
    window = QtWidgets.QWidget()
    edit = QtWidgets.QLineEdit(window)
    completer = QtWidgets.QCompleter(["alpha", "alpine", "beta"], edit)
    edit.setCompleter(completer)
    window.show()
    completer.setCompletionPrefix("al")
    completer.complete()
    assert completer.popup().isVisible()
    _LEFT_OPEN.append(window)
"""

_MENU_POPUP_MODULE = """\
from qtpy import QtCore, QtWidgets

_LEFT_OPEN = []


def test_opens_a_menu(isolated_qapp):
    window = QtWidgets.QWidget()
    menu = QtWidgets.QMenu(window)
    menu.addAction("first")
    window.show()
    menu.popup(window.mapToGlobal(QtCore.QPoint(10, 10)))
    assert menu.isVisible()
    _LEFT_OPEN.append(window)
"""

# Both modules keep the QApplication across their two tests: the shared application the drain exists for. With
# a fresh one per test, the windows the old one outlived are freed when it goes, which would hide a teardown that
# drains nothing (measured: no drain, a fresh application, both freed; no drain, the application kept, neither).
_LEFT_OPEN_WINDOWS_MODULE = """\
from qtpy import QtCore, QtWidgets, sip

_KEPT = []
_CLOSED = []


class _Window(QtWidgets.QWidget):
    def closeEvent(self, event):
        # The organization name when the drain closes the window: the fixture's own is test-org-<directory>.
        _CLOSED.append((type(self).__name__, QtCore.QCoreApplication.organizationName()))
        super().closeEvent(event)


def test_leaves_a_window_and_a_dialog_open(isolated_qapp):
    window = _Window()
    dialog = QtWidgets.QDialog()
    window.show()
    dialog.show()
    _KEPT.extend([isolated_qapp, window, dialog])


def test_the_teardown_closed_and_freed_them(isolated_qapp):
    app, window, dialog = _KEPT
    assert isolated_qapp is app
    assert [name for name, _ in _CLOSED] == ["_Window"]
    assert not _CLOSED[0][1].startswith("test-org-"), _CLOSED  # the identity was restored before the drain
    assert [sip.isdeleted(window), sip.isdeleted(dialog)] == [True, True]


def test_no_fixture_identity_is_left_after_the_fixture():
    org = QtCore.QCoreApplication.organizationName()
    domain = QtCore.QCoreApplication.organizationDomain()
    app = QtCore.QCoreApplication.applicationName()
    assert not org.startswith("test-org-") and domain != "example.test" and not app.startswith("test-app-")
"""

_FREED_DURING_DRAIN_MODULE = """\
from qtpy import QtWidgets, sip

_KEPT = []


class _Partner(QtWidgets.QWidget):
    # Closing either window frees the other at once, so the one the drain reaches second is already destroyed.
    partner = None

    def closeEvent(self, event):
        if self.partner is not None and not sip.isdeleted(self.partner):
            sip.delete(self.partner)
        super().closeEvent(event)


def test_leaves_two_windows_that_free_each_other(isolated_qapp):
    first, second = _Partner(), _Partner()
    first.partner, second.partner = second, first
    first.show()
    second.show()
    _KEPT.extend([isolated_qapp, first, second])


def test_the_teardown_got_past_the_destroyed_one(isolated_qapp):
    app, first, second = _KEPT
    assert isolated_qapp is app
    assert [sip.isdeleted(first), sip.isdeleted(second)] == [True, True]
"""


def _repeat_inner_pytest(tmp_path, module_name, source, runs, wait):
    """`_run_inner_pytest` `runs` times, each in a directory and a process of its own, as many at once as there
    are CPUs. Whether a run crashes is decided inside its process, so the count is over processes (plan H3)."""
    with concurrent.futures.ThreadPoolExecutor(max_workers=min(runs, os.cpu_count() or 1)) as pool:
        futures = [
            pool.submit(_run_inner_pytest, tmp_path / f"run-{index:02d}", module_name, source, wait)
            for index in range(runs)
        ]
        return [future.result() for future in futures]


def _excerpt(output):
    """Where a crashed run says why (faulthandler's header and the frames under it), else the output's tail."""
    text = output.decode(errors="replace")
    start = text.find("Fatal Python error")
    return text[start : start + 1500] if start >= 0 else text[-1500:]


def _assert_every_run_passed(results, tests):
    failed = [(status, output) for status, output in results if status != 0]
    detail = _excerpt(failed[0][1]) if failed else ""
    assert not failed, f"{len(failed)} of {len(results)} runs failed, exit statuses {[s for s, _ in failed]}:\n{detail}"
    for _, output in results:
        assert f"{tests} passed".encode() in output, _excerpt(output)


@pytest.mark.timeout(300)
def test_teardown_survives_an_open_completer_popup(tmp_path):
    """H3: a test that leaves a QCompleter's pop-up open survives the teardown, in every run of _TEARDOWN_RUNS.

    The pop-up is a parentless top-level window that the completer owns through a raw pointer and deletes in
    its own destructor. A teardown that also deletes it frees it twice. gdb put the crash in
    QCompleter::~QCompleter, under the fixture's DeferredDelete flush. A run crashes when Qt lists the pop-up
    before its window, which it does in about half the runs, so each run is a process of its own.
    """
    results = _repeat_inner_pytest(tmp_path, "test_completer_popup.py", _COMPLETER_POPUP_MODULE, _TEARDOWN_RUNS, 120)
    _assert_every_run_passed(results, tests=1)


@pytest.mark.timeout(300)
def test_teardown_survives_an_open_menu(tmp_path):
    """H3: the same for a QMenu its window owns, opened with popup()."""
    results = _repeat_inner_pytest(tmp_path, "test_menu_popup.py", _MENU_POPUP_MODULE, _TEARDOWN_RUNS, 120)
    _assert_every_run_passed(results, tests=1)


@pytest.mark.timeout(120)
def test_teardown_frees_a_window_and_a_dialog_left_open(tmp_path):
    """H2, through the fixture: a plain window and a dialog left open by one test are closed (the window's
    closeEvent runs) and freed before the next test, on the same QApplication. The identity is restored first:
    when the drain closes the window, the organization name is no longer the fixture's, and after the fixture
    none of the three names it set is left."""
    status, output = _run_inner_pytest(tmp_path, "test_left_open.py", _LEFT_OPEN_WINDOWS_MODULE, 100)
    assert status == 0 and b"3 passed" in output, _excerpt(output)


@pytest.mark.timeout(120)
def test_teardown_gets_past_a_window_already_freed_during_the_drain(tmp_path):
    """H2: a window freed while the drain closes another one is already destroyed when the drain reaches it.
    The RuntimeError that raises is swallowed and the drain goes on, so the teardown does not error."""
    status, output = _run_inner_pytest(tmp_path, "test_freed_during_drain.py", _FREED_DURING_DRAIN_MODULE, 100)
    assert status == 0 and b"2 passed" in output, _excerpt(output)


@pytest.mark.parametrize("name", _FREED_TYPES + _OWNED_TYPES)
def test_the_teardown_frees_a_window_only_when_nothing_owns_it(isolated_qapp, name):
    """H1, the rule as a table. A parentless top-level widget of a type a test makes as a window is the
    teardown's to free. A pop-up, a tooltip, a splash screen and the rarer kinds belong to whoever opened them,
    and a window with a parent belongs to its parent. The trade: a parentless pop-up that nothing owns is left
    to die with the process, a leak rather than a crash. (A parentless SubWindow becomes type 0x13, measured.)
    """
    from launcher.tests.conftest import _is_test_window

    flag = getattr(QtCore.Qt, name)
    owner = QtWidgets.QWidget()
    alone, parented = QtWidgets.QWidget(None, flag), QtWidgets.QWidget(owner, flag)
    assert alone.isWindow()
    assert _is_test_window(alone) is (name in _FREED_TYPES)
    assert _is_test_window(parented) is False


def test_the_drain_frees_the_test_windows_now_and_leaves_the_owned_ones(isolated_qapp):
    """H1 and H2 through the drain itself. When it returns, every parentless test window is already freed,
    shown or hidden (by the DeferredDelete flush, not by a later event loop). Every owned top-level is
    untouched, shown or hidden, a dialog whose owner survives included."""
    from launcher.tests.conftest import _drain_test_windows

    freed = [QtWidgets.QWidget(None, getattr(QtCore.Qt, name)) for name in _FREED_TYPES]
    owner = QtWidgets.QWidget(None, QtCore.Qt.Popup)
    owned = [QtWidgets.QWidget(None, getattr(QtCore.Qt, name)) for name in _OWNED_TYPES]
    owned += [owner, QtWidgets.QDialog(owner)]
    for index, widget in enumerate(freed + owned):
        if index % 2 == 0:
            widget.show()
    assert any(not widget.isVisible() for widget in freed) and any(widget.isVisible() for widget in freed)
    _drain_test_windows(isolated_qapp)
    assert [sip.isdeleted(widget) for widget in freed] == [True] * len(freed)
    assert [sip.isdeleted(widget) for widget in owned] == [False] * len(owned)
