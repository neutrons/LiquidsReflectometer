import atexit
import json
import os
import shutil
import signal
import subprocess
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from qtpy import QtCore, QtWidgets

# --- Process-wide settings redirect, installed at IMPORT ---------------------
# Qt caches the settings root at the first QSettings construction in a process.
# A fixture cannot win that race: any module constructing QSettings at import or
# collection time would already have bound the developer's real ~/.config for
# the whole session. setPath is consulted per construction, so installing it
# here — before collection — is what actually isolates.
#
# Both formats are redirected: setDefaultFormat(IniFormat) governs new
# QSettings(), but anything constructing NativeFormat explicitly would escape a
# single-format redirect.
_SCRATCH_ROOT = tempfile.mkdtemp(prefix="launcher-tests-settings-")
# Only this process's own root — a concurrent run owns a different one.
atexit.register(shutil.rmtree, _SCRATCH_ROOT, True)
QtCore.QSettings.setDefaultFormat(QtCore.QSettings.IniFormat)
for _format in (QtCore.QSettings.IniFormat, QtCore.QSettings.NativeFormat):
    QtCore.QSettings.setPath(_format, QtCore.QSettings.UserScope, _SCRATCH_ROOT)

# --- matplotlib's cache, preset at IMPORT ------------------------------------
# Importing launcher.new_launcher chooses MPLCONFIGDIR (launcher/runtime_env.py)
# and leaves a value already set alone. Presetting it here keeps a test that
# imports the launcher from creating the per-user directory under /var/tmp, and
# keeps this run's font cache in the scratch root. setdefault: a value the
# developer or CI set wins, as it does for the launcher.
os.environ.setdefault("MPLCONFIGDIR", os.path.join(_SCRATCH_ROOT, "matplotlib"))


# --- What the teardown frees -------------------------------------------------
# Only the windows a test made and nothing else owns: a parentless top-level
# widget of one of these types. Everything else is freed by its owner: a pop-up
# (a QCompleter's list, a QMenu), a tooltip, a splash screen, the rarer types,
# and any window with a parent. A QCompleter holds its pop-up through a raw
# pointer and deletes it in its own destructor. When the teardown freed every
# top-level widget, the pop-up was freed twice whenever Qt listed it before its
# window. The run then crashed in QCompleter::~QCompleter (gdb) under the flush
# at the end of _drain_test_windows. That happened in 5, 8, 7, 7 and 8 of 12
# subprocess runs of
# test_harness.py::test_teardown_survives_an_open_completer_popup; under this
# rule, 0 of 12 in each of four runs. The cost: a parentless pop-up that nothing
# owns is left to die with the process, a leak rather than a crash. Types are
# compared exactly, since Tool is Popup|Dialog as bits.
_TEST_WINDOW_TYPES = frozenset(
    {QtCore.Qt.Window, QtCore.Qt.Dialog, QtCore.Qt.Tool, QtCore.Qt.Sheet, QtCore.Qt.Drawer}
)


def _is_test_window(widget):
    """True when the teardown frees `widget`: nothing owns it, and a test made it as a window."""
    return widget.parent() is None and widget.windowType() in _TEST_WINDOW_TYPES


def _drain_test_windows(app):
    """Close and free the windows the test left open, now rather than at some later event loop."""
    for widget in QtWidgets.QApplication.topLevelWidgets():
        try:
            if _is_test_window(widget):
                widget.close()
                widget.deleteLater()
        except RuntimeError:
            # Already destroyed on the C++ side; nothing to drain.
            pass
    app.processEvents()
    # deleteLater only queues; without this the objects are never freed and the
    # "drain" drains nothing.
    app.sendPostedEvents(None, QtCore.QEvent.DeferredDelete)


@pytest.fixture
def isolated_qapp(tmp_path, monkeypatch):
    """QApplication with a per-test QSettings root.

    The import-scope redirect above keeps every test inside a scratch tree;
    this narrows it to one directory per test so stores cannot bleed between
    them. Teardown restores the org/app/domain names but deliberately NOT
    defaultFormat: setPath has no undo, so restoring the format would re-open
    exactly the window this fixture exists to close.
    """
    prev_org = QtCore.QCoreApplication.organizationName()
    prev_domain = QtCore.QCoreApplication.organizationDomain()
    prev_app = QtCore.QCoreApplication.applicationName()

    for fmt in (QtCore.QSettings.IniFormat, QtCore.QSettings.NativeFormat):
        QtCore.QSettings.setPath(fmt, QtCore.QSettings.UserScope, str(tmp_path))
    QtCore.QCoreApplication.setOrganizationName(f"test-org-{tmp_path.name}")
    QtCore.QCoreApplication.setOrganizationDomain("example.test")
    QtCore.QCoreApplication.setApplicationName(f"test-app-{tmp_path.name}")
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))

    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    try:
        yield app
    finally:
        # Identity first: a widget that raises during the drain must not cost
        # the next test its isolation.
        #
        # On this environment the fixture drops its reference to the
        # QApplication, so the common case is a fresh one per test rather than
        # one shared across the session. The drain is therefore cheap
        # defense-in-depth for the session-scoped case, not the load-bearing
        # cleanup an earlier version of this comment claimed it was.
        QtCore.QCoreApplication.setOrganizationName(prev_org)
        QtCore.QCoreApplication.setOrganizationDomain(prev_domain)
        QtCore.QCoreApplication.setApplicationName(prev_app)
        _drain_test_windows(app)


@pytest.fixture
def no_qmessagebox(monkeypatch):
    """Neutralize modal dialogs that block forever under offscreen Qt.

    The patch targets QDialog, not QMessageBox: exec_/exec are defined on
    QDialog and merely inherited by QMessageBox, so patching the subclass
    shadows an attribute the production call sites never reach — every blocking
    site in launcher/ is a plain QDialog(...).exec_().

    The return value is QDialog.Accepted (1), never QMessageBox.Ok (1024):
    call sites compare against 1 / QDialog.Accepted, so returning Ok from a
    base-class patch would silently turn accept paths into cancel paths.
    """
    for name in ("warning", "information", "critical", "question"):
        monkeypatch.setattr(
            QtWidgets.QMessageBox,
            name,
            staticmethod(lambda *_args, **_kwargs: None),
        )
    for name in ("exec_", "exec", "open"):
        if hasattr(QtWidgets.QDialog, name):
            monkeypatch.setattr(
                QtWidgets.QDialog,
                name,
                lambda *_args, **_kwargs: QtWidgets.QDialog.Accepted,
            )


@pytest.fixture(autouse=True)
def no_qfiledialog(monkeypatch):
    """Neutralize modal file dialogs (26 call sites in launcher/).

    autouse because nothing opts in and every one of those sites blocks; a test
    that wants a specific chosen path patches the method itself, which overrides
    this.
    """
    for name, result in (
        ("getOpenFileName", ("", "")),
        ("getSaveFileName", ("", "")),
        ("getOpenFileNames", None),  # fresh list per call — see below
        ("getExistingDirectory", ""),
    ):
        if result is None:
            monkeypatch.setattr(
                QtWidgets.QFileDialog,
                name,
                staticmethod(lambda *_args, **_kwargs: ([], "")),
            )
        else:
            monkeypatch.setattr(
                QtWidgets.QFileDialog,
                name,
                staticmethod(lambda *_args, _result=result, **_kwargs: _result),
            )


def pytest_collection_modifyitems(items):
    """Arm a timeout for the tests this conftest governs, and, in a group's process, keep only that group's tests.

    Scoped two ways: only items under this directory (a combined invocation
    such as `pytest tests launcher/tests` must not time-limit the reduction
    suite), and only when no explicit --timeout was passed (the documented RED
    technique of running at a short timeout has to win).
    """
    here = Path(__file__).parent
    for item in items:
        if item.config.getoption("timeout", None) is not None:
            continue
        try:
            governed = Path(item.path).is_relative_to(here)
        except (AttributeError, ValueError):
            continue
        if governed and item.get_closest_marker("timeout") is None:
            item.add_marker(pytest.mark.timeout(120))
    group = os.environ.get(_GROUP_VARIABLE)
    config = items[0].config if items else None
    if group and config is not None and _in_the_repository(config):
        deselected = [item for item in items if _group_of(item) != group]
        if deselected:
            config.hook.pytest_deselected(items=deselected)
            items[:] = [item for item in items if _group_of(item) == group]


# --- One pytest process per process-risk group --------------------------------
# A native crash (a segmentation fault in Qt or in Mantid) or a hang ends the
# process it happens in. With the whole launcher suite in one process, one
# crash ended every test after it, and the record named neither the test nor
# the signal. So the repository's own launcher session runs each group in a
# pytest process of its own, always all of them:
#   model  the Qt-free model tests: no QApplication;
#   data   the real-data tests, marked `realdata` (Mantid's LoadEventNexus on
#          the fixture runs);
#   ui     everything else, the Qt session; a new test module lands here.
# A group's process is the same pytest invocation with the group named in the
# environment: it keeps that group's tests, deselects the rest, and writes each
# test's start and reports to a file. The session folds those reports into its
# own, so its summary, its counts and its exit status read as if every test
# had run in it; a test that started and never reported is the one its
# process died in, and is reported failed with the exit status or signal. In
# the group's process, faulthandler prints the stack of a native crash and
# pytest-timeout's thread method the stacks of a hang. Nothing changes when the
# session holds one group only, under --collect-only or --pdb, or outside the
# repository: the harness's inner runs copy this file to the root of a scratch
# directory, where it does nothing.
_GROUP_VARIABLE = "LR_LAUNCHER_TEST_GROUP"
_RECORD_VARIABLE = "LR_LAUNCHER_TEST_RECORD"
_GROUPS = ("model", "data", "ui")
_MODEL_MODULES = frozenset({"test_settings_document.py", "test_roi_estimate.py"})
_GROUP_RESULTS = []
_RECORDING = []  # the record file's descriptor, in a group's process of the repository's session
# The record is written while other tests' monkeypatches are in force (json.dumps, open, os.replace are all patched
# somewhere in this suite), so it uses only what it took at import.
_ENCODE = json.JSONEncoder().encode
_OS_OPEN, _OS_WRITE = os.open, os.write


def _in_the_repository(config):
    """True in the repository's launcher session: this file is the root's launcher/tests/conftest.py."""
    rootpath = getattr(config, "rootpath", None)
    if rootpath is None:
        return False
    return Path(__file__).resolve() == (Path(rootpath) / "launcher" / "tests" / "conftest.py").resolve()


def _group_of(item):
    if item.get_closest_marker("realdata") is not None:
        return "data"
    return "model" if Path(item.path).name in _MODEL_MODULES else "ui"


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "realdata: reads the real REF_L event files through Mantid; its own process in the launcher suite"
    )
    # The harness's inner runs inherit both variables from the ui group's process; outside the repository they
    # must not write into its record.
    if os.environ.get(_GROUP_VARIABLE) and os.environ.get(_RECORD_VARIABLE) and _in_the_repository(config):
        _RECORDING.append(_OS_OPEN(os.environ[_RECORD_VARIABLE], os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600))


def _record(entry):
    """In a group's process: one line per event, flushed, so that a crash loses nothing written before it."""
    if _RECORDING:
        _OS_WRITE(_RECORDING[0], (_ENCODE(entry) + "\n").encode("utf-8"))


def pytest_runtest_logstart(nodeid, location):
    _record({"start": nodeid, "location": list(location)})


def pytest_runtest_logreport(report):
    if _RECORDING:
        longrepr = report.longrepr
        if longrepr is not None and not isinstance(longrepr, tuple):
            longrepr = str(longrepr)
        _record({"nodeid": report.nodeid, "location": list(report.location), "when": report.when,
                 "outcome": report.outcome, "longrepr": longrepr, "wasxfail": getattr(report, "wasxfail", None),
                 "duration": report.duration})


def _replay(config, terminal, session, group, record_path, status):
    """Fold a group's reports into this session's; returns the test the process died in, if it died in one."""
    started, finished = [], set()
    if record_path.exists():
        for line in record_path.read_text(encoding="utf-8").splitlines():
            entry = json.loads(line)
            if "start" in entry:
                started.append((entry["start"], tuple(entry["location"])))
                continue
            longrepr = entry["longrepr"]
            report = pytest.TestReport(
                nodeid=entry["nodeid"], location=tuple(entry["location"]), keywords={}, outcome=entry["outcome"],
                longrepr=tuple(longrepr) if isinstance(longrepr, list) else longrepr, when=entry["when"],
                duration=entry["duration"],
            )
            if entry["wasxfail"] is not None:
                report.wasxfail = entry["wasxfail"]
            _count(config, terminal, session, report)
            if entry["when"] == "teardown":
                finished.add(entry["nodeid"])
    died_in = next(((nodeid, location) for nodeid, location in reversed(started) if nodeid not in finished), None)
    if status != 0 and died_in is not None:
        nodeid, location = died_in
        _count(config, terminal, session, pytest.TestReport(
            nodeid=nodeid, location=location, keywords={}, outcome="failed", when="call",
            longrepr=f"the {group} group's process {_describe(status)} in this test; its output above has the stack",
        ))
    return died_in[0] if died_in else None


def _count(config, terminal, session, report):
    category = config.hook.pytest_report_teststatus(report=report, config=config)[0]
    terminal.stats.setdefault(category, []).append(report)
    if report.failed:
        session.testsfailed += 1


def pytest_runtestloop(session):
    """Run each group present in a pytest process of its own, and report every group's tests as this session's."""
    config = session.config
    if os.environ.get(_GROUP_VARIABLE) or not _in_the_repository(config):
        return None
    if config.option.collectonly or config.option.usepdb:
        return None
    counts = {group: sum(_group_of(item) == group for item in session.items) for group in _GROUPS}
    groups = [group for group in _GROUPS if counts[group]]
    if len(groups) < 2:
        return None
    terminal = config.pluginmanager.getplugin("terminalreporter")
    capture = config.pluginmanager.getplugin("capturemanager")
    command = [sys.executable, "-m", "pytest", "-p", "faulthandler", *config.invocation_params.args]
    records = Path(tempfile.mkdtemp(prefix="launcher-groups-"))
    try:
        for group in groups:
            record_path = records / f"{group}.jsonl"
            environment = dict(os.environ, PYTHONFAULTHANDLER="1", **{_GROUP_VARIABLE: group,
                                                                    _RECORD_VARIABLE: str(record_path)})
            with capture.global_and_fixture_disabled():
                terminal.write_sep("=", f"launcher group {group}: {counts[group]} tests, in a process of its own")
                terminal.flush()
                status = subprocess.run(command, cwd=config.invocation_params.dir, env=environment,
                                        check=False).returncode
            died_in = _replay(config, terminal, session, group, record_path, status)
            if status != 0 and died_in is None:
                session.testsfailed += 1  # the process failed outside any test: collection, a session fixture
            _GROUP_RESULTS.append((group, counts[group], status, died_in))
    finally:
        shutil.rmtree(records, ignore_errors=True)
    return True


def _describe(status):
    if status == 0:
        return "passed"
    if status < 0:
        try:
            return f"died: {signal.Signals(-status).name}"
        except ValueError:
            return f"died: signal {-status}"
    return f"failed: exit status {status}"


def pytest_terminal_summary(terminalreporter):
    if _GROUP_RESULTS:
        terminalreporter.write_sep("=", "launcher groups")
        for group, count, status, died_in in _GROUP_RESULTS:
            where = f" (in {died_in})" if status != 0 and died_in else ""
            terminalreporter.write_line(f"{group:6} {count:5} tests  {_describe(status)}{where}")
