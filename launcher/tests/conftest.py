import atexit
import os
import shutil
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
    """Arm a timeout for the tests this conftest governs.

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
