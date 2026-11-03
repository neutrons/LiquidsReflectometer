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
        for widget in QtWidgets.QApplication.topLevelWidgets():
            try:
                widget.close()
                widget.deleteLater()
            except RuntimeError:
                # Already destroyed on the C++ side; nothing to drain.
                pass
        app.processEvents()
        # deleteLater only queues; without this the objects are never freed and
        # the "drain" drains nothing.
        app.sendPostedEvents(None, QtCore.QEvent.DeferredDelete)


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
