import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from qtpy import QtCore, QtWidgets


@pytest.fixture
def isolated_qapp(tmp_path, monkeypatch):
    """Fresh QApplication with throwaway QSettings backing.

    Isolation is enforced by QSettings.setPath, which is consulted on every
    construction. The XDG_CONFIG_HOME redirect alone is not enough: Qt caches
    the settings root at the first QSettings construction in a process, so a
    later test would inherit the first test's directory — or, if a launcher
    module constructs QSettings before any fixture runs, the developer's real
    ~/.config for the whole session.
    """
    prev_org = QtCore.QCoreApplication.organizationName()
    prev_domain = QtCore.QCoreApplication.organizationDomain()
    prev_app = QtCore.QCoreApplication.applicationName()
    prev_format = QtCore.QSettings.defaultFormat()

    QtCore.QSettings.setDefaultFormat(QtCore.QSettings.IniFormat)
    QtCore.QSettings.setPath(QtCore.QSettings.IniFormat, QtCore.QSettings.UserScope, str(tmp_path))
    QtCore.QCoreApplication.setOrganizationName(f"test-org-{tmp_path.name}")
    QtCore.QCoreApplication.setOrganizationDomain("example.test")
    QtCore.QCoreApplication.setApplicationName(f"test-app-{tmp_path.name}")
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    if QtWidgets.QApplication.instance() is None:
        app = QtWidgets.QApplication([])
    else:
        app = QtWidgets.QApplication.instance()
    yield app

    # One QApplication is shared by the whole session, so widgets a test leaves
    # behind outlive it; drain them before restoring the identity.
    for widget in QtWidgets.QApplication.topLevelWidgets():
        widget.close()
        widget.deleteLater()
    app.processEvents()
    QtCore.QCoreApplication.setOrganizationName(prev_org)
    QtCore.QCoreApplication.setOrganizationDomain(prev_domain)
    QtCore.QCoreApplication.setApplicationName(prev_app)
    QtCore.QSettings.setDefaultFormat(prev_format)


@pytest.fixture
def no_qmessagebox(monkeypatch):
    """Neutralize modal QMessageBox calls that hang under offscreen Qt.

    Covers both the static convenience methods and instance .exec_()/.exec(),
    which route around a statics-only patch.
    """
    from qtpy import QtWidgets

    for name in ("warning", "information", "critical", "question"):
        monkeypatch.setattr(
            QtWidgets.QMessageBox,
            name,
            staticmethod(lambda *_args, **_kwargs: None),
        )
    for name in ("exec_", "exec"):
        if hasattr(QtWidgets.QMessageBox, name):
            monkeypatch.setattr(
                QtWidgets.QMessageBox,
                name,
                lambda *_args, **_kwargs: QtWidgets.QMessageBox.Ok,
            )


@pytest.fixture
def no_qfiledialog(monkeypatch):
    """Neutralize modal QFileDialog calls (26 sites in launcher/).

    Opt-in: a test that wants a chosen path patches the specific method
    itself; this fixture only guarantees nothing blocks.
    """
    from qtpy import QtWidgets

    for name, result in (
        ("getOpenFileName", ("", "")),
        ("getSaveFileName", ("", "")),
        ("getOpenFileNames", ([], "")),
        ("getExistingDirectory", ""),
    ):
        monkeypatch.setattr(
            QtWidgets.QFileDialog,
            name,
            staticmethod(lambda *_args, _result=result, **_kwargs: _result),
        )


def pytest_collection_modifyitems(items):
    """Arm a timeout for launcher tests regardless of how pytest was invoked.

    The pixi task passes --timeout explicitly, but the inner loop a slug author
    actually types (`pytest launcher/tests/test_x.py -k ...`) would otherwise
    have no timeout at all — the 10-hour-orphan class. This conftest governs
    only launcher/tests, so the reduction suite's timings are untouched.
    """
    for item in items:
        if item.get_closest_marker("timeout") is None:
            item.add_marker(pytest.mark.timeout(120))
