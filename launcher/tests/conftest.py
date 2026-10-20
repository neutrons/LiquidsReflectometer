import os

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import pytest
from qtpy import QtCore, QtWidgets


@pytest.fixture
def isolated_qapp(tmp_path, monkeypatch):
    """Fresh QApplication with throwaway QSettings backing.

    Uses an isolated org/app so the test never reads or writes the
    user's real new_launcher config.
    """
    QtCore.QCoreApplication.setOrganizationName(f"test-org-{tmp_path.name}")
    QtCore.QCoreApplication.setOrganizationDomain("example.test")
    QtCore.QCoreApplication.setApplicationName(f"test-app-{tmp_path.name}")
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    if QtWidgets.QApplication.instance() is None:
        app = QtWidgets.QApplication([])
    else:
        app = QtWidgets.QApplication.instance()
    yield app


@pytest.fixture
def no_qmessagebox(monkeypatch):
    """Neutralize modal QMessageBox calls that hang under offscreen Qt."""
    from qtpy import QtWidgets

    for name in ("warning", "information", "critical", "question"):
        monkeypatch.setattr(
            QtWidgets.QMessageBox,
            name,
            staticmethod(lambda *_args, **_kwargs: None),
        )
