# launcher/tests/test_harness.py
import pytest
from qtpy import QtCore, QtWidgets


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_no_qmessagebox_fixture_neutralizes_modals():
    # Without the fixture this call blocks forever under offscreen Qt
    # (the 10-hour-orphan class); with it, it returns immediately.
    assert QtWidgets.QMessageBox.warning(None, "t", "must not block") is None


@pytest.mark.usefixtures("isolated_qapp")
def test_isolated_qsettings_roundtrip():
    settings = QtCore.QSettings()
    settings.setValue("harness/probe", "x")
    settings.sync()
    assert settings.value("harness/probe") == "x"
    assert "test-org" in QtCore.QCoreApplication.organizationName()
