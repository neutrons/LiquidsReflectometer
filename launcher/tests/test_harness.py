# launcher/tests/test_harness.py
import subprocess
import sys

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


def test_qsettings_file_lands_in_tmp(isolated_qapp, tmp_path):
    """The roundtrip test passes even with isolation broken; this one does not.

    Qt caches the settings root at the first QSettings construction in a
    process, so an XDG_CONFIG_HOME redirect alone is honoured only for the
    first fixture use. Assert the actual file path.
    """
    settings = QtCore.QSettings()
    settings.setValue("harness/probe", "x")
    settings.sync()
    assert QtCore.QSettings().fileName().startswith(str(tmp_path))


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_instance_exec_neutralized():
    """Instance .exec_() routes around a fixture that patches statics only."""
    box = QtWidgets.QMessageBox()
    assert box.exec_() == QtWidgets.QMessageBox.Ok


def test_timeout_backstop_fires(tmp_path):
    """Self-test: if pytest-timeout is missing or the method regresses, this
    test goes red instead of a future slug hanging for hours."""
    hanging = tmp_path / "test_hangs.py"
    hanging.write_text("def test_hangs():\n    while True:\n        pass\n")
    proc = subprocess.run(
        [sys.executable, "-m", "pytest", str(hanging), "--timeout=5", "--timeout-method=thread", "-p", "no:cacheprovider"],
        capture_output=True,
        timeout=60,
        check=False,
    )
    assert proc.returncode != 0
    assert b"Timeout" in proc.stdout + proc.stderr
