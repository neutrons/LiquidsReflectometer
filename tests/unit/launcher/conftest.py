"""Fixtures for launcher UI tests using offscreen rendering."""

import os
import sys

import pytest

# Force offscreen rendering before any Qt imports
os.environ["QT_QPA_PLATFORM"] = "offscreen"


@pytest.fixture(scope="session")
def qapp():
    """Create a QApplication instance for the test session."""
    from qtpy.QtCore import QCoreApplication
    from qtpy.QtWidgets import QApplication

    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    # Use a test-specific settings scope to avoid polluting/reading real launcher settings
    QCoreApplication.setOrganizationName("lr_reduction_tests")
    QCoreApplication.setApplicationName("test_launcher")
    yield app


@pytest.fixture(autouse=True)
def clear_tmpl_settings(qapp):
    """Clear template-related QSettings before and after each test to avoid pollution."""
    from qtpy import QtCore

    settings = QtCore.QSettings()
    for key in settings.allKeys():
        if key.startswith("tmpl_"):
            settings.remove(key)
    settings.sync()
    yield
    for key in settings.allKeys():
        if key.startswith("tmpl_"):
            settings.remove(key)
    settings.sync()
