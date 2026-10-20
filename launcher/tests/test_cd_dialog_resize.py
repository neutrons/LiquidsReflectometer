# launcher/tests/test_cd_dialog_resize.py
import pytest


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_cd_dialog_is_qdialog():
    from qtpy.QtWidgets import QDialog, QMessageBox

    from launcher.apps.direct_beam import CdSettingsDialog

    dlg = CdSettingsDialog()
    assert isinstance(dlg, QDialog)
    assert not isinstance(dlg, QMessageBox)


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_cd_edit_minimum_width():
    from launcher.apps.direct_beam import CdSettingsDialog

    dlg = CdSettingsDialog(defaults={"mu_file": "", "Cd": [5, 126.5, 249.5, 499.0], "flip_atten": False})
    assert dlg.cd_edit.minimumWidth() >= 320


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_cd_dialog_size_hint():
    from launcher.apps.direct_beam import CdSettingsDialog

    dlg = CdSettingsDialog(defaults={"mu_file": "/x/y", "Cd": [5, 126.5, 249.5, 499.0], "flip_atten": False})
    fm = dlg.cd_edit.fontMetrics()
    needed = fm.horizontalAdvance("5, 126.5, 249.5, 499.0") + 16
    assert dlg.sizeHint().width() >= needed


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_moderator_dialog_is_qdialog():
    from qtpy.QtWidgets import QDialog, QMessageBox

    from launcher.apps.direct_beam import ModeratorDialog

    dlg = ModeratorDialog()
    assert isinstance(dlg, QDialog)
    assert not isinstance(dlg, QMessageBox)


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_get_values_round_trip():
    from launcher.apps.direct_beam import CdSettingsDialog

    dlg = CdSettingsDialog(defaults={"mu_file": "/x/y", "Cd": [5, 126.5, 249.5, 499.0], "flip_atten": True})
    vals = dlg.get_values()
    assert vals["mu_file"] == "/x/y"
    assert vals["Cd"] == [5.0, 126.5, 249.5, 499.0]
    assert vals["flip_atten"] is True


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_reset_defaults():
    from launcher.apps.direct_beam import CdSettingsDialog

    initial = {"mu_file": "/initial", "Cd": [1.0], "flip_atten": False}
    working = {"mu_file": "/working", "Cd": [9.0], "flip_atten": True}
    dlg = CdSettingsDialog(defaults=working, initial_defaults=initial)
    dlg._reset_defaults()
    vals = dlg.get_values()
    assert vals["mu_file"] == "/initial"
    assert vals["Cd"] == [1.0]
    assert vals["flip_atten"] is False
