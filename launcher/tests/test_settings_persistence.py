import os
from pathlib import Path

import pytest
from qtpy import QtCore


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_direct_beam_ipts_persists():
    from launcher.apps.direct_beam import DirectBeamTab

    tab = DirectBeamTab()
    tab.ipts_edit.setText("36776")
    tab.save_settings()
    tab.deleteLater()

    tab2 = DirectBeamTab()
    assert tab2.ipts_edit.text() == "36776"


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_direct_beam_cd_vals_round_trip():
    from launcher.apps.direct_beam import DirectBeamTab

    tab = DirectBeamTab()
    tab.cd_vals = {"mu_file": "/x/y", "Cd": [5.0, 126.5, 249.5, 499.0], "flip_atten": True}
    tab.save_settings()
    tab.deleteLater()

    tab2 = DirectBeamTab()
    assert tab2.cd_vals == {
        "mu_file": "/x/y",
        "Cd": [5.0, 126.5, 249.5, 499.0],
        "flip_atten": True,
    }


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_direct_beam_corrupt_cd_vals_tolerated():
    from launcher.apps.direct_beam import DirectBeamTab

    s = QtCore.QSettings()
    s.setValue("direct_beam_cd_vals", "not valid json {{")
    s.sync()
    tab = DirectBeamTab()
    assert isinstance(tab.cd_vals, dict)


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_direct_beam_tof_spins_persist():
    """The #155/#156 TOF-rebin spinboxes postdate PR #11; guard them too."""
    from launcher.apps.direct_beam import DirectBeamTab

    tab = DirectBeamTab()
    tab.tofbin_spin.setValue(125.0)
    tab.tofmin_spin.setValue(3000.0)
    tab.tofmax_spin.setValue(55000.0)
    tab.save_settings()
    tab.deleteLater()

    tab2 = DirectBeamTab()
    assert tab2.tofbin_spin.value() == pytest.approx(125.0)
    assert tab2.tofmin_spin.value() == pytest.approx(3000.0)
    assert tab2.tofmax_spin.value() == pytest.approx(55000.0)


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_overplot_saves_ytransform():
    from launcher.apps.overplot import Overplot

    tab = Overplot()
    tab.ytransform_combo.setCurrentText("R*Q^4")
    tab.save_settings()
    s = QtCore.QSettings()
    assert s.value("overplot_ytransform") == "R*Q^4"


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_direct_beam_first_launch_defaults():
    from launcher.apps.direct_beam import DirectBeamTab

    tab = DirectBeamTab()
    assert tab.ipts_edit.text() == ""
    assert tab.tofbin_spin.value() == pytest.approx(50.0)
    assert tab.tofmin_spin.value() == pytest.approx(0.0)
    assert tab.tofmax_spin.value() == pytest.approx(100000.0)
    # exp's deployed first-launch behaviour: the toggle starts unchecked and
    # the manual path fields are therefore editable.
    assert tab.ipts_toggle.isChecked() is False
    assert tab.nexus_edit.isEnabled() is True
    assert tab.savepath_edit.isEnabled() is True


# Every persisted Direct-beam field: a setter, a reader, and a non-default
# value. Parametrized so a broken key is named individually rather than hidden
# behind one composite assertion.
FIELD_SPECS = {
    "run_list": (lambda t, v: t.run_list_edit.setText(v), lambda t: t.run_list_edit.text(), "12345, 12346"),
    "ipts": (lambda t, v: t.ipts_edit.setText(v), lambda t: t.ipts_edit.text(), "36776"),
    "use_ipts_path_structure": (
        lambda t, v: t.ipts_toggle.setChecked(v),
        lambda t: t.ipts_toggle.isChecked(),
        True,
    ),
    "nexus_path": (lambda t, v: t.nexus_edit.setText(v), lambda t: t.nexus_edit.text(), "/tmp/nexus/"),
    "save_path": (lambda t, v: t.savepath_edit.setText(v), lambda t: t.savepath_edit.text(), "/tmp/out/"),
    "save_name": (lambda t, v: t.savename_edit.setText(v), lambda t: t.savename_edit.text(), "db_custom"),
    "DTCcut": (lambda t, v: t.DTCcut_spin.setValue(v), lambda t: t.DTCcut_spin.value(), 7.5),
    "DTCcut_config1": (lambda t, v: t.DTCcut1_spin.setValue(v), lambda t: t.DTCcut1_spin.value(), 8.25),
    "Icut": (lambda t, v: t.Icut_spin.setValue(v), lambda t: t.Icut_spin.value(), 3.5),
    "chopper_cut_offset": (lambda t, v: t.CutOffset_spin.setValue(v), lambda t: t.CutOffset_spin.value(), 2.25),
    "tofbin": (lambda t, v: t.tofbin_spin.setValue(v), lambda t: t.tofbin_spin.value(), 125.0),
    "tofmin": (lambda t, v: t.tofmin_spin.setValue(v), lambda t: t.tofmin_spin.value(), 3000.0),
    "tofmax": (lambda t, v: t.tofmax_spin.setValue(v), lambda t: t.tofmax_spin.value(), 55000.0),
    "y_ROI": (lambda t, v: t.yroi_edit.setText(v), lambda t: t.yroi_edit.text(), "110,160"),
    "x_ROI": (lambda t, v: t.lowres_edit.setText(v), lambda t: t.lowres_edit.text(), "60,200"),
    "plot": (lambda t, v: t.plot_cb.setChecked(v), lambda t: t.plot_cb.isChecked(), False),
    "cd_vals": (
        lambda t, v: setattr(t, "cd_vals", v),
        lambda t: t.cd_vals,
        {"mu_file": "/x", "Cd": [5.0], "flip_atten": True},
    ),
    "mod_vals": (
        lambda t, v: setattr(t, "mod_vals", v),
        lambda t: t.mod_vals,
        {"Chop2_cut_fn": [1.0, 2.0], "dMod": 15500.0, "t0": [3.0, 4.0]},
    ),
}


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
@pytest.mark.parametrize("field", sorted(FIELD_SPECS))
def test_direct_beam_field_round_trip(field):
    """Every persisted key survives a tab restart, named individually."""
    from launcher.apps.direct_beam import DirectBeamTab

    setter, getter, value = FIELD_SPECS[field]
    tab = DirectBeamTab()
    setter(tab, value)
    tab.save_settings()
    tab.deleteLater()

    tab2 = DirectBeamTab()
    read_back = getter(tab2)
    if isinstance(value, float):
        assert read_back == pytest.approx(value)
    else:
        assert read_back == value


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_direct_beam_saves_without_explicit_call():
    """The defense-in-depth wiring, not just save_settings() itself."""
    from launcher.apps.direct_beam import DirectBeamTab

    tab = DirectBeamTab()
    tab.ipts_edit.setText("36776")
    tab.ipts_edit.editingFinished.emit()
    tab.deleteLater()

    tab2 = DirectBeamTab()
    assert tab2.ipts_edit.text() == "36776"


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_direct_beam_settings_reach_disk():
    """A store that never lands on disk persists nothing across a restart."""
    import os

    from launcher.apps.direct_beam import DirectBeamTab

    tab = DirectBeamTab()
    tab.ipts_edit.setText("36776")
    tab.save_settings()
    path = tab.settings.fileName()
    assert os.path.isfile(path), f"settings file not written: {path}"
    with open(path, encoding="utf-8") as handle:
        assert "36776" in handle.read()


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
@pytest.mark.parametrize("stored,expect_enabled", [(True, False), (False, True)])
def test_ipts_mode_applied_on_restore(stored, expect_enabled):
    """B1: the restored toggle must also apply its enable/disable mode.

    Without this the checkbox reads "use IPTS path structure" while the two
    manual path fields stay editable — the scientist types paths that
    _run_create_db never reads and the values are silently discarded.
    """
    from launcher.apps.direct_beam import DirectBeamTab

    s = QtCore.QSettings()
    s.setValue("direct_beam_use_ipts_path_structure", stored)
    s.sync()
    tab = DirectBeamTab()
    assert tab.ipts_toggle.isChecked() is stored
    assert tab.nexus_edit.isEnabled() is expect_enabled
    assert tab.savepath_edit.isEnabled() is expect_enabled


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_ensure_identity_does_not_clobber_existing():
    """Idempotent and non-clobbering: a test fixture's throwaway identity (or a
    second call) must survive, or tabs would drag every store to the
    production identity mid-session."""
    from launcher.app_identity import ensure_identity

    before = (
        QtCore.QCoreApplication.organizationName(),
        QtCore.QCoreApplication.applicationName(),
    )
    ensure_identity()
    ensure_identity()
    after = (
        QtCore.QCoreApplication.organizationName(),
        QtCore.QCoreApplication.applicationName(),
    )
    assert after == before


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_ensure_identity_installs_when_unset():
    """Restores the identity itself: installing the production org/app
    process-wide and leaving it there is the very leak B2 was rejected for —
    a later test would resolve QSettings to the developer's real config."""
    from launcher import app_identity

    saved = (
        QtCore.QCoreApplication.organizationName(),
        QtCore.QCoreApplication.organizationDomain(),
        QtCore.QCoreApplication.applicationName(),
    )
    try:
        QtCore.QCoreApplication.setOrganizationName("")
        QtCore.QCoreApplication.setOrganizationDomain("")
        QtCore.QCoreApplication.setApplicationName("")
        app_identity.ensure_identity()
        assert QtCore.QCoreApplication.organizationName() == app_identity.ORG_NAME
        assert QtCore.QCoreApplication.organizationDomain() == app_identity.ORG_DOMAIN
        assert QtCore.QCoreApplication.applicationName() == app_identity.APP_NAME
    finally:
        QtCore.QCoreApplication.setOrganizationName(saved[0])
        QtCore.QCoreApplication.setOrganizationDomain(saved[1])
        QtCore.QCoreApplication.setApplicationName(saved[2])


def _seed_legacy_store(values):
    """Write *values* into the identity-less store the launcher used before.

    Blanking the identity is process-global, so the restore is in a finally —
    a failure here would otherwise leave every later test nameless.
    """
    saved = (
        QtCore.QCoreApplication.organizationName(),
        QtCore.QCoreApplication.organizationDomain(),
        QtCore.QCoreApplication.applicationName(),
    )
    try:
        QtCore.QCoreApplication.setOrganizationName("")
        QtCore.QCoreApplication.setOrganizationDomain("")
        QtCore.QCoreApplication.setApplicationName("")
        legacy = QtCore.QSettings()
        for key, value in values.items():
            legacy.setValue(key, value)
        legacy.sync()
    finally:
        QtCore.QCoreApplication.setOrganizationName(saved[0])
        QtCore.QCoreApplication.setOrganizationDomain(saved[1])
        QtCore.QCoreApplication.setApplicationName(saved[2])


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_migrate_legacy_settings_copies_and_restores_identity():
    """The migration reads an org-less store, so it mutates the process-global
    identity; it must put it back even though it succeeded."""
    from launcher import app_identity

    identity_before = (
        QtCore.QCoreApplication.organizationName(),
        QtCore.QCoreApplication.organizationDomain(),
        QtCore.QCoreApplication.applicationName(),
    )
    _seed_legacy_store({"overplot_folder": "/legacy/folder", "overplot_xscale": "log"})

    app_identity.migrate_legacy_settings()

    current = QtCore.QSettings()
    assert current.value("overplot_folder") == "/legacy/folder"
    assert current.value("overplot_xscale") == "log"
    identity_after = (
        QtCore.QCoreApplication.organizationName(),
        QtCore.QCoreApplication.organizationDomain(),
        QtCore.QCoreApplication.applicationName(),
    )
    assert identity_after == identity_before


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_migrate_legacy_settings_never_overwrites():
    """One-shot: a value the user has already set in the new store wins."""
    from launcher import app_identity

    _seed_legacy_store({"overplot_folder": "/legacy/folder"})
    current = QtCore.QSettings()
    current.setValue("overplot_folder", "/current/folder")
    current.sync()
    app_identity.migrate_legacy_settings()
    assert QtCore.QSettings().value("overplot_folder") == "/current/folder"


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_migration_carries_every_launcher_module():
    """Installing the identity repoints the store for the whole process, so a
    migration narrower than the key set is a one-time mass-forget across tabs
    that have nothing to do with this slug. One key per affected module."""
    from launcher import app_identity

    probes = {
        "30Hz_template": "/legacy/30hz.xml",          # dynamic_30Hz
        "60Hz_template": "/legacy/60hz.xml",          # dynamic_60Hz
        "settings_Spath": "/legacy/batch/out",        # file_batch
        "offspec_output_dir": "/legacy/offspec",      # off_spec
        "overplot_folder": "/legacy/overplot",        # overplot
        "quick_output_dir": "/legacy/quick",          # quick_reduce
        "reduction_template": "/legacy/reduction.xml",  # reduction
        "refracted_material": "Si",                   # refracted
        "sld_composition": "H2O",                     # sld_calculator
        "template_Spath": "/legacy/template/out",     # template_batch
        "xrr_data_file": "/legacy/xrr.dat",           # xrr
    }
    assert set(probes) <= set(app_identity.LEGACY_KEYS)
    _seed_legacy_store(probes)

    app_identity.migrate_legacy_settings()
    current = QtCore.QSettings()
    for key, value in probes.items():
        assert current.value(key) == value, f"{key} lost across the identity switch"

    # Second run is a no-op: the sentinel records that the copy happened, so a
    # value the user clears afterwards does not come back from the dead.
    current.remove("overplot_folder")
    current.sync()
    app_identity.migrate_legacy_settings()
    assert QtCore.QSettings().value("overplot_folder") is None


# The 16 mutating signals wired for defense-in-depth saves. Table-driven so a
# dropped connection names itself; one explicit save_settings() test cannot see
# a missing connection at all.
SAVE_CONNECTIONS = [
    ("run_list_edit", "editingFinished", lambda t: t.run_list_edit.setText("999"), "direct_beam_run_list", "999"),
    ("ipts_edit", "editingFinished", lambda t: t.ipts_edit.setText("36776"), "direct_beam_ipts", "36776"),
    ("nexus_edit", "editingFinished", lambda t: t.nexus_edit.setText("/n/"), "direct_beam_nexus_path", "/n/"),
    ("savepath_edit", "editingFinished", lambda t: t.savepath_edit.setText("/s/"), "direct_beam_save_path", "/s/"),
    ("savename_edit", "editingFinished", lambda t: t.savename_edit.setText("nm"), "direct_beam_save_name", "nm"),
    ("yroi_edit", "editingFinished", lambda t: t.yroi_edit.setText("1,2"), "direct_beam_y_ROI", "1,2"),
    ("lowres_edit", "editingFinished", lambda t: t.lowres_edit.setText("3,4"), "direct_beam_x_ROI", "3,4"),
    ("DTCcut_spin", None, lambda t: t.DTCcut_spin.setValue(7.5), "direct_beam_DTCcut", 7.5),
    ("DTCcut1_spin", None, lambda t: t.DTCcut1_spin.setValue(8.5), "direct_beam_DTCcut_config1", 8.5),
    ("Icut_spin", None, lambda t: t.Icut_spin.setValue(3.5), "direct_beam_Icut", 3.5),
    ("CutOffset_spin", None, lambda t: t.CutOffset_spin.setValue(2.5), "direct_beam_chopper_cut_offset", 2.5),
    ("tofbin_spin", None, lambda t: t.tofbin_spin.setValue(125.0), "direct_beam_tofbin", 125.0),
    ("tofmin_spin", None, lambda t: t.tofmin_spin.setValue(3000.0), "direct_beam_tofmin", 3000.0),
    ("tofmax_spin", None, lambda t: t.tofmax_spin.setValue(55000.0), "direct_beam_tofmax", 55000.0),
    ("plot_cb", None, lambda t: t.plot_cb.setChecked(False), "direct_beam_plot", False),
    ("ipts_toggle", None, lambda t: t.ipts_toggle.setChecked(True), "direct_beam_use_ipts_path_structure", True),
]


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
@pytest.mark.parametrize("widget,signal,mutate,key,expected", SAVE_CONNECTIONS, ids=[c[0] for c in SAVE_CONNECTIONS])
def test_mutating_signal_persists_without_explicit_save(widget, signal, mutate, key, expected):
    """Each wired signal must reach save_settings on its own."""
    from launcher.apps.direct_beam import DirectBeamTab

    tab = DirectBeamTab()
    mutate(tab)
    if signal is not None:
        getattr(getattr(tab, widget), signal).emit()
    stored = QtCore.QSettings().value(key)
    if isinstance(expected, float):
        assert float(stored) == pytest.approx(expected)
    elif isinstance(expected, bool):
        assert stored in (expected, str(expected).lower())
    else:
        assert stored == expected


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_dialog_accept_persists_cd_and_mod_values(monkeypatch):
    """The dialog saves are the only persistence path for cd_vals/mod_vals.

    exec_ is patched here rather than left to no_qmessagebox: on this branch's
    base the two dialogs still wrap a QMessageBox around an inner QDialog and
    forward exec_ to it, so the fixture's statics-only patch does not reach the
    call and the test hangs until the timeout backstop kills it (observed).
    The harness-hardening slug moves that patch to QDialog, but it has not
    merged into exp yet — a test-local patch works on either base.
    """
    import json

    from launcher.apps import direct_beam as direct_beam_module
    from launcher.apps.direct_beam import DirectBeamTab

    monkeypatch.setattr(direct_beam_module.CdSettingsDialog, "exec_", lambda _self: 1)
    monkeypatch.setattr(direct_beam_module.ModeratorDialog, "exec_", lambda _self: 1)

    tab = DirectBeamTab()
    tab._open_cd_settings()
    tab._open_mod_settings()
    stored_cd = json.loads(QtCore.QSettings().value("direct_beam_cd_vals"))
    stored_mod = json.loads(QtCore.QSettings().value("direct_beam_mod_vals"))
    assert isinstance(stored_cd, dict) and "Cd" in stored_cd
    assert isinstance(stored_mod, dict) and "dMod" in stored_mod


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_live_toggle_applies_ipts_mode():
    """The click path, not just the restore path: _ipts_toggled must apply the
    enable/disable mode itself."""
    from launcher.apps.direct_beam import DirectBeamTab

    tab = DirectBeamTab()
    assert tab.nexus_edit.isEnabled() is True
    tab.ipts_toggle.setChecked(True)
    assert tab.nexus_edit.isEnabled() is False
    assert tab.savepath_edit.isEnabled() is False
    tab.ipts_toggle.setChecked(False)
    assert tab.nexus_edit.isEnabled() is True
    assert tab.savepath_edit.isEnabled() is True


def test_settings_survive_a_real_process_restart(tmp_path):
    """Two processes, because an in-process 'restart' proves nothing.

    Qt serves a second in-process QSettings from the first's cached QConfFile,
    so every same-process round-trip skips the read-back coercion entirely:
    inverting the string branch of _bool leaves the whole suite green while a
    real restart returns booleans inverted. Only a fresh interpreter reading
    the file from disk exercises that path — do not "simplify" this back into
    the isolated_qapp fixture.
    """
    import json
    import subprocess
    import sys

    home = tmp_path / "home"
    home.mkdir()
    env = dict(os.environ, HOME=str(home), XDG_CONFIG_HOME=str(home / ".config"), QT_QPA_PLATFORM="offscreen")

    save = (
        "from qtpy.QtWidgets import QApplication\n"
        "app = QApplication([])\n"
        "from launcher.app_identity import ensure_identity\n"
        "ensure_identity()\n"
        "from launcher.apps.direct_beam import DirectBeamTab\n"
        "tab = DirectBeamTab()\n"
        "tab.run_list_edit.setText('101,102')\n"
        "tab.ipts_edit.setText('36776')\n"
        "tab.ipts_toggle.setChecked(True)\n"
        "tab.savename_edit.setText('db_custom')\n"
        "tab.DTCcut_spin.setValue(7.5)\n"
        "tab.tofbin_spin.setValue(125.0)\n"
        "tab.plot_cb.setChecked(False)\n"
        "tab.cd_vals = {'Cd': [5.0], 'flip_atten': True}\n"
        "tab.save_settings()\n"
    )
    read = (
        "import json\n"
        "from qtpy.QtWidgets import QApplication\n"
        "app = QApplication([])\n"
        "from launcher.app_identity import ensure_identity\n"
        "ensure_identity()\n"
        "from launcher.apps.direct_beam import DirectBeamTab\n"
        "tab = DirectBeamTab()\n"
        "print('RESULT' + json.dumps({\n"
        "    'run_list': tab.run_list_edit.text(),\n"
        "    'ipts': tab.ipts_edit.text(),\n"
        "    'toggle': tab.ipts_toggle.isChecked(),\n"
        "    'savename': tab.savename_edit.text(),\n"
        "    'dtc': tab.DTCcut_spin.value(),\n"
        "    'tofbin': tab.tofbin_spin.value(),\n"
        "    'plot': tab.plot_cb.isChecked(),\n"
        "    'cd': tab.cd_vals,\n"
        "    'nexus_enabled': tab.nexus_edit.isEnabled(),\n"
        "}))\n"
    )
    root = str(Path(__file__).resolve().parents[2])
    for source in (save, read):
        proc = subprocess.run(
            [sys.executable, "-c", source], cwd=root, env=env, capture_output=True, timeout=180, check=False
        )
        assert proc.returncode == 0, proc.stderr.decode(errors="replace")[-2000:]
        last = proc.stdout

    payload = json.loads(last.decode().split("RESULT", 1)[1].strip())
    assert payload == {
        "run_list": "101,102",
        "ipts": "36776",
        "toggle": True,
        "savename": "db_custom",
        "dtc": 7.5,
        "tofbin": 125.0,
        "plot": False,
        "cd": {"Cd": [5.0], "flip_atten": True},
        "nexus_enabled": False,
    }


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_direct_beam_plot_does_not_destroy_stored_ytransform(tmp_path):
    """S2<->S3 interaction: plotting direct-beam data must not erase the
    user's stored R*Q^4 preference.

    #17 disables the R*Q^4 entry off-reflectivity and resets the combo to
    "None" — correct for the axes. But save_settings then writes that reset
    value, so one direct-beam plot silently destroys the stored preference.
    The display reset is intended; the persistence of it is not.
    """
    from launcher.apps.overplot import Overplot

    (tmp_path / "db.txt").write_text("# columns = lambda intensity error\n1.0 100.0 1.0\n2.0 200.0 1.4\n")
    s = QtCore.QSettings()
    s.setValue("overplot_ytransform", "R*Q^4")
    s.sync()

    tab = Overplot()
    assert tab.ytransform_combo.currentText() == "R*Q^4", "precondition: stored value restored"
    tab.folder_edit.setText(str(tmp_path))
    tab.populate_file_list(str(tmp_path))
    for i in range(tab.file_list.count()):
        tab.file_list.item(i).setCheckState(QtCore.Qt.Checked)
    tab.plot_selected()

    assert QtCore.QSettings().value("overplot_ytransform") == "R*Q^4", (
        "a direct-beam plot destroyed the stored transform preference"
    )


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_stored_ytransform_survives_a_round_trip_through_direct_beam(tmp_path):
    """The same, one plot later — the case a disabled-only guard misses.

    Re-enabling R*Q^4 does not restore the combo, so the next reflectivity
    plot saves the residual "None" with the entry *enabled*. A guard that only
    skips the write while disabled therefore delays the destruction by one
    plot rather than preventing it.
    """
    from launcher.apps.overplot import Overplot

    (tmp_path / "db.txt").write_text("# columns = lambda intensity error\n1.0 100.0 1.0\n2.0 200.0 1.4\n")
    (tmp_path / "refl.txt").write_text("# Q [1/Angstrom] R dR dQ\n0.01 0.5 0.05 0.001\n0.02 0.4 0.04 0.001\n")
    s = QtCore.QSettings()
    s.setValue("overplot_ytransform", "R*Q^4")
    s.sync()

    tab = Overplot()
    tab.folder_edit.setText(str(tmp_path))
    tab.populate_file_list(str(tmp_path))

    def _check_only(name):
        for i in range(tab.file_list.count()):
            item = tab.file_list.item(i)
            item.setCheckState(QtCore.Qt.Checked if item.text() == name else QtCore.Qt.Unchecked)

    _check_only("db.txt")
    tab.plot_selected()
    _check_only("refl.txt")
    tab.plot_selected()

    assert QtCore.QSettings().value("overplot_ytransform") == "R*Q^4", (
        "the stored transform did not survive a direct-beam plot followed by a reflectivity plot"
    )


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_user_choice_of_none_is_not_resurrected():
    """The set-aside value must not override a choice the user actually makes.

    Restoring a displaced R*Q^4 is right after a programmatic reset; it would be
    wrong after the user selects None themselves, so an explicit pick clears it.
    """
    from launcher.apps.overplot import Overplot

    s = QtCore.QSettings()
    s.setValue("overplot_ytransform", "R*Q^4")
    s.sync()

    tab = Overplot()
    tab._set_rq4_enabled(False)
    assert tab._rq4_displaced_choice == "R*Q^4", "precondition: the choice was set aside"

    tab._on_ytransform_chosen(tab.ytransform_combo.currentIndex())  # the user picks None
    tab.save_settings()

    assert QtCore.QSettings().value("overplot_ytransform") == "None"
    tab._set_rq4_enabled(True)
    assert tab.ytransform_combo.currentText() == "None", "a user's own choice must not be overridden"
