import os
import shutil
from pathlib import Path

import pytest


def _seed_folder(folder, *fixture_names):
    """Copy named fixture files from launcher/tests/data into folder."""
    fixture_root = os.path.join(os.path.dirname(__file__), "data")
    for name in fixture_names:
        shutil.copy(os.path.join(fixture_root, name), os.path.join(folder, name))


def _select_files(widget, *filenames):
    from qtpy import QtCore

    targets = set(filenames)
    for i in range(widget.file_list.count()):
        item = widget.file_list.item(i)
        if item.text() in targets:
            item.setCheckState(QtCore.Qt.Checked)


def _rq4_is_enabled(combo):
    from qtpy import QtCore

    idx = combo.findText("R*Q^4")
    if idx < 0:
        return False
    item = combo.model().item(idx)
    return bool(item.flags() & QtCore.Qt.ItemIsEnabled)


def test_classify_db_fixture(tmp_path):
    from launcher.apps.overplot import classify_file

    p = tmp_path / "DB_test.txt"
    p.write_text("# Header\n# columns = lambda intensity error\n1.0 100.0 1.0\n2.0 200.0 1.4\n")
    assert classify_file(str(p)) == "direct_beam"


def test_classify_refl_fixture(tmp_path):
    from launcher.apps.overplot import classify_file

    p = tmp_path / "REFL_test.txt"
    p.write_text("# header line A\n# Q [1/Angstrom] R dR dQ [FWHM]\n0.01 0.5 0.05 0.001\n0.02 0.4 0.04 0.001\n")
    assert classify_file(str(p)) == "reflectivity"


def test_classify_unknown_no_header(tmp_path):
    from launcher.apps.overplot import classify_file

    p = tmp_path / "garbled.txt"
    p.write_text("0.01 0.5\n0.02 0.4\n")
    assert classify_file(str(p)) == "unknown"


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_overplot_db_labels(tmp_path):
    from launcher.apps.overplot import Overplot

    _seed_folder(tmp_path, "db_fixture.txt")
    w = Overplot()
    w.folder_edit.setText(str(tmp_path))
    w.populate_file_list(str(tmp_path))
    _select_files(w, "db_fixture.txt")
    w.plot_selected()
    ax = w.figure.axes[0]
    assert ax.get_xlabel() == "λ [Å]"
    assert ax.get_ylabel() == "I"


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_overplot_refl_labels(tmp_path):
    from launcher.apps.overplot import Overplot

    _seed_folder(tmp_path, "refl_fixture.txt")
    w = Overplot()
    w.folder_edit.setText(str(tmp_path))
    w.populate_file_list(str(tmp_path))
    _select_files(w, "refl_fixture.txt")
    w.plot_selected()
    ax = w.figure.axes[0]
    assert ax.get_xlabel() == "Q [1/Å]"
    assert ax.get_ylabel() == "R"


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_rq4_disabled_for_db(tmp_path):
    from launcher.apps.overplot import Overplot

    _seed_folder(tmp_path, "db_fixture.txt")
    w = Overplot()
    w.folder_edit.setText(str(tmp_path))
    w.populate_file_list(str(tmp_path))
    _select_files(w, "db_fixture.txt")
    w.plot_selected()
    assert not _rq4_is_enabled(w.ytransform_combo)


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_mixed_selection(tmp_path):
    from launcher.apps.overplot import Overplot

    _seed_folder(tmp_path, "db_fixture.txt", "refl_fixture.txt")
    w = Overplot()
    w.folder_edit.setText(str(tmp_path))
    w.populate_file_list(str(tmp_path))
    _select_files(w, "db_fixture.txt", "refl_fixture.txt")
    w.plot_selected()
    assert not _rq4_is_enabled(w.ytransform_combo)
    ax = w.figure.axes[0]
    assert ax.get_xlabel() == "x"
    assert ax.get_ylabel() == "y"


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_mode_persists():
    from launcher.apps.overplot import Overplot

    w = Overplot()
    w.plot_mode_combo.setCurrentText("Direct Beam")
    w.save_settings()
    w.deleteLater()

    w2 = Overplot()
    assert w2.plot_mode_combo.currentText() == "Direct Beam"


# The reduction writers put the reflectivity marker after a metadata preamble,
# so these tracked files carry it at lines 13-20 — outside any fixed-size
# window. Globbed rather than hand-listed: every reference_*.txt is tracked,
# and the hand-list had already drifted (it missed reference_fbck.txt). Their
# REFL_*.txt siblings are deliberately excluded — those are gitignored
# test-run byproducts and would collect zero cases on a fresh clone.
def _data_dir():
    return Path(__file__).resolve().parents[2] / "tests" / "data"


def _repo_data(name):
    return _data_dir() / name


REAL_REFLECTIVITY_FILES = tuple(sorted(p.name for p in _data_dir().glob("reference_*.txt")))


@pytest.mark.parametrize("name", REAL_REFLECTIVITY_FILES)
def test_classify_real_reflectivity_corpus(name):
    """B1: the formats the tab actually consumes, not a synthetic stand-in."""
    from launcher.apps.overplot import classify_file

    path = _repo_data(name)
    assert path.is_file(), f"tracked fixture missing: {path}"
    assert classify_file(str(path)) == "reflectivity"


def test_classify_save_reduced_data_format(tmp_path):
    """B1: save_reduced_data.py is a second writer with a different marker.

    The header comes from the writer itself rather than a hand-written
    imitation, so a change to its format shows up here instead of silently
    diverging from what the tab must classify.
    """
    from types import SimpleNamespace

    import numpy as np

    from launcher.apps.overplot import classify_file
    from lr_reduction.save_reduced_data import _build_header

    config = SimpleNamespace(
        RBnum="1234",
        DBname="db.nxs",
        method_per_run="m",
        Normalize=True,
        AutoScale=False,
        ScaleFactor=1.0,
        LambdaMinUse=2.0,
        LambdaMaxUse=15.0,
    )
    logs = {"ths": 1.0, "thi": 1.0, "ThCen": 1.0, "title": "fixture"}
    header = _build_header(config_header=config, log_values=logs)
    p = tmp_path / "second_writer.dat"
    # Written exactly as the writer does — savetxt is what prefixes '# '.
    np.savetxt(str(p), np.array([[0.01, 0.5, 0.05, 0.001]]), header=header, delimiter="\t")
    assert classify_file(str(p)) == "reflectivity"


def test_classify_deep_preamble(tmp_path):
    """Pins the window removal: the marker sits well below any fixed count."""
    from launcher.apps.overplot import classify_file

    preamble = "".join(f"# metadata line {i} — Å units\n" for i in range(40))
    p = tmp_path / "deep.txt"
    p.write_text(preamble + "# Q [1/Angstrom] R dR dQ [FWHM]\n0.01 0.5 0.05 0.001\n")
    assert classify_file(str(p)) == "reflectivity"


def test_classify_marker_in_data_body_is_unknown(tmp_path):
    """Pins the termination rule: a marker below the header block is not a
    header, and scanning on into the data would misclassify it."""
    from launcher.apps.overplot import classify_file

    p = tmp_path / "late_marker.txt"
    p.write_text("# plain header\n0.01 0.5\n0.02 0.4\n# Q [1/Angstrom] R dR dQ\n0.03 0.3\n")
    assert classify_file(str(p)) == "unknown"


def test_classify_non_utf8_bytes_still_classifies(tmp_path):
    """The preamble carries Å; a locale-dependent decode would reproduce B1."""
    from launcher.apps.overplot import classify_file

    p = tmp_path / "latin1.txt"
    p.write_bytes(
        "# Datafile: lambda range 2\u00c5 to 15\u00c5\n".encode("latin-1")
        + b"# Q [1/Angstrom] R dR dQ [FWHM]\n0.01 0.5 0.05 0.001\n"
    )
    assert classify_file(str(p)) == "reflectivity"


def test_classify_binary_file_is_unknown(tmp_path):
    """B1: a non-text file must classify unknown, not raise."""
    from launcher.apps.overplot import classify_file

    p = tmp_path / "blob.dat"
    p.write_bytes(b"\x00\x01\x02\xff\xfe")
    assert classify_file(str(p)) == "unknown"


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_unknown_is_a_kind_in_homogeneity(tmp_path):
    """B2: one recognized DB file must not drag headerless files onto λ/I."""
    from launcher.apps.overplot import Overplot

    _seed_folder(tmp_path, "db_fixture.txt")
    (tmp_path / "headerless.dat").write_text("0.01 0.5\n0.02 0.4\n")
    w = Overplot()
    w.folder_edit.setText(str(tmp_path))
    w.populate_file_list(str(tmp_path))
    _select_files(w, "db_fixture.txt", "headerless.dat")
    w.plot_selected()
    assert w._last_sel_mode == "mixed"
    ax = w.figure.axes[0]
    assert ax.get_xlabel() == "x"
    assert ax.get_ylabel() == "y"


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_rq4_reenabled_for_real_reflectivity(tmp_path):
    """B1 consequence: the primary case must not be worse than before the slug."""
    from launcher.apps.overplot import Overplot

    shutil.copy(_repo_data("reference_rq.txt"), tmp_path / "reference_rq.txt")
    w = Overplot()
    w.folder_edit.setText(str(tmp_path))
    w.populate_file_list(str(tmp_path))
    _select_files(w, "reference_rq.txt")
    w.plot_selected()
    assert w._last_sel_mode == "reflectivity"
    assert _rq4_is_enabled(w.ytransform_combo)
    ax = w.figure.axes[0]
    assert ax.get_xlabel() == "Q [1/Å]"


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_popout_path_uses_resolved_labels(tmp_path, monkeypatch):
    """Kills the named revert: popout labels hardcoded back to 'x'/'y'."""
    from launcher.apps import overplot as overplot_module
    from launcher.apps.overplot import Overplot

    class _FakeAx:
        def __init__(self):
            self.xlabel = None
            self.ylabel = None

        def errorbar(self, *_a, **_k):
            pass

        def plot(self, *_a, **_k):
            pass

        def set_yscale(self, *_a, **_k):
            pass

        def set_xscale(self, *_a, **_k):
            pass

        def legend(self, *_a, **_k):
            pass

        def set_xlabel(self, value):
            self.xlabel = value

        def set_ylabel(self, value):
            self.ylabel = value

    class _FakePlt:
        def __init__(self):
            self.ax = _FakeAx()

        def figure(self, *_a, **_k):
            pass

        def gca(self):
            return self.ax

        def tight_layout(self, *_a, **_k):
            pass

        def show(self, *_a, **_k):
            pass

        def close(self, *_a, **_k):
            pass

    fake = _FakePlt()
    monkeypatch.setattr(overplot_module, "plt", fake)
    shutil.copy(_repo_data("reference_rq.txt"), tmp_path / "reference_rq.txt")
    w = Overplot()
    w.folder_edit.setText(str(tmp_path))
    w.populate_file_list(str(tmp_path))
    _select_files(w, "reference_rq.txt")
    w.canvas = None
    w.plot_selected()
    assert fake.ax.xlabel == "Q [1/Å]"
    assert fake.ax.ylabel == "R"


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_unreachable_saved_folder_is_not_erased():
    """B4: starting before /SNS mounts must not destroy the stored path."""
    from qtpy import QtCore

    from launcher.apps.overplot import Overplot

    s = QtCore.QSettings()
    s.setValue("overplot_folder", "/nonexistent/xyz")
    s.sync()
    Overplot()
    assert QtCore.QSettings().value("overplot_folder") == "/nonexistent/xyz"


@pytest.mark.usefixtures("isolated_qapp")
def test_override_mismatch_warns(tmp_path, monkeypatch):
    """An explicit mode that disagrees with the headers must say so rather than
    draw the files under the wrong axes silently."""
    from qtpy.QtWidgets import QMessageBox

    from launcher.apps.overplot import Overplot

    _seed_folder(tmp_path, "db_fixture.txt")
    captured = {}
    monkeypatch.setattr(QMessageBox, "warning", lambda *a, **_k: captured.setdefault("warn", a))
    monkeypatch.setattr(QMessageBox, "information", lambda *_a, **_k: None)
    w = Overplot()
    w.folder_edit.setText(str(tmp_path))
    w.populate_file_list(str(tmp_path))
    _select_files(w, "db_fixture.txt")
    w.plot_mode_combo.setCurrentText("Reflectivity")
    w.plot_selected()
    assert "warn" in captured
    assert "db_fixture.txt" in captured["warn"][2]


@pytest.mark.usefixtures("isolated_qapp", "no_qmessagebox")
def test_mixed_labels_name_the_unrecognized_file(tmp_path):
    """In a mixed selection every curve is labelled with what it was taken to
    be, including the unrecognized file that broke homogeneity."""
    from launcher.apps.overplot import Overplot

    _seed_folder(tmp_path, "db_fixture.txt")
    (tmp_path / "headerless.dat").write_text("0.01 0.5\n0.02 0.4\n")
    w = Overplot()
    w.folder_edit.setText(str(tmp_path))
    w.populate_file_list(str(tmp_path))
    _select_files(w, "db_fixture.txt", "headerless.dat")
    w.plot_selected()
    labels = [line.get_label() for line in w.figure.axes[0].lines]
    labels += [container.get_label() for container in w.figure.axes[0].containers]
    joined = " ".join(labels)
    assert "[direct_beam]" in joined
    assert "[unknown]" in joined
