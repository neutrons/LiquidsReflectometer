"""Model-level tests for the settings editor (T2).

Deliberately Qt-free and fast: `SettingsDocument` and `FIELD_SPEC` import no
Qt, so the whole editing model — load, validate, normalize, add/remove angle —
is exercised here in milliseconds without a display. The view's tests
(`launcher/tests/test_settings_editor.py`) then only have to cover wiring.
"""

import ast
import json
import os
import pathlib
import re
import subprocess
import sys
import textwrap

import pytest

from lr_reduction import field_spec as fs
from lr_reduction.new_reduction_from_file import json_to_config, save_config_json
from lr_reduction.nr_reduction_config import NRReductionConfig
from lr_reduction.save_reduced_data import make_json_safe
from lr_reduction.settings_document import SettingsDocument

# --------------------------------------------------------------------------
# FIELD_SPEC must mirror the real class, in both directions
# --------------------------------------------------------------------------


def test_field_spec_names_are_real_config_attributes():
    """Every FIELD_SPEC name must be a real ``__dict__`` key of the config.

    Stricter than the rationale the first version of this docstring gave. It
    described `json_to_config` raising `AttributeError`, but that gates on
    `hasattr`, which also passes for properties — including `base_path`, which
    has no setter and raises on assignment. The subset check the body actually
    performs is against `__dict__`, which excludes properties entirely, and that
    is the property worth pinning.
    """
    attributes = set(NRReductionConfig().__dict__)
    assert {f.name for f in fs.FIELD_SPEC} <= attributes


def test_field_spec_covers_every_config_attribute():
    """The other direction: a field nobody tabled is a field the editor hides."""
    attributes = set(NRReductionConfig().__dict__)
    assert attributes <= {f.name for f in fs.FIELD_SPEC}


def test_field_spec_defaults_match_the_config():
    """Defaults are copied from the class; drift here misinforms every prompt."""
    config = NRReductionConfig()
    mismatched = {
        f.name: (f.default, getattr(config, f.name))
        for f in fs.FIELD_SPEC
        if f.default != getattr(config, f.name)
    }
    assert not mismatched


def test_field_spec_excludes_base_path():
    """`base_path` passes `hasattr` and raises on `setattr` — it must not be tabled.

    This is the trap the name-coverage guard alone would not catch: it is a
    property with no setter, so a "mirror every attribute" loop that used
    `dir()` instead of `__dict__` would include it and then fail at load.
    """
    assert "base_path" not in fs.BY_NAME
    with pytest.raises(AttributeError):
        setattr(NRReductionConfig(), "base_path", "/tmp")


def test_normalized_document_loads_back_through_json_to_config():
    """The end-to-end contract: what we save, the reducer can read."""
    doc = SettingsDocument()
    doc.add_angle(DBname="db.dat", RB_Ymin=100, RB_Ymax=150)
    reloaded = json_to_config(doc.normalize())
    assert isinstance(reloaded, NRReductionConfig)


def test_per_angle_names_match_the_config_shape():
    """The 13 per-angle fields, pinned against the class rather than a comment."""
    config = NRReductionConfig()
    list_valued = {k for k, v in config.__dict__.items() if isinstance(v, list)}
    # data_x_range is a two-element detector range, not one entry per angle.
    assert set(fs.PER_ANGLE_NAMES) - set(fs.OPTIONAL_LIST_NAMES) == list_valued - {"data_x_range"}
    assert len(fs.PER_ANGLE_NAMES) == 13


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------


def test_defaults_match_a_fresh_config():
    assert SettingsDocument().to_dict() == NRReductionConfig().__dict__


def test_json_seed_round_trips(tmp_path):
    seed = tmp_path / "settings.json"
    doc = SettingsDocument()
    doc.set("Sname", "my_reduction")
    doc.set("qmax", 0.25)
    doc.save(seed)

    assert SettingsDocument.from_file(seed).to_dict() == doc.to_dict()
    assert json.loads(seed.read_text())["Sname"] == "my_reduction"


def test_dat_header_seed(tmp_path):
    """A pre-reduced .dat carries its config in a `# Config:` header line."""
    dat = tmp_path / "reduced.dat"
    config = {"Sname": "from_dat", "qmax": 0.42}
    dat.write_text(
        "# y\n# Config: %s\n# columns = Q, R, dR, dQ\n0.01 1.0 0.1 0.001\n" % json.dumps(config)
    )
    doc = SettingsDocument.from_file(dat)
    assert doc.get("Sname") == "from_dat"
    assert doc.get("qmax") == 0.42


def test_unknown_key_in_a_seed_is_reported_by_name(tmp_path):
    seed = tmp_path / "bad.json"
    seed.write_text(json.dumps({"Snam": "typo"}))
    with pytest.raises(ValueError, match="Snam"):
        SettingsDocument.from_file(seed)


# --------------------------------------------------------------------------
# Angles
# --------------------------------------------------------------------------


def test_add_angle_grows_every_per_angle_field():
    """Every list an angle is read from moves together, or the arrays silently desynchronise.

    Growing only the obvious eight was the original T2 defect: `RBnum` is easy to
    forget, and a short array shifts every subsequent angle's settings by one.

    Since editor-angle-count (G6), a list the reducer fills, broadcasts or derives is
    left compact when no value is given: empty, one broadcast entry, or None. The
    reducer expands it to every angle itself (nr_reduction_calc.py:42-43, :77-79,
    :99-110), and the unset entries this test used to require were what made an
    editor-authored file unreducible (F5). Every angle-defining list still grows.
    """
    doc = SettingsDocument()
    doc.add_angle()
    doc.add_angle()
    for name in fs.PER_ANGLE_NAMES:
        if name in fs.ANGLE_DEFINING_NAMES:
            assert len(doc.get(name)) == 2, name
        else:
            assert doc.get(name) in ([], None), name
    assert doc.n_angles == 2


def test_add_angle_leaves_an_unset_optional_list_unset():
    """`LambdaMin=None` means "derive from the choppers" — a real state, not a gap.

    `nr_reduction_calc.py:381-383` derives the bound when it is None, and
    `:70-71` raises if a supplied list is shorter than the run count. So
    materialising `[None, None]` on the first add would convert a valid
    "derive it" config into one that reports a length but carries no values,
    which `web_report.py:547` then indexes.
    """
    doc = SettingsDocument()
    doc.add_angle()
    assert doc.get("LambdaMin") is None
    assert doc.get("LambdaMax") is None


def test_add_angle_materializes_an_optional_list_when_given_a_value():
    doc = SettingsDocument()
    doc.add_angle()
    doc.add_angle(LambdaMin=2.5)
    assert doc.get("LambdaMin") == [None, 2.5]
    assert len(doc.get("LambdaMin")) == doc.n_angles


def test_remove_angle_shrinks_every_per_angle_field():
    doc = SettingsDocument()
    doc.add_angle(DBname="a.dat")
    doc.add_angle(DBname="b.dat")
    doc.add_angle(DBname="c.dat")
    doc.remove_angle(1)
    assert doc.get("DBname") == ["a.dat", "c.dat"]
    assert doc.n_angles == 2


def test_set_angle_field_uses_the_index_it_is_given():
    """The active-row trap, pinned at the model layer.

    The view must pass the row being edited, never the row that happens to be
    selected. Enforced here by construction: the model has no notion of a
    current row to fall back on.
    """
    doc = SettingsDocument()
    for name in ("a.dat", "b.dat", "c.dat"):
        doc.add_angle(DBname=name)
    # Index 1 of 3, deliberately. Editing index 0 is satisfied by an
    # implementation that ignores the index entirely and hard-codes 0 — which
    # is what the first version of this test asserted, and it passed against
    # exactly that mutation.
    doc.set_angle_field(1, "DBname", "edited.dat")
    assert doc.get("DBname") == ["a.dat", "edited.dat", "c.dat"]


def test_set_angle_field_rejects_an_out_of_range_index():
    doc = SettingsDocument()
    doc.add_angle()
    with pytest.raises(IndexError):
        doc.set_angle_field(5, "DBname", "x.dat")


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------


def test_validate_accepts_a_fresh_document():
    assert SettingsDocument().validate() == []


def test_validate_rejects_an_unknown_method():
    doc = SettingsDocument()
    doc.add_angle(method_per_run="constantBanana")
    assert any("constantBanana" in m for m in doc.validate())


def test_validate_accepts_methods_case_insensitively():
    """`nr_reduction_calc.py:81` lowercases before checking `:84`."""
    doc = SettingsDocument()
    doc.add_angle(method_per_run="CONSTANTQ")
    assert doc.validate() == []


def test_validate_accepts_a_broadcast_method_per_run():
    """A single entry is broadcast to every angle (`nr_reduction_calc.py:76-78`).

    A strict equal-length rule over all per-angle fields would reject this
    valid configuration.
    """
    doc = SettingsDocument()
    doc.add_angle()
    doc.add_angle()
    doc.set("method_per_run", ["meanTheta"])
    assert doc.validate() == []


def test_validate_flags_a_short_non_broadcast_array():
    doc = SettingsDocument()
    doc.add_angle()
    doc.add_angle()
    doc.set("DBname", ["only_one.dat"])
    assert any("DBname" in m for m in doc.validate())


def test_validate_flags_a_partially_specified_optional_list():
    """Detection is complete even though the fix is the scientist's."""
    doc = SettingsDocument()
    doc.add_angle()
    doc.add_angle(LambdaMin=2.5)
    messages = doc.validate()
    assert any("LambdaMin" in m for m in messages)


def test_validate_flags_a_value_outside_its_range():
    doc = SettingsDocument()
    doc.set("Qline_threshold", 4.0)
    assert any("Qline_threshold" in m for m in doc.validate())


def test_validate_rejects_an_unknown_enumerated_value():
    doc = SettingsDocument()
    doc.set("peak_type", "triangle")
    assert any("peak_type" in m for m in doc.validate())


def test_set_rejects_an_unknown_field_name():
    with pytest.raises(KeyError):
        SettingsDocument().set("Snam", "typo")


# --------------------------------------------------------------------------
# Normalization and diffing
# --------------------------------------------------------------------------


def test_normalize_drops_runtime_owned_fields():
    """Literal names, not the tuple that drives the filter.

    Iterating RUNTIME_OWNED_NAMES here asserted that the filter agrees with
    itself: emptying the tuple, or dropping runtime_owned from RBnum, left this
    green while normalize() happily emitted an authored RBnum into the
    reduction — the precise thing its docstring forbids.
    """
    doc = SettingsDocument()
    doc.add_angle(RBnum=197912)
    normalized = doc.normalize()
    assert {"RBnum", "LambdaMinUse", "LambdaMaxUse"}.isdisjoint(normalized)
    assert "Sname" in normalized


def test_runtime_owned_names_are_exactly_the_three_expected():
    """Pin the tuple itself, separately from the behaviour it drives."""
    assert set(fs.RUNTIME_OWNED_NAMES) == {"RBnum", "LambdaMinUse", "LambdaMaxUse"}


def test_changed_vs_seed_reports_only_what_moved():
    doc = SettingsDocument()
    assert doc.changed_vs_seed() == {}
    doc.set("Sname", "changed")
    assert doc.changed_vs_seed() == {"Sname": ("reduction_output", "changed")}


def test_changed_vs_seed_is_relative_to_the_loaded_seed(tmp_path):
    seed = tmp_path / "seed.json"
    seed.write_text(json.dumps({"Sname": "seeded"}))
    doc = SettingsDocument.from_file(seed)
    assert doc.changed_vs_seed() == {}
    doc.set("Sname", "edited")
    assert doc.changed_vs_seed() == {"Sname": ("seeded", "edited")}


# --------------------------------------------------------------------------
# The seam
# --------------------------------------------------------------------------


@pytest.mark.parametrize("module", ["settings_document", "field_spec"])
def test_model_modules_import_no_qt(module):
    """T3 builds on this seam; a stray Qt import would cost it the fast tests.

    Parsed rather than grepped. A substring search reports the word "qtpy"
    wherever it appears — including in the docstring that explains the module
    is Qt-free — so the first version of this guard failed on prose. `ast` sees
    only actual import statements.
    """
    source = (pathlib.Path(__file__).parents[3] / "src" / "lr_reduction" / f"{module}.py").read_text()
    imported = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.Import):
            imported.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module)
    offenders = {name for name in imported if name.split(".")[0] in {"qtpy", "PyQt5", "PyQt6", "PySide2", "PySide6"}}
    assert not offenders


# --------------------------------------------------------------------------
# C1 — a short per-angle column must not raise out of a Qt slot
# --------------------------------------------------------------------------


def test_editing_a_short_per_angle_column_pads_instead_of_raising():
    """The abort path: IndexError in a Qt slot reaches qFatal() and kills the app.

    A short column is not exotic. The reducer sanctions a length-1
    `method_per_run` broadcast to every angle, and `normalize()` drops the
    runtime-owned `RBnum`, so the editor's own save/reload round trip produces
    one.
    """
    doc = SettingsDocument()
    doc.add_angle()
    doc.add_angle()
    doc.set("DBname", ["only_one.dat"])
    doc.set_angle_field(1, "DBname", "second.dat")
    assert doc.get("DBname") == ["only_one.dat", "second.dat"]


def test_a_string_where_a_per_angle_list_belongs_is_not_exploded():
    """`len()` on a string succeeds, so a len-based guard ran `list("abc")`.

    The paired short-column test uses a real list, so it exercises the padding
    branch but never the wrong-type one — it stayed green with the isinstance
    guard reverted.
    """
    doc = SettingsDocument.from_dict({"DBname": "abc", "tof_min": [1.0, 2.0]})
    assert doc.n_angles == 2
    doc.set_angle_field(1, "DBname", "x.dat")
    assert doc.get("DBname") == [None, "x.dat"]


def test_the_editors_own_round_trip_is_editable():
    """Save, reload, edit — the crash cycle, using only the editor's artifacts."""
    doc = SettingsDocument()
    doc.add_angle(RBnum=197912)
    doc.add_angle(RBnum=197913)
    reloaded = SettingsDocument.from_dict(doc.normalize())
    assert reloaded.get("RBnum") == []
    reloaded.set_angle_field(1, "RBnum", 197913)
    assert reloaded.get("RBnum")[1] == 197913


@pytest.mark.parametrize(
    "payload",
    [
        pytest.param({"tof_min": 5}, id="int-where-a-list-belongs"),
        pytest.param({"LambdaMin": 3.5}, id="float-where-a-list-belongs"),
        pytest.param({"DBname": "one.dat"}, id="str-where-a-list-belongs"),
    ],
)
def test_a_non_sequence_per_angle_value_is_reported_not_iterated(payload):
    """validate() used to len()/enumerate() whatever the file contained."""
    doc = SettingsDocument.from_dict(payload)
    problems = doc.validate()
    assert any(next(iter(payload)) in message for message in problems)


# --------------------------------------------------------------------------
# C2 — values must arrive as their declared type
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name, text, expected",
    [
        pytest.param("RB_Ymin", "150", 150, id="list[int]-cell"),
        pytest.param("tof_min", "150", 150.0, id="list[float]-cell"),
        pytest.param("useBS", "False", False, id="list[bool]-cell-false"),
        pytest.param("useBS", "1", True, id="list[bool]-cell-true"),
        pytest.param("ScaleFactor", "1.05", 1.05, id="list[float]-scale"),
    ],
)
def test_a_per_angle_cell_is_coerced_to_its_element_type(name, text, expected):
    """`useBS` holding the string "False" is truthy: background gets subtracted
    when the scientist switched it off (`nr_reduction_calc`, `if useBS[i]`)."""
    value = fs.get(name).coerce_element(text)
    assert value == expected
    assert type(value) is type(expected)


def test_a_boolean_is_never_parsed_with_bool():
    """bool("False") is True — the whole point."""
    assert fs.get("useBS").coerce_element("False") is False
    assert fs.get("useBS").coerce_element("false") is False
    assert fs.get("useBS").coerce_element("0") is False


def test_a_two_value_scalar_list_parses_both_values():
    """`data_x_range` reaching the writer as "60, 210" produced x_min_pixel=6."""
    assert fs.get("data_x_range").coerce("60, 210") == [60, 210]
    assert fs.get("data_x_range").coerce("60 210") == [60, 210]


def test_a_value_whose_type_contradicts_the_field_is_reported():
    doc = SettingsDocument()
    doc.add_angle()
    doc.set("RB_Ymin", ["150"])
    assert any("RB_Ymin" in m for m in doc.validate())


def test_bounds_apply_to_per_angle_entries_too():
    """An out-of-range value, which the earlier version of this test lacked.

    It set ScaleFactor=1.0 — in range — and asserted validate() == [], so
    neutering the whole per-angle check loop left it green. It was the named
    guard for "the bounds finally apply to per-angle entries", and it never
    tested a bound.
    """
    doc = SettingsDocument()
    doc.add_angle(tof_min=-5.0)
    assert any("tof_min" in m and "below" in m for m in doc.validate())


def test_a_scalar_below_its_minimum_is_reported():
    doc = SettingsDocument()
    doc.set("dead_time", -1.0)
    assert any("dead_time" in m and "below" in m for m in doc.validate())


# --------------------------------------------------------------------------
# C4 — saving must not destroy the previous good file
# --------------------------------------------------------------------------


def test_save_leaves_the_previous_file_intact_when_the_write_fails(tmp_path, monkeypatch):
    """open(path, "w") truncates before a byte is produced."""
    target = tmp_path / "settings.json"
    target.write_text('{"good": "settings"}')

    def explode(*_args, **_kwargs):
        raise OSError(28, "No space left on device")

    monkeypatch.setattr(json, "dumps", explode)
    doc = SettingsDocument()
    with pytest.raises(OSError):
        doc.save(target)
    assert json.loads(target.read_text()) == {"good": "settings"}


def test_save_leaves_no_temporary_file_behind(tmp_path, monkeypatch):
    """Inject at os.replace, where a temp file exists to be cleaned up.

    The earlier version patched json.dumps, which runs BEFORE mkstemp — so no
    temp file was ever created and the assertion was true by construction.
    Deleting the whole cleanup block left it green.
    """
    target = tmp_path / "settings.json"

    def explode(*_args, **_kwargs):
        raise OSError(28, "No space left on device")

    monkeypatch.setattr(os, "replace", explode)
    with pytest.raises(OSError):
        SettingsDocument().save(target)
    assert list(tmp_path.iterdir()) == []


def test_save_refuses_to_write_through_a_symlink(tmp_path):
    """The overwrite dialog names the link, not the file that would be destroyed."""
    victim = tmp_path / "victim.txt"
    victim.write_text("precious")
    link = tmp_path / "settings.json"
    link.symlink_to(victim)
    with pytest.raises(ValueError, match="symbolic link"):
        SettingsDocument().save(link)
    assert victim.read_text() == "precious"


def test_saved_file_is_group_readable(tmp_path):
    """A shared IPTS directory: 0644 deliberately, not 0600."""
    target = tmp_path / "settings.json"
    SettingsDocument().save(target)
    assert target.stat().st_mode & 0o777 == 0o644


# --------------------------------------------------------------------------
# C5 — useCalcTheta is an enum, not a checkbox
# --------------------------------------------------------------------------


def test_theta_source_is_declared_as_a_choice_not_a_boolean():
    """The model-level half of C5: the DECLARATION, which is what the view reads.

    An earlier version of this test round-tripped 'sample_angle' through the
    document and passed even with the field declared `bool` — because the
    document stores whatever it is given and never consulted the type. The bug
    was always in the view, which built a checkbox from the declaration, so the
    declaration is the only part of it the model can actually pin. The widget
    itself is pinned in launcher/tests/test_settings_editor.py.
    """
    field = fs.get("useCalcTheta")
    assert field.type == "str"
    assert field.allowed == ("detector_angle", "sample_angle")


def test_sample_angle_survives_a_round_trip(tmp_path):
    seed = tmp_path / "s.json"
    seed.write_text(json.dumps({"useCalcTheta": "sample_angle"}))
    doc = SettingsDocument.from_file(seed)
    assert doc.get("useCalcTheta") == "sample_angle"
    assert doc.validate() == []
    out = tmp_path / "out.json"
    doc.save(out)
    assert json.loads(out.read_text())["useCalcTheta"] == "sample_angle"


def test_the_legacy_true_is_migrated_to_the_reducers_meaning():
    """The reducer maps True -> detector_angle; reporting it would cry wolf."""
    doc = SettingsDocument.from_dict({"useCalcTheta": True})
    assert doc.get("useCalcTheta") == "detector_angle"
    assert doc.validate() == []


def test_the_choice_lists_are_the_reducers_own():
    """Imported from one definition, not hand-mirrored."""
    from lr_reduction import reduction_domains

    assert fs.get("useCalcTheta").allowed is reduction_domains.CALC_THETA_CHOICES
    assert fs.get("method_per_run").allowed is reduction_domains.METHOD_CHOICES
    assert fs.get("peak_type").allowed is reduction_domains.PEAK_TYPE_CHOICES
    assert fs.get("DetResFn").allowed is reduction_domains.DET_RES_CHOICES


def test_the_reducer_validates_against_the_shared_domains():
    """Kept as a cheap structural pin; the BEHAVIOURAL guards are below.

    On its own this is a source grep — a hand-copy in double quotes alongside a
    dead reference to the constant would satisfy it. It stays because it names
    the intent at the point of the seam, but the tests that actually drive
    `_validate_config`, `fit_peak` and the domains are what make drift fail.
    """
    import inspect

    from lr_reduction import nr_reduction_calc

    source = inspect.getsource(nr_reduction_calc.NR_Reduction._validate_config)
    assert "domains.METHOD_CHOICES" in source
    assert "domains.CALC_THETA_CHOICES" in source
    assert "'meantheta'" not in source


# --------------------------------------------------------------------------
# C6 — the interface T3 consumes
# --------------------------------------------------------------------------


def test_overrides_reports_only_what_this_layer_contributes():
    doc = SettingsDocument()
    assert doc.overrides() == {}
    doc.set("Sname", "mine")
    assert doc.overrides() == {"Sname": "mine"}


def test_field_default_is_never_the_shared_object():
    """frozen=True freezes the binding, not the list behind it."""
    field = fs.get("method_per_run")
    borrowed = field.default_value()
    borrowed.append("meanTheta")
    assert field.default_value() == []


def test_every_declared_type_is_in_the_vocabulary():
    """A typo like "flaot" would silently yield an uncoerced text box."""
    assert {f.type for f in fs.FIELD_SPEC} <= set(fs.TYPES)


def test_field_names_are_unique():
    assert len(fs.BY_NAME) == len(fs.FIELD_SPEC)


# --------------------------------------------------------------------------
# Validation must not cry wolf
# --------------------------------------------------------------------------


def test_a_valid_three_angle_document_reports_no_problems():
    """Six of the nine problems v1 reported on a valid file were false.

    `RBnum` is runtime-owned and the five `default_if_empty` arrays are filled
    in by the reducer when empty. A panel that cries wolf on a good file teaches
    the scientist to ignore it, which is how a real problem goes unread.
    """
    # Loaded, NOT built with add_angle: the auto-defaulted arrays have to be
    # genuinely ABSENT for the exemption to be exercised. Growing them with
    # add_angle gives every column length 3, so the earlier version of this test
    # passed with the exemption disabled — it never reached it.
    doc = SettingsDocument.from_dict(
        {
            "DBname": ["db_0.dat", "db_1.dat", "db_2.dat"],
            "RB_Ymin": [100, 100, 100],
            "RB_Ymax": [150, 150, 150],
            "method_per_run": ["meanTheta"],
            # Supplied because it genuinely IS required per angle — the reducer
            # does not auto-default it and web_report indexes it directly
            # (BkgROI[idx]). validate() reporting an empty one for 3 angles is a
            # TRUE positive, so this test supplies it rather than adding an
            # exemption to make itself pass.
            "BkgROI": [[10, 20, 30, 40]] * 3,
        }
    )
    assert doc.n_angles == 3
    # The five the reducer fills in when empty, plus runtime-owned RBnum, are
    # the ones that must stay silent.
    assert doc.get("useBS") == []
    assert doc.get("ThetaShift") == []
    assert doc.get("RBnum") == []
    assert doc.validate() == []


def test_an_emptied_optional_list_returns_to_derive_from_choppers():
    """Otherwise touching one Lambda cell is a one-way door."""
    doc = SettingsDocument()
    doc.add_angle()
    doc.set_angle_field(0, "LambdaMin", 2.5)
    assert doc.get("LambdaMin") == [2.5]
    doc.set_angle_field(0, "LambdaMin", None)
    assert doc.get("LambdaMin") is None


@pytest.mark.parametrize(
    "name, value",
    [
        pytest.param("experiment_id", "/etc/passwd", id="absolute-path"),
        pytest.param("_Spath_override", "../../elsewhere", id="parent-traversal"),
        pytest.param("Sname", "../../../.bashrc", id="separator-in-a-file-name"),
    ],
)
def test_a_path_that_escapes_the_experiment_directory_is_reported(name, value):
    """These are joined verbatim into the output location by save_reduced_data."""
    doc = SettingsDocument()
    doc.set(name, value)
    assert any(name in message for message in doc.validate())


def test_the_model_pulls_no_qt_into_a_fresh_interpreter():
    """The property T3 depends on, asserted at the level it actually holds.

    The AST guard above is correct for what it names — this file's own import
    statements — but it cannot see a TRANSITIVE one. Adding
    `from launcher.app_identity import ensure_identity` to settings_document
    keeps every AST check passing while genuinely loading PyQt5.QtCore.

    A subprocess is the only honest test: this session's other tests import Qt,
    so `sys.modules` in-process proves nothing.
    """
    program = textwrap.dedent(
        """
        import sys
        from lr_reduction import field_spec, settings_document
        qt = sorted(m for m in sys.modules
                    if m.split(".")[0] in {"qtpy", "PyQt5", "PyQt6", "PySide2", "PySide6"})
        assert not qt, f"model imports pulled in Qt: {qt}"
        doc = settings_document.SettingsDocument()
        doc.add_angle(DBname="a.dat")
        assert doc.n_angles == 1
        assert doc.validate() == []
        print("clean")
        """
    )
    result = subprocess.run(
        [sys.executable, "-c", program],
        capture_output=True, text=True, timeout=300,
    )
    assert result.returncode == 0, result.stderr
    assert "clean" in result.stdout


# --------------------------------------------------------------------------
# C3 — the handoff to the reducer
# --------------------------------------------------------------------------


def test_config_is_the_object_the_reduction_receives():
    """The seam T3 builds on, named by the v1 review and still unpinned in v2.

    Renaming the property left the whole suite green, because nothing — tests
    included — consumed it.
    """
    doc = SettingsDocument()
    config = doc.config
    assert isinstance(config, NRReductionConfig)
    doc.set("Sname", "handed_off")
    assert config.Sname == "handed_off"
    doc.add_angle(DBname="db.dat")
    assert config.DBname == ["db.dat"]
    assert doc.config is config


# --------------------------------------------------------------------------
# The domains, driven through the code that enforces them
# --------------------------------------------------------------------------


def _bare_reduction(**config_values):
    """An NR_Reduction with a config, skipping __init__'s heavy setup."""
    from lr_reduction import nr_reduction_calc

    instance = nr_reduction_calc.NR_Reduction.__new__(nr_reduction_calc.NR_Reduction)
    instance.config = NRReductionConfig()
    for key, value in config_values.items():
        setattr(instance.config, key, value)
    return instance


@pytest.mark.parametrize("method", fs.METHOD_CHOICES)
def test_every_declared_method_is_accepted_by_the_reducer(method):
    """Behavioural, not a source grep.

    The previous guard used inspect.getsource and asserted a substring, which a
    hand-copy in double quotes would satisfy. This drives the actual validator.
    """
    reduction = _bare_reduction(
        RBnum=[1], DBname=["db.dat"], method_per_run=[method],
        RB_Ymin=[1], RB_Ymax=[2],
    )
    reduction._validate_config()
    assert reduction.config.method_per_run == [method.lower()]


def test_a_method_outside_the_domain_is_rejected_by_the_reducer():
    reduction = _bare_reduction(
        RBnum=[1], DBname=["db.dat"], method_per_run=["constantBanana"],
        RB_Ymin=[1], RB_Ymax=[2],
    )
    with pytest.raises(ValueError, match="Invalid method"):
        reduction._validate_config()


@pytest.mark.parametrize("choice", fs.CALC_THETA_CHOICES)
def test_every_declared_theta_source_is_accepted_by_the_reducer(choice):
    reduction = _bare_reduction(
        RBnum=[1], DBname=["db.dat"], RB_Ymin=[1], RB_Ymax=[2],
        useCalcTheta=choice,
    )
    reduction._validate_config()
    assert reduction.config.useCalcTheta == choice


@pytest.mark.parametrize("peak_type", fs.PEAK_TYPE_CHOICES)
def test_every_declared_peak_type_is_accepted_by_the_fitter(peak_type):
    """fit_peak raises for an unknown peaktype; a declared one must get past it."""
    import numpy as np

    from lr_reduction import nr_tools

    ypix = np.arange(40.0)
    iY = 100.0 * np.exp(-0.5 * ((ypix - 20.0) / 3.0) ** 2) + 1.0
    nr_tools.fit_peak(ypix, iY, peaktype=peak_type, bkgtype="none")


def test_a_peak_type_outside_the_domain_is_rejected_by_the_fitter(monkeypatch):
    """Rejection AND derivation, in one test.

    Matching "peaktype must be" alone proves only that it raises — a hardcoded
    message satisfies it, so it could not tell a derived error from a copied
    one. Extending the domain and asserting the message follows is what
    actually pins the derivation.
    """
    import numpy as np

    from lr_reduction import nr_tools, reduction_domains

    ypix = np.arange(40.0)
    iY = 100.0 * np.exp(-0.5 * ((ypix - 20.0) / 3.0) ** 2) + 1.0

    monkeypatch.setattr(
        reduction_domains, "PEAK_TYPE_CHOICES", ("gauss", "supergauss", "hexagauss")
    )
    with pytest.raises(ValueError) as raised:
        nr_tools.fit_peak(ypix, iY, peaktype="sombrero", bkgtype="none")
    assert "hexagauss" in str(raised.value)


def test_a_detector_resolution_outside_the_domain_is_rejected_by_nr_tools(monkeypatch):
    """The same derivation guard for the other nr_tools domain."""
    import numpy as np

    from lr_reduction import nr_tools, reduction_domains

    monkeypatch.setattr(reduction_domains, "DET_RES_CHOICES", ("rectangular", "hexbox"))
    with pytest.raises(ValueError) as raised:
        nr_tools.calc_beam_on_detector(
            Ypix=np.arange(64.0), CenPix=32.0, Si=0.25, S1=0.39, dS1Si=1000.0,
            dSiSam=100.0, dSamDet=1500.0, mmpix=0.7, DetRes=0.8,
            DetResFn="sombrero",
        )
    assert "hexbox" in str(raised.value)


def test_the_detector_resolution_disagreement_is_recorded_not_hidden():
    """'none' is accepted by one consumer and crashes the other.

    reduction_domains records that rather than picking a side, and the editor
    reports the reason instead of a bare "not one of".
    """
    from lr_reduction import reduction_domains

    assert "none" not in reduction_domains.DET_RES_CHOICES
    assert "none" in reduction_domains.DET_RES_TOLERATED
    message = fs.get("DetResFn").check_element("none")
    assert "UnboundLocalError" in message


def test_default_if_empty_names_are_exactly_the_reducers_optional_arrays():
    """Pin the tuple, and give the export a consumer.

    These are the five the reducer fills in under "Set defaults for optional
    arrays". BkgROI is deliberately NOT among them — it is indexed per angle by
    web_report and has no auto-default, so an empty one is a real problem.
    """
    assert set(fs.DEFAULT_IF_EMPTY_NAMES) == {
        "ThetaShift", "useBS", "ScaleFactor", "tof_min", "tof_max",
    }
    assert "BkgROI" not in fs.DEFAULT_IF_EMPTY_NAMES


# --------------------------------------------------------------------------
# The domain contract: what derives, and what only a pin can catch
# --------------------------------------------------------------------------
#
# The validators derive from reduction_domains, so they cannot drift. The
# DISPATCH cannot derive — each branch computes something different — so adding
# a value to a domain does NOT teach the reducer to compute it. These pins are
# what catch that, and the positive drivers below assert the other direction:
# every value the editor offers is one a consumer actually accepts.


def test_method_choices_are_exactly_what_the_reducer_handles():
    assert fs.METHOD_CHOICES == ("meanTheta", "constantQ", "constantTOF")


def test_peak_type_choices_are_exactly_what_the_fitter_dispatches_on():
    assert fs.PEAK_TYPE_CHOICES == ("gauss", "supergauss")


def test_det_res_choices_are_exactly_what_the_consumers_dispatch_on():
    assert fs.DET_RES_CHOICES == ("rectangular", "gaussian")


def test_theta_dispatch_is_a_subset_of_the_declared_methods():
    """`constantTOF` is validated but has no theta branch — it routes elsewhere.

    Recorded as its own tuple so the dispatch's error message names what that
    function can actually compute, instead of the full domain.
    """
    from lr_reduction import reduction_domains

    assert set(reduction_domains.THETA_DISPATCH_CHOICES) < set(fs.METHOD_CHOICES)
    assert reduction_domains.THETA_DISPATCH_CHOICES == ("constantQ", "meanTheta")


def _beam_kwargs(**overrides):
    import numpy as np

    kwargs = dict(
        Ypix=np.arange(64.0), CenPix=32.0, Si=0.25, S1=0.39, dS1Si=1000.0,
        dSiSam=100.0, dSamDet=1500.0, mmpix=0.7, DetRes=0.8,
    )
    kwargs.update(overrides)
    return kwargs


@pytest.mark.parametrize("resolution_function", fs.DET_RES_CHOICES)
def test_every_offered_detector_resolution_is_accepted(resolution_function):
    """The direction a contents pin cannot cover: the consumer must accept it.

    A pin says the domain equals a literal; only a call says the reducer can
    actually do something with each value.
    """
    from lr_reduction import nr_tools

    nr_tools.calc_beam_on_detector(**_beam_kwargs(DetResFn=resolution_function))


def test_the_tolerated_detector_resolution_is_accepted_by_nr_tools():
    """'none' is the recorded disagreement: nr_tools skips, the other consumer raises."""
    from lr_reduction import nr_tools

    nr_tools.calc_beam_on_detector(**_beam_kwargs(DetResFn="none"))
    nr_tools.calc_beam_on_detector(**_beam_kwargs(DetResFn=None))


def test_the_theta_dispatch_message_names_the_methods_it_handles(monkeypatch):
    """Derivation pin for the one message v3 left literal."""
    from lr_reduction import reduction_domains

    monkeypatch.setattr(
        reduction_domains, "THETA_DISPATCH_CHOICES", ("constantQ", "meanTheta", "hexTheta")
    )
    reduction = _bare_reduction()
    with pytest.raises(ValueError) as raised:
        reduction._calculate_theta_and_bins(ypix := None, 0.0, "sombrero")  # noqa: F841
    assert "hexTheta" in str(raised.value)


# --------------------------------------------------------------------------
# Pins for v3 repairs that had none (each was "suite green after mutation")
# --------------------------------------------------------------------------


def test_check_recurses_into_list_elements():
    """The container being a list said nothing about what was in it."""
    problem = fs.get("data_x_range").check(["[50", "200]"])
    assert "data_x_range" in problem
    assert "[50" in problem


def test_type_problem_recurses_into_nested_lists():
    problem = fs.get("BkgROI").check([["120", 130]])
    assert "BkgROI" in problem


@pytest.mark.parametrize(
    "name", ["subname", "DTCsubname", "BINsubname", "errBINsubname"]
)
def test_the_subname_siblings_reject_a_traversal(name):
    """All four share Sname's f-string join; only Sname was guarded in v2.

    Measured then: subname="/../../../../../../tmp/pwn" validated clean and the
    sink resolved to /tmp/pwn.dat.
    """
    doc = SettingsDocument()
    doc.set(name, "/../../tmp/pwn")
    assert any(name in message for message in doc.validate())


@pytest.mark.parametrize(
    "name",
    ["_Spath_override", "_NEXUSpathRB_override", "_DBpath_override", "_BINpath_override"],
)
def test_an_absolute_override_path_is_accepted(name):
    """These fields ARE the location, so absolute is their legitimate shape.

    It is what QFileDialog.getExistingDirectory returns. An earlier guard
    rejected exactly that while accepting a relative value, which resolves
    against the process working directory instead.
    """
    doc = SettingsDocument()
    doc.set(name, "/SNS/REF_L/IPTS-30101/shared/reduced")
    assert doc.validate() == []


def test_a_traversal_in_an_override_path_is_still_rejected():
    doc = SettingsDocument()
    doc.set("_Spath_override", "/SNS/REF_L/../../etc")
    assert any("_Spath_override" in message for message in doc.validate())


@pytest.mark.parametrize("name", ["DetResFn", "peak_type"])
def test_a_case_variant_is_normalised_on_entry(name):
    """nr_tools compares these exactly, so a stored 'Gaussian' never matches."""
    field = fs.get(name)
    variant = field.allowed[-1].upper()
    assert field.coerce(variant) == field.allowed[-1]


@pytest.mark.parametrize("name", ["DetResFn", "peak_type"])
def test_a_loaded_case_variant_is_reported(name):
    """from_dict does not coerce, so validation is the only thing that can catch it."""
    field = fs.get(name)
    doc = SettingsDocument.from_dict({name: field.allowed[-1].upper()})
    assert any("differs in case" in message for message in doc.validate())


def test_a_case_variant_of_a_case_insensitive_field_is_accepted():
    """method_per_run IS lower-cased by the reducer, so it must stay tolerant."""
    doc = SettingsDocument()
    doc.add_angle(method_per_run="CONSTANTQ")
    assert doc.validate() == []


def test_a_per_angle_enumerated_cell_is_normalised_on_entry():
    """The cell path, which the scalar test does not reach.

    `coerce` delegates to `coerce_element` for non-list fields, so both are
    pinned by one implementation — but a per-angle enumerated field only ever
    goes through `coerce_element`, and that path needs its own driver.
    """
    doc = SettingsDocument()
    doc.add_angle(method_per_run=fs.get("method_per_run").coerce_element("CONSTANTQ"))
    assert doc.get("method_per_run") == ["constantQ"]
    assert doc.validate() == []


# --------------------------------------------------------------------------
# editor-load-fidelity — a reducer-written file loads quietly, saves as it was
# written, and its runtime record is quiet
# --------------------------------------------------------------------------
#
# True == 1 and isinstance(True, int), so every assertion below about encoding
# or canonicalization compares type(), a repr, or the saved TEXT.
# `doc.get("useBS") == [True, True, False]` passes on [1, 1, 0] and guards
# nothing.


def _reducer_shaped_config(record=True):
    """A config shaped by the reduction's own writers, not by this editor.

    `useBS` holds integers because that is what the reducer writes:
    `nr_reduction_calc.py:103` fills an empty one with `[1] * n` and
    `new_reduction_from_template.py:182` writes 0 for "off". The runtime record
    is one scalar per call (`nr_reduction_calc.py:385-391`), not the
    `list[float]` FIELD_SPEC declares for it.
    """
    config = NRReductionConfig()
    config.experiment_id = "IPTS-00000"
    config.RBnum = [201282, 201283, 201284]
    config.DBname = ["db_a.dat", "db_b.dat", "db_c.dat"]
    config.method_per_run = ["meanTheta"]
    config.RB_Ymin = [140, 141, 142]
    config.RB_Ymax = [150, 151, 152]
    config.BkgROI = [[120, 130], [121, 131], [122, 132]]
    config.useBS = [1] * 3
    config.useBS[2] = 0
    config.useCalcTheta = "detector_angle"
    if record:
        config.LambdaMinUse = 2.95
        config.LambdaMaxUse = 6.1
    return config


def _reducer_written_json(directory):
    """The fixture as a settings file, written by the reduction's own saver."""
    path = pathlib.Path(directory) / "run_settings.json"
    save_config_json(path, _reducer_shaped_config())
    return path


def _reducer_written_dat(directory):
    """The same config as a reduced .dat carries it, on its `# Config:` line."""
    path = pathlib.Path(directory) / "run.dat"
    config_json = json.dumps(make_json_safe(_reducer_shaped_config().__dict__))
    path.write_text(f"# Config: {config_json}\n# columns = Q, R, dR, dQ\n0.01 1.0 0.1 0.001\n")
    return path


_REDUCER_WRITTEN = [
    pytest.param(_reducer_written_json, id="json"),
    pytest.param(_reducer_written_dat, id="dat-header"),
]


def _saved_block(name, value):
    """The exact text `save()` writes for one top-level key (indent=2)."""
    return json.dumps({name: value}, indent=2)[len("{\n  ") : -len("\n}")]


@pytest.mark.parametrize("write", _REDUCER_WRITTEN)
def test_a_reducer_written_file_validates_clean(tmp_path, write):
    """Three useBS integers and the scalar record each used to be a problem line."""
    assert SettingsDocument.from_file(write(tmp_path)).validate() == []


@pytest.mark.parametrize(
    "load",
    [
        pytest.param(lambda d: SettingsDocument.from_file(_reducer_written_json(d)), id="json"),
        pytest.param(lambda d: SettingsDocument.from_file(_reducer_written_dat(d)), id="dat-header"),
        pytest.param(lambda _d: SettingsDocument.from_dict({"useBS": [1, True, 0]}), id="mixed"),
    ],
)
def test_a_loaded_useBS_is_held_as_booleans(tmp_path, load):
    held = load(tmp_path).get("useBS")
    assert [type(v) for v in held] == [bool, bool, bool]
    assert repr(held) == "[True, True, False]"


def test_an_injected_config_with_integer_useBS_validates_clean():
    """No load and so no canonicalization: acceptance has to hold on its own."""
    doc = SettingsDocument(_reducer_shaped_config(record=False))
    assert [type(v) for v in doc.get("useBS")] == [int, int, int]
    assert doc.validate() == []


@pytest.mark.parametrize(
    "build, expected",
    [
        pytest.param(lambda: SettingsDocument.from_dict({"useBS": [1, 1, 0]}), [1, 1, 0], id="ints"),
        pytest.param(
            lambda: SettingsDocument.from_dict({"useBS": [True, True, False]}), [1, 1, 0], id="bools"
        ),
        pytest.param(lambda: SettingsDocument.from_dict({"useBS": [1, True, 0]}), [1, 1, 0], id="mixed"),
        pytest.param(lambda: SettingsDocument(_reducer_shaped_config()), [1, 1, 0], id="injected-ints"),
        pytest.param(lambda: SettingsDocument.from_dict({"useBS": [False, None]}), [0, None], id="unset-entry"),
    ],
)
def test_useBS_is_saved_as_ones_and_zeros(tmp_path, build, expected):
    """The reduction's own spelling, whatever the document was built from."""
    text = build().save(tmp_path / "out.json").read_text()
    saved = json.loads(text)["useBS"]
    assert [type(v) for v in saved] == [type(v) for v in expected]
    assert _saved_block("useBS", expected) in text


def test_normalize_writes_useBS_as_ones_and_zeros():
    normalized = SettingsDocument.from_dict({"useBS": [True, True, False]}).normalize()["useBS"]
    assert [type(v) for v in normalized] == [int, int, int]
    assert normalized == [1, 1, 0]


def test_a_scalar_true_is_saved_as_true(tmp_path):
    """Only the declared list is integer-encoded; a scalar boolean is written as held.

    (v2: the `1` leg that asserted a loaded scalar 1 saves as true is withdrawn — the reducer reads
    useGravity with `is True`, so that rewrite changed the reduction. The scalar matrix below replaces it.)
    """
    doc = SettingsDocument.from_dict({"Normalize": True})
    assert doc.validate() == []
    assert '"Normalize": true' in doc.save(tmp_path / "out.json").read_text()


@pytest.mark.parametrize(
    "entry",
    [
        pytest.param(2, id="two"),
        pytest.param("0", id="string-zero"),
        pytest.param(1.0, id="float-one"),
        pytest.param([1], id="nested"),
    ],
)
def test_a_non_boolean_useBS_entry_is_reported_and_kept(tmp_path, entry):
    """Never coerced on load: "0" is truthy to the reducer (`if useBS[i]:`), so
    accepting it silently would subtract a background its author switched off."""
    doc = SettingsDocument.from_dict({"useBS": [True, entry]})
    problems = [m for m in doc.validate() if "(useBS)" in m]
    assert len(problems) == 1
    assert "at angle 1" in problems[0]
    assert "true/false (or 1/0)" in problems[0]
    assert repr(doc.get("useBS")[1]) == repr(entry)
    saved = json.loads(doc.save(tmp_path / "out.json").read_text())["useBS"]
    assert repr(saved[1]) == repr(entry)


def test_a_non_list_useBS_is_reported_once_and_saved_as_held(tmp_path):
    doc = SettingsDocument.from_dict({"useBS": 1})
    assert len([m for m in doc.validate() if "(useBS)" in m]) == 1
    assert repr(doc.get("useBS")) == "1"
    assert repr(json.loads(doc.save(tmp_path / "out.json").read_text())["useBS"]) == "1"


@pytest.mark.parametrize("name", ["LambdaMinUse", "LambdaMaxUse"])
@pytest.mark.parametrize(
    "value",
    [
        pytest.param(2.95, id="scalar"),
        pytest.param(None, id="unset"),
        pytest.param([2.95, 3.1], id="list"),
        pytest.param("n/a", id="text"),
    ],
)
def test_the_runtime_record_is_never_a_problem(name, value):
    """Not an input: `nr_reduction_calc.py:385-391` overwrites it before its first
    use, and a problem line on a field nobody can edit has no remedy."""
    assert SettingsDocument.from_dict({name: value}).validate() == []


def test_the_round_trip_is_idempotent_and_keeps_useBS_as_written(tmp_path):
    source = _reducer_written_json(tmp_path)
    first = SettingsDocument.from_file(source).save(tmp_path / "first.json")
    second = SettingsDocument.from_file(first).save(tmp_path / "second.json")
    assert first.read_bytes() == second.read_bytes()
    # json.dumps tells 1 from true, so this compares the written spelling.
    assert json.dumps(json.loads(first.read_text())["useBS"]) == json.dumps(
        json.loads(source.read_text())["useBS"]
    )


def test_the_seed_is_canonical_too(tmp_path):
    """Compared by repr: [1, 1, 0] == [True, True, False], so equality cannot
    tell a seed taken before canonicalization from one taken after."""
    doc = SettingsDocument.from_file(_reducer_written_json(tmp_path))
    assert doc.changed_vs_seed() == {}
    doc.set_angle_field(0, "useBS", False)
    before, after = doc.changed_vs_seed()["useBS"]
    assert repr(before) == "[True, True, False]"
    assert repr(after) == "[False, True, False]"


def test_the_integer_encoded_names_are_exactly_useBS():
    """Pin the declaration separately from the behaviour it drives."""
    assert set(fs.INT_ENCODED_NAMES) == {"useBS"}


# --------------------------------------------------------------------------
# editor-load-fidelity v2 — load -> save never changes what the reduction does with a
# declared boolean (B8). The reducer reads useGravity with `is True`
# (nr_reduction_calc.py:1079), so a scalar is left exactly as loaded: 1 stays 1.
# --------------------------------------------------------------------------

#: Derived, not typed: the matrix below must grow when a boolean is added.
SCALAR_BOOLEANS = tuple(f.name for f in fs.FIELD_SPEC if f.type == "bool")
BOOLEAN_SPELLINGS = [
    pytest.param(1, id="one"),
    pytest.param(0, id="zero"),
    pytest.param(True, id="true"),
    pytest.param(False, id="false"),
]


def test_the_scalar_booleans_are_the_seven_whose_readers_were_checked():
    """A pin on the derivation the matrix iterates. Before updating it for an eighth boolean, read
    how the reduction reads that one: truthiness, identity (`is True`), or formatting into a header."""
    assert set(SCALAR_BOOLEANS) == {
        "Normalize", "AutoScale", "plotON", "plotQ4", "save8col", "useGravity", "use_emission_time",
    }


@pytest.mark.parametrize("value", BOOLEAN_SPELLINGS)
@pytest.mark.parametrize("name", SCALAR_BOOLEANS)
def test_load_then_save_keeps_a_scalar_boolean_exactly_as_written(tmp_path, name, value):
    source = tmp_path / "source.json"
    source.write_text(json.dumps({name: value}))
    saved = json.loads(SettingsDocument.from_file(source).save(tmp_path / "saved.json").read_text())[name]
    assert type(saved) is type(value)
    assert saved == value


@pytest.mark.parametrize("value", BOOLEAN_SPELLINGS)
@pytest.mark.parametrize("name", SCALAR_BOOLEANS)
def test_an_integer_scalar_boolean_is_reported_without_offering_1_or_0(name, value):
    lines = [m for m in SettingsDocument.from_dict({name: value}).validate() if f"({name})" in m]
    if type(value) is bool:
        assert lines == []
    else:
        assert len(lines) == 1
        assert "true/false" in lines[0]
        assert "1/0" not in lines[0]


@pytest.mark.parametrize("value", BOOLEAN_SPELLINGS)
def test_the_reduction_reads_a_saved_useGravity_as_it_read_the_source(tmp_path, value):
    """The rejection's reproduction (review 8b62952): `nr_reduction_calc.py:1079` tests
    `useGravity is True`, and json_to_config does no conversion, so 1 means gravity OFF."""
    source = tmp_path / "source.json"
    source.write_text(json.dumps({"useGravity": value}))
    saved = SettingsDocument.from_file(source).save(tmp_path / "saved.json")

    def gravity_on(path):
        return json_to_config(json.loads(path.read_text())).useGravity is True

    assert gravity_on(saved) == gravity_on(source)


@pytest.mark.parametrize("value", [pytest.param(1, id="one"), pytest.param(0, id="zero")])
@pytest.mark.parametrize("name", SCALAR_BOOLEANS)
def test_an_injected_integer_scalar_boolean_is_reported_and_kept(tmp_path, name, value):
    config = NRReductionConfig()
    setattr(config, name, value)
    doc = SettingsDocument(config)
    assert len([m for m in doc.validate() if f"({name})" in m]) == 1
    assert repr(doc.get(name)) == repr(value)
    assert repr(json.loads(doc.save(tmp_path / "out.json").read_text())[name]) == repr(value)


@pytest.mark.parametrize("value", BOOLEAN_SPELLINGS)
def test_a_saved_useBS_entry_is_read_as_the_source_entry(tmp_path, value):
    """Its readers: truthiness (nr_reduction_calc.py:509, :979) and `== 1`
    (new_reduction_from_template.py:224) — both must see what they saw in the source."""
    source = tmp_path / "source.json"
    source.write_text(json.dumps({"useBS": [value]}))
    saved = json.loads(SettingsDocument.from_file(source).save(tmp_path / "saved.json").read_text())["useBS"][0]
    assert type(saved) is int
    assert bool(saved) == bool(value)
    assert (saved == 1) == (value == 1)


# --------------------------------------------------------------------------
# editor-angle-count — the editor counts, adds and saves angles the way the
# reduction reads them (nr_reduction_calc.py:61-110)
# --------------------------------------------------------------------------

_THREE_ANGLES = {
    "RBnum": [201282, 201283, 201284],
    "DBname": ["db_a.dat", "db_b.dat", "db_c.dat"],
    "RB_Ymin": [140, 141, 142],
    "RB_Ymax": [150, 151, 152],
    "BkgROI": [[120, 130], [121, 131], [122, 132]],
}


def _surplus_document(extra=1):
    """Three angles by RBnum, with `useBS` longer — the real files' shape
    (IPTS-36119 reduce_settings.json holds useBS [1, 1, 1, 1] for three runs)."""
    return SettingsDocument.from_dict({**_THREE_ANGLES, "useBS": [1] * (3 + extra)})


def test_a_surplus_entry_is_not_a_problem():
    assert _surplus_document().validate() == []


@pytest.mark.parametrize("extra", [1, 5], ids=["one-extra", "five-extra"])
def test_a_surplus_entry_is_one_note_naming_the_field_and_the_count(extra):
    doc = _surplus_document(extra)
    notes = [line for line in doc.notes() if "(useBS)" in line]
    assert len(notes) == 1  # one line per field, never per entry
    assert f"{extra} extra" in notes[0]
    assert not any("(useBS)" in line for line in doc.validate())


def test_the_reduction_counts_its_angles_and_the_table_shows_every_row():
    doc = _surplus_document()
    assert doc.reduction_angles == 3
    assert doc.n_angles == 4
    # every angle-defining list empty: the count falls back to the longest list, so a
    # document holding only optional columns still shows them, without a note
    fallback = SettingsDocument.from_dict({"LambdaMin": [2.5, 2.6], "ThetaShift": [0.1, 0.2]})
    assert fallback.reduction_angles == 2 == fallback.n_angles
    assert fallback.notes() == []


def test_a_short_angle_defining_list_is_still_reported():
    """A pin against over-relaxing: the reducer raises on a short DBname (nr_reduction_calc.py:67-68)."""
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, "DBname": ["db_a.dat", "db_b.dat"]})
    assert "Direct-beam file (DBname) has 2 entries for 3 angles" in doc.validate()


def test_the_angle_defining_names_are_derived_and_every_per_angle_field_has_one_kind():
    assert set(fs.ANGLE_DEFINING_NAMES) == {"DBname", "RBnum", "RB_Ymin", "RB_Ymax", "BkgROI"}
    for name in fs.PER_ANGLE_NAMES:
        field = fs.get(name)
        kinds = [name in fs.ANGLE_DEFINING_NAMES, field.default_if_empty, field.broadcast_ok, field.optional_list]
        assert sum(kinds) == 1, name


def test_an_added_angle_has_one_index_in_every_list():
    """v2 (G6 revised): the new angle is the reduction's next one, index m, directly after
    the last real angle; the surplus entry shifts down one row and stays surplus. (v1 put
    it after the surplus rows, which turned them into angles with no required values.)"""
    doc = _surplus_document()
    doc.add_angle(DBname="d.dat", useBS=False)
    assert doc.reduction_angles == 4
    for name in ("DBname", "RBnum", "RB_Ymin", "RB_Ymax", "BkgROI"):
        assert len(doc.get(name)) == 4, name
    assert doc.get("DBname")[3] == "d.dat"
    assert repr(doc.get("useBS")) == "[True, True, True, False, True]"
    assert doc.n_angles == 5
    # the compact lists the reducer fills or derives are untouched
    assert doc.get("method_per_run") == []
    assert doc.get("ThetaShift") == []
    assert doc.get("LambdaMin") is None


def test_an_added_angle_leaves_a_broadcast_method_broadcast():
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, "method_per_run": ["meanTheta"]})
    doc.add_angle()
    assert doc.get("method_per_run") == ["meanTheta"]
    doc.add_angle(method_per_run="constantQ")
    assert doc.get("method_per_run") == ["meanTheta"] * 4 + ["constantQ"]


def test_setting_one_angle_of_a_broadcast_method_keeps_it_for_the_others():
    """The §5 row "set method_per_run for the new angle only": the earlier angles keep
    the broadcast value the reducer would have used (nr_reduction_calc.py:77-79)."""
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, "method_per_run": ["meanTheta"]})
    doc.set_angle_field(0, "method_per_run", "constantQ")
    assert doc.get("method_per_run") == ["constantQ", "meanTheta", "meanTheta"]


def test_an_added_angle_leaves_an_empty_default_list_empty():
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, "ThetaShift": []})
    doc.add_angle()
    assert doc.get("ThetaShift") == []


def test_the_reducer_accepts_what_the_editor_writes():
    """F5: a file authored in the editor validated clean and could not be reduced.

    Runs the reducer's own path: NR_Reduction.__init__ (which defaults an empty
    method_per_run, nr_reduction_calc.py:42-43) and _validate_config (:55-110). The
    module's _bare_reduction skips __init__, so it would not show that default.
    """
    from lr_reduction.nr_reduction_calc import NR_Reduction

    doc = SettingsDocument()
    for k in range(2):
        doc.add_angle(DBname=f"db_{k}.dat", RB_Ymin=140, RB_Ymax=150, BkgROI=[120, 130])
    config = json_to_config(doc.normalize())
    config.RBnum = [201282, 201283]  # as reduce_from_file sets it from the runs (new_reduction_from_file.py:65)
    NR_Reduction(config)
    assert config.ThetaShift == [0, 0]
    assert config.method_per_run == ["meantheta", "meantheta"]


def _all_unset_document():
    return SettingsDocument.from_dict({
        **_THREE_ANGLES,
        "DBname": [None, None, None],
        "ThetaShift": [None, None, None],
        "method_per_run": [None, None, None],
        "useBS": [None, None, None],
    })


@pytest.mark.parametrize("output", ["save", "normalize"])
def test_an_all_unset_default_or_broadcast_list_is_written_empty(tmp_path, output):
    doc = _all_unset_document()
    if output == "save":
        written = json.loads(doc.save(tmp_path / "out.json").read_text())
    else:
        written = doc.normalize()
    for name in ("ThetaShift", "method_per_run", "useBS"):
        assert written[name] == [], name
    # an angle-defining list is never collapsed: unset there is not a default
    assert written["DBname"] == [None, None, None]


#: One set value per list the reducer fills or broadcasts — the dimension G7's rule must cover
#: (v2, review 1568397 B-3: a ThetaShift-only test let the broadcast leg go unguarded).
_FILLS_ITSELF = [
    pytest.param(f.name, value, id=f.name)
    for f, value in (
        (fs.get("ThetaShift"), 0.1), (fs.get("useBS"), True), (fs.get("ScaleFactor"), 1.05),
        (fs.get("tof_min"), 10.0), (fs.get("tof_max"), 50000.0), (fs.get("method_per_run"), "constantQ"),
    )
]


def test_the_fills_itself_matrix_covers_every_default_and_broadcast_list():
    """A pin: when a field gains default_if_empty or broadcast_ok, add it to _FILLS_ITSELF."""
    covered = {p.values[0] for p in _FILLS_ITSELF}
    assert covered == {f.name for f in fs.FIELD_SPEC if f.default_if_empty or f.broadcast_ok}


@pytest.mark.parametrize("name, value", _FILLS_ITSELF)
def test_a_partly_set_default_list_is_a_problem_naming_the_unset_angles_and_is_saved_as_held(tmp_path, name, value):
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, name: [value, None, None]})
    lines = [m for m in doc.validate() if f"({name})" in m]
    assert len(lines) == 1
    assert "[1, 2]" in lines[0]
    saved = json.loads(doc.save(tmp_path / "out.json").read_text())[name]
    assert saved[1:] == [None, None]


def test_an_all_null_default_list_loads_as_unset_and_saves_empty(tmp_path):
    source = tmp_path / "older_editor_save.json"
    source.write_text(json.dumps({**_THREE_ANGLES, "ThetaShift": [None, None, None]}))
    doc = SettingsDocument.from_file(source)
    assert doc.get("ThetaShift") == [None, None, None]
    assert doc.validate() == []
    assert '"ThetaShift": []' in doc.save(tmp_path / "out.json").read_text()


def test_removing_the_surplus_angle_clears_the_note_and_trims_only_useBS():
    doc = _surplus_document()
    doc.remove_angle(3)
    assert doc.notes() == []
    assert doc.get("useBS") == [True, True, True]
    assert doc.get("DBname") == _THREE_ANGLES["DBname"]


def test_an_unset_background_switch_is_noted_as_the_reductions_default():
    """Plan A2: an all-unset useBS is written [], which the reducer fills with 1 (on) for
    every angle (nr_reduction_calc.py:102-103), so the panel says so."""
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, "useBS": []})
    notes = [line for line in doc.notes() if "(useBS)" in line]
    assert len(notes) == 1
    assert re.search(r"\bon\b", notes[0]), notes[0]


def test_the_boolean_default_lists_are_exactly_useBS():
    """A pin on the derivation the A2 note uses: its "on" is the reducer's fill for useBS
    (nr_reduction_calc.py:103). Before adding another, read that field's default fill."""
    assert {f.name for f in fs.FIELD_SPEC if f.default_if_empty and f.element_type == "bool"} == {"useBS"}


def test_a_per_angle_value_that_is_not_a_list_is_kept_out_of_the_counts():
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, "tof_min": 5})
    assert doc.reduction_angles == 3
    assert any("(tof_min)" in m for m in doc.validate())
    assert not any("(tof_min)" in line for line in doc.notes())


def test_an_all_unset_default_list_shorter_than_the_angles_is_not_reported_short(tmp_path):
    """Added when mutation F7 survived the battery. An all-unset default list is written []
    (G7), so the reducer fills it whatever its length, and reporting it short would cry
    wolf on a file that reduces. No test held one shorter than the angle count."""
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, "ThetaShift": [None]})
    assert not any("(ThetaShift)" in m for m in doc.validate())
    assert json.loads(doc.save(tmp_path / "out.json").read_text())["ThetaShift"] == []


# --------------------------------------------------------------------------
# editor-angle-count v2 — G8: an edit changes exactly the entry edited (review 1568397, B-1:
# one cell edit on Aug2026/REFL_231105 padded DBname to the table's 7 rows, the count went
# 3 -> 7, the notes vanished and four false problems appeared)
# --------------------------------------------------------------------------


def _aug2026_document():
    """Aug2026/REFL_231105_settings.json's shape: three angles, useBS 6 long, method_per_run 7."""
    return SettingsDocument.from_dict({**_THREE_ANGLES, "useBS": [1] * 6, "method_per_run": ["meanTheta"] * 7})


def _filled_document():
    """Every list the reducer reads holds an entry per angle, useBS one longer: the reducer-written
    reduce_settings.json shape (ThetaShift, ScaleFactor, method_per_run for the three angles), λ set too."""
    return SettingsDocument.from_dict({
        **_THREE_ANGLES, "useBS": [1] * 4, "ThetaShift": [0, 0, 0], "ScaleFactor": [1, 1, 1],
        "method_per_run": ["meanTheta"] * 3, "LambdaMin": [2.5, 2.6, 2.7],
    })


_SHAPES = {"surplus": _surplus_document, "aug2026": _aug2026_document, "filled": _filled_document}
_EDITS = {  # one list of each kind, with a value of its type
    "DBname": "x.dat", "ThetaShift": 0.01, "useBS": False, "method_per_run": "constantQ", "LambdaMin": 2.5,
}
#: The (shape, list) cells where the list holds entries the reducer reads as they are. v3 moved the
#: cells of a compact list (empty, unset, one broadcast entry, None) to G9's matrix below, which
#: asserts what the reducer reads after the edit; the "filled" shape gives every kind a cell here.
_NON_COMPACT = [
    pytest.param(shape, name, id=f"{shape}-{name}")
    for shape, names in (
        ("surplus", ("DBname", "useBS")),
        ("aug2026", ("DBname", "useBS", "method_per_run")),
        ("filled", ("DBname", "ThetaShift", "useBS", "method_per_run", "LambdaMin")),
    )
    for name in names
]


def _angles_named(lines):
    """Every angle index a problem line names."""
    return [int(i) for line in lines for group in re.findall(r"angles \[([0-9, ]*)\]", line)
            for i in group.split(",") if i.strip()]


@pytest.mark.parametrize("where", ["first", "last-angle", "surplus-row"])
@pytest.mark.parametrize("shape, name", _NON_COMPACT)
def test_an_edit_changes_exactly_the_entry_edited(tmp_path, shape, name, where):
    doc = _SHAPES[shape]()
    value = _EDITS[name]
    m, notes_before, problems_before = doc.reduction_angles, doc.notes(), doc.validate()
    index = {"first": 0, "last-angle": m - 1, "surplus-row": m}[where]
    before = doc.get(name)
    # a cell G8 decides: the list holds entries at the angles and reaches the edited row, or one short of it
    assert len(before) >= m and any(entry is not None for entry in before[:m])
    expected = list(before) + [None] * (index + 1 - len(before))
    expected[index] = value

    doc.set_angle_field(index, name, value)

    assert repr(doc.get(name)) == repr(expected)
    defines = name in fs.ANGLE_DEFINING_NAMES and index >= m
    assert doc.reduction_angles == (index + 1 if defines else m)
    problems = doc.validate()
    assert all(i < doc.reduction_angles for i in _angles_named(problems))
    if not defines:
        # nothing else changes: an edit to one list never alters what is said about another
        def others(lines):
            return [line for line in lines if f"({name})" not in line]

        assert others(problems) == others(problems_before)
    field = fs.get(name)
    if index >= m and (field.default_if_empty or field.broadcast_ok):
        # a value in a surplus row of a list the reducer fills is never read: no problem about it
        assert not any(f"({name})" in line for line in problems)
    if index < m:
        assert doc.notes() == notes_before
    first = doc.save(tmp_path / "first.json")
    # the written list, not only a round trip (which agrees with any consistent wrong encoding)
    written = json.loads(first.read_text())[name]
    assert repr(written) == repr([int(entry) if isinstance(entry, bool) else entry for entry in expected])
    second = SettingsDocument.from_file(first).save(tmp_path / "second.json")
    assert first.read_bytes() == second.read_bytes()


def test_an_edit_below_the_count_reports_nothing_about_a_surplus_row():
    """B-1's second reproduction: the reducer-written reduce_settings.json shape (ThetaShift,
    ScaleFactor and method_per_run filled for the three angles, useBS one longer)."""
    doc = SettingsDocument.from_dict({
        **_THREE_ANGLES, "useBS": [1] * 4, "ThetaShift": [0, 0, 0], "ScaleFactor": [1, 1, 1],
        "method_per_run": ["meanTheta"] * 3,
    })
    doc.set_angle_field(0, "ThetaShift", 0.01)
    assert doc.validate() == []
    assert doc.get("ThetaShift") == [0.01, 0, 0]


def test_a_broadcast_entry_expands_to_the_reductions_angles_not_the_tables_rows():
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, "useBS": [1] * 4, "method_per_run": ["meanTheta"]})
    doc.set_angle_field(0, "method_per_run", "constantQ")
    assert doc.get("method_per_run") == ["constantQ", "meanTheta", "meanTheta"]


def test_an_angle_defining_edit_in_a_surplus_row_makes_it_an_angle_and_nothing_else_grows():
    doc = _surplus_document()
    doc.set_angle_field(3, "DBname", "d.dat")
    assert doc.reduction_angles == 4
    assert doc.get("DBname") == ["db_a.dat", "db_b.dat", "db_c.dat", "d.dat"]
    assert doc.get("RB_Ymin") == [140, 141, 142]
    assert len(doc.get("useBS")) == 4


def test_a_value_for_an_empty_default_list_lands_on_the_new_angle():
    """B-3: G6's "compact, given a value -> expanded" for an EMPTY list when n > 0. Asserts the
    position: a bare [0.1] would put the new angle's value on angle 0. v3 (G9): the existing angles
    hold the reducer's own value (0, nr_reduction_calc.py:101), not unset entries."""
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, "ThetaShift": []})
    doc.add_angle(DBname="d.dat", ThetaShift=0.1)
    assert repr(doc.get("ThetaShift")) == "[0, 0, 0, 0.1]"


def test_an_unset_default_list_with_a_surplus_value_is_written_empty(tmp_path):
    """Added when frame row V2 (save() deciding "all unset" over the table's rows) and plan row N3
    (the some-unset rule over the whole list) survived the first battery. v3: an edit of an empty
    list fills the angles (G9), so the state is reached by clearing them again. The angles are all
    unset, so the reducer's default applies, and the surplus value is one it never reads."""
    doc = _surplus_document()
    doc.set_angle_field(3, "ThetaShift", 0.01)
    for angle in range(3):
        doc.set_angle_field(angle, "ThetaShift", None)
    assert repr(doc.get("ThetaShift")) == "[None, None, None, 0.01]"
    assert doc.validate() == []
    assert any("(ThetaShift)" in line and "1 extra" in line for line in doc.notes())
    assert json.loads(doc.save(tmp_path / "out.json").read_text())["ThetaShift"] == []
    assert doc.normalize()["ThetaShift"] == []


def test_an_unset_background_switch_with_a_surplus_value_is_still_noted_as_the_default(tmp_path):
    """Added when frame row V3 (the A2 note deciding "unset" over the whole list) survived. v3:
    reached by clearing the angles an edit filled (G9). Unset at every angle is what the reducer
    sees: the file is written [] and it fills 1 (on)."""
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, "method_per_run": ["meanTheta"] * 4})
    doc.set_angle_field(3, "useBS", False)
    for angle in range(3):
        doc.set_angle_field(angle, "useBS", None)
    assert repr(doc.get("useBS")) == "[None, None, None, False]"
    assert any("(useBS)" in line and re.search(r"\bon\b", line) for line in doc.notes())
    assert json.loads(doc.save(tmp_path / "out.json").read_text())["useBS"] == []


# --------------------------------------------------------------------------
# editor-angle-count v3 — G9: an edit of one angle of a compact list (one the reducer fills,
# broadcasts or derives itself) never changes what the reducer reads at another angle, and never
# leaves a file it refuses or crashes on without a line saying so. Review c286e9b: C-1
# method_per_run [] -> ['constantQ'], broadcast to every angle with no line; C-2 ThetaShift [] ->
# [0.01], IndexError at nr_reduction_calc.py:413; C-3 a λ in a surplus row turned derivation off.
# --------------------------------------------------------------------------


def _reducer_reading(path):
    """What the reducer reads at each angle of a saved file, by the reducer's own code.

    NR_Reduction runs __init__ (the method default, nr_reduction_calc.py:42-43) and _validate_config
    (the broadcast and lower-casing, :77-82; the default fills, :99-110), and raises on a file it
    refuses. Per per-angle list, one entry per angle (len(RBnum)): the repr of what it reads, so 1
    and True differ; "derived" for an optional list left None (:381-383); "IndexError" where the
    list stops short (:413 indexes ThetaShift unchecked).
    """
    from lr_reduction.nr_reduction_calc import NR_Reduction

    config = json_to_config(json.loads(pathlib.Path(path).read_text()))
    NR_Reduction(config)
    angles = range(len(config.RBnum))
    reading = {}
    for name in fs.PER_ANGLE_NAMES:
        value = getattr(config, name)
        if value is None:
            reading[name] = ["derived" for _ in angles]
        else:
            reading[name] = [repr(value[j]) if j < len(value) else "IndexError" for j in angles]
    return reading


def _m1_document():
    """One angle; useBS and method_per_run three long, so it has two surplus rows."""
    return SettingsDocument.from_dict({
        **{name: values[:1] for name, values in _THREE_ANGLES.items()},
        "useBS": [1] * 3, "method_per_run": ["meanTheta"] * 3,
    })


_G9_DOCUMENTS = {"m1": _m1_document, "surplus": _surplus_document, "aug2026": _aug2026_document}

#: Each list the reducer fills or broadcasts: the value an edit types, as the file holds it; and
#: the reducer's own value at an angle the list does not give (nr_reduction_calc.py:42-43, :99-110),
#: as the document holds it and as the file does. Stated here from the reducer, not read from
#: Field.reducer_default: the S2 tests pin that declaration to the reducer.
_G9_FILLS = {
    "ThetaShift": (0.01, 0.01, 0, 0),
    "useBS": (False, 0, True, 1),
    "ScaleFactor": (1.05, 1.05, 1, 1),
    "tof_min": (10.0, 10.0, 0, 0),
    "tof_max": (50000.0, 50000.0, 100000, 100000),
    "method_per_run": ("constantQ", "constantQ", "meanTheta", "meanTheta"),
}
_G9_LAMBDAS = {"LambdaMin": 2.5, "LambdaMax": 9.5}


def _g9_states(name, m):
    """{compact state: (the list as held, what the reducer reads at every angle while it is so —
    as held, as written)}. Unset at every angle and None are written as the reducer's own value
    too ([] and null: _encode_for_file; `if not ...` at :99-110), so they are compact as well."""
    if name in _G9_LAMBDAS:
        return {"derived": (None, None, None)}
    fill_held, fill_written = _G9_FILLS[name][2:]
    states = {
        "empty": ([], fill_held, fill_written),
        "null": (None, fill_held, fill_written),
        "unset": ([None] * m, fill_held, fill_written),
    }
    if fs.get(name).broadcast_ok:
        states["one-entry"] = (["constantTOF"], "constantTOF", "constantTOF")
    return states


def _g9_cells():
    """Every document x list x compact state x edited row: the first angle, the last, the first
    surplus row and the last one (with a gap between them where there is one)."""
    cells = []
    for document, build in _G9_DOCUMENTS.items():
        for name in [*_G9_FILLS, *_G9_LAMBDAS]:
            m = build().reduction_angles
            for state, (held, _, _) in _g9_states(name, m).items():
                doc = build()
                doc.set(name, held)
                rows = sorted({0, m - 1, m, doc.n_angles - 1} & set(range(doc.n_angles)))
                cells += [pytest.param(document, name, state, row, id=f"{document}-{name}-{state}-row{row}")
                          for row in rows]
    return cells


@pytest.mark.parametrize("document, name, state, index", _g9_cells())
def test_an_edit_of_a_compact_list_keeps_what_the_reducer_reads_at_every_other_angle(
        tmp_path, document, name, state, index):
    """S1 — plan v3's per-cell table (G9). Each cell asserts the list held and the list written,
    the lines, and the invariant itself, read by the reducer from the saved file: every angle but
    the edited one reads as it did, or a problem line names it."""
    doc = _G9_DOCUMENTS[document]()
    m = doc.reduction_angles
    held, fill_held, fill_written = _g9_states(name, m)[state]
    doc.set(name, held)
    problems_before, notes_before = doc.validate(), doc.notes()
    reading_before = _reducer_reading(doc.save(tmp_path / "before.json"))
    if name in _G9_LAMBDAS:
        value = value_written = _G9_LAMBDAS[name]
    else:
        value, value_written = _G9_FILLS[name][:2]

    doc.set_angle_field(index, name, value)

    problems, notes = doc.validate(), doc.notes()
    saved = doc.save(tmp_path / "after.json")
    written = json.loads(saved.read_text())[name]
    reading = _reducer_reading(saved)
    lines = [line for line in problems if f"({name})" in line]
    others = [j for j in range(m) if j != index]

    # the invariant: another angle reads as it did, or a problem line names it
    changed = [j for j in others if reading[name][j] != reading_before[name][j]]
    if changed:
        assert any(f"angles {changed}" in line for line in lines), (changed, problems)
    # an edit of one list says nothing new about another, and changes nothing another list reads
    assert [line for line in problems if f"({name})" not in line] == [
        line for line in problems_before if f"({name})" not in line]
    assert {k: v for k, v in reading.items() if k != name} == {
        k: v for k, v in reading_before.items() if k != name}

    if name in _G9_LAMBDAS and index >= m:
        # (b) refused: a λ only in a row the reducer never reads would end derivation at every angle
        assert doc.get(name) is None and written is None
        assert len(lines) == 1, lines
        assert f"surplus row {index + 1}" in lines[0] and "leave it derived" in lines[0]
        assert reading == reading_before and notes == notes_before
    elif name in _G9_LAMBDAS:
        # (b): written out to the reduction's angles with unset entries, which one line names
        expected = [None] * m
        expected[index] = value
        assert repr(doc.get(name)) == repr(expected) == repr(written)
        assert changed == others and len(lines) == (1 if others else 0)
        assert notes == notes_before
    else:
        # (a): every other angle holds what the reducer would have used. A broadcast list holds no
        # unset entry anywhere, surplus rows included: the reducer lower-cases all of it (:82)
        gap_held, gap_written = (fill_held, fill_written) if fs.get(name).broadcast_ok else (None, None)
        expected = [fill_held] * m + [gap_held] * (index + 1 - m)
        expected_written = [fill_written] * m + [gap_written] * (index + 1 - m)
        expected[index], expected_written[index] = value, value_written
        assert repr(doc.get(name)) == repr(expected)
        assert repr(written) == repr(expected_written)
        assert lines == [] and changed == []
        new, gone = set(notes) - set(notes_before), set(notes_before) - set(notes)
        if index >= m:
            assert len(new) == 1, new
            note = new.pop()
            assert f"({name})" in note and f"{index + 1 - m} extra" in note
        else:
            assert new == set()
        assert gone == ({line for line in notes_before if "(useBS)" in line} if name == "useBS" else set())
    if index < m:
        read_as = value_written.lower() if isinstance(value_written, str) else value_written
        assert reading[name][index] == repr(read_as)
    again = SettingsDocument.from_file(saved).save(tmp_path / "again.json")
    assert again.read_bytes() == saved.read_bytes()


@pytest.mark.parametrize("name", fs.DEFAULT_IF_EMPTY_NAMES)
def test_the_declared_reducer_default_is_what_the_reducer_fills_and_what_an_edit_writes(name):
    """S2: G9 writes Field.reducer_default at the angles an edit does not touch, so it must be the
    reducer's own fill exactly (nr_reduction_calc.py:99-110) — by repr, so 1 and True, 0 and 0.0
    differ — and it must be what the edit writes."""
    from lr_reduction.nr_reduction_calc import NR_Reduction

    config = json_to_config({**_THREE_ANGLES, name: []})
    NR_Reduction(config)
    declared = fs.get(name).reducer_default
    assert repr(getattr(config, name)) == repr([declared] * 3)
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, name: []})
    doc.set_angle_field(0, name, _G9_FILLS[name][0])
    assert repr(doc.normalize()[name][1:]) == repr([declared] * 2)


def test_the_declared_method_default_is_the_reducers():
    """S2: the reducer defaults an empty method_per_run to 'meantheta' (nr_reduction_calc.py:42-43);
    the declared default is the editor's spelling of that choice, and what an edit writes."""
    from lr_reduction.nr_reduction_calc import NR_Reduction

    config = json_to_config({**_THREE_ANGLES, "method_per_run": []})
    NR_Reduction(config)
    declared = fs.get("method_per_run").reducer_default
    assert declared in fs.get("method_per_run").allowed
    assert config.method_per_run == [declared.lower()] * 3
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, "method_per_run": []})
    doc.set_angle_field(0, "method_per_run", "constantQ")
    assert doc.normalize()["method_per_run"][1:] == [declared] * 2


def test_a_reducer_default_is_declared_for_exactly_the_lists_the_reducer_fills():
    assert {f.name for f in fs.FIELD_SPEC if f.reducer_default is not None} == {
        f.name for f in fs.FIELD_SPEC if f.default_if_empty or f.broadcast_ok}


def test_a_lambda_typed_into_a_surplus_row_of_a_derived_list_is_refused_with_a_line():
    """S3 (C-3): the list stays None — derived at every angle — and one line names the field, the
    row as its header shows it, and the remedy. Writing it out ([None] * 4 + [3.0]) ended derivation
    for the three angles, and the reducer raised at nr_reduction_calc.py:452."""
    doc = _surplus_document()
    doc.set_angle_field(3, "LambdaMin", 3.0)
    assert doc.get("LambdaMin") is None
    lines = [line for line in doc.validate() if "(LambdaMin)" in line]
    assert len(lines) == 1
    assert "surplus row 4" in lines[0] and "leave it derived" in lines[0]


@pytest.mark.parametrize(
    "gesture", ["set-an-angle", "clear-the-cell", "remove-the-row", "add-an-angle", "make-it-an-angle"])
def test_a_refused_lambda_is_reported_only_while_it_still_applies(gesture):
    """The refusal is about one gesture on one row. It goes when the field is given a value, when
    the cell is cleared, when rows move, and when the row becomes an angle: a line naming a row that
    no longer holds what the user typed into, or is no longer surplus, would be false."""
    doc = _aug2026_document()
    doc.set_angle_field(5, "LambdaMin", 3.0)
    assert any("surplus row 6" in line for line in doc.validate())
    if gesture == "set-an-angle":
        doc.set_angle_field(0, "LambdaMin", 2.5)
    elif gesture == "clear-the-cell":
        doc.set_angle_field(5, "LambdaMin", None)
    elif gesture == "remove-the-row":
        doc.remove_angle(5)
    elif gesture == "add-an-angle":
        doc.add_angle(DBname="d.dat")
    else:
        doc.set_angle_field(5, "DBname", "d.dat")
    assert not any("surplus row" in line for line in doc.validate())


def test_values_given_to_add_land_on_the_new_angle_and_the_others_read_as_before():
    """S4 (B-T1): on the surplus document (m = 3, n = 4) every value lands at index m, the
    reduction's next angle — not n, the table's next row — and each compact list's angles hold
    what G9 writes: the reducer's own value, or unset and named by a line (λ). Every earlier
    with-value test ran where m == n, so reverting m to n left them green."""
    doc = _surplus_document()
    doc.add_angle(DBname="d.dat", method_per_run="constantQ", ThetaShift=0.1, LambdaMin=2.5)
    assert doc.get("DBname")[3] == "d.dat"
    assert repr(doc.get("method_per_run")) == repr(["meanTheta"] * 3 + ["constantQ"])
    assert repr(doc.get("ThetaShift")) == "[0, 0, 0, 0.1]"
    assert repr(doc.get("LambdaMin")) == "[None, None, None, 2.5]"
    assert repr(doc.get("useBS")) == "[True, True, True, None, True]"
    problems = doc.validate()
    assert any("(LambdaMin)" in line and "angles [0, 1, 2]" in line for line in problems)
    assert not any("(method_per_run)" in line or "(ThetaShift)" in line for line in problems)


def test_an_unset_list_keeps_its_surplus_entry_behind_an_added_angle():
    """Add with a value on a list unset at every angle that holds a surplus entry (a loaded file):
    the angles get the reducer's value, the new angle the value, and the entry stays surplus."""
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, "ThetaShift": [None, None, None, 0.02]})
    doc.add_angle(DBname="d.dat", ThetaShift=0.1)
    assert repr(doc.get("ThetaShift")) == "[0, 0, 0, 0.1, 0.02]"


def test_the_unset_rule_for_an_optional_list_looks_only_at_the_reductions_angles():
    """S5 (B-T2): clearing a surplus entry of a λ list leaves every angle set, so no line; clearing
    an angle's entry is reported, by that angle."""
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, "useBS": [1] * 4, "LambdaMin": [2.5] * 4})
    doc.set_angle_field(3, "LambdaMin", None)
    assert not any("(LambdaMin)" in line for line in doc.validate())
    doc.set_angle_field(1, "LambdaMin", None)
    lines = [line for line in doc.validate() if "(LambdaMin)" in line]
    assert len(lines) == 1 and "angles [1]" in lines[0]


def test_a_broadcast_list_edited_past_its_end_holds_no_unset_entry(tmp_path):
    """The reducer lower-cases every entry of method_per_run, surplus rows included
    (nr_reduction_calc.py:82): a gap left unset there is a file it refuses. Surplus rows are no
    angle's, so the gap holds the reducer's own default, and the file reduces."""
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, "useBS": [1] * 6, "method_per_run": ["constantQ"] * 3})
    doc.set_angle_field(5, "method_per_run", "constantTOF")
    assert doc.get("method_per_run") == ["constantQ"] * 3 + ["meanTheta"] * 2 + ["constantTOF"]
    assert doc.validate() == []
    assert _reducer_reading(doc.save(tmp_path / "out.json"))["method_per_run"] == ["'constantq'"] * 3


def test_an_unset_surplus_entry_of_a_broadcast_list_is_a_problem(tmp_path):
    """Clearing a surplus method cell of an Aug2026-shaped file: the reducer never uses the row, but
    it lower-cases every entry (nr_reduction_calc.py:82) and fails on the unset one. Nothing said so
    before v3; the reducer's refusal, run below, is what makes the line true."""
    from lr_reduction.nr_reduction_calc import NR_Reduction

    doc = _aug2026_document()
    doc.set_angle_field(5, "method_per_run", None)
    lines = [line for line in doc.validate() if "(method_per_run)" in line]
    assert len(lines) == 1 and "surplus row 6" in lines[0]
    with pytest.raises(AttributeError):
        NR_Reduction(json_to_config(json.loads(doc.save(tmp_path / "out.json").read_text())))


def test_clearing_a_cell_that_holds_no_entry_changes_nothing():
    """A single broadcast entry shows on row 0 only; the rows below hold no entry. Clearing one of
    them is not an edit: expanding the list around an unset entry made a file the reducer refuses
    (nr_reduction_calc.py:82) out of one it reduced."""
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, "method_per_run": ["constantQ"]})
    doc.set_angle_field(1, "method_per_run", None)
    assert doc.get("method_per_run") == ["constantQ"]
    assert doc.validate() == []


# Added when rows of the v3 battery survived (ledger scripts/mutations-editor-angle-count.py).


def test_clearing_the_one_angle_a_broadcast_entry_shows_on_keeps_it_for_the_others():
    """F6, N2: the entry shows on row 0 only, but the reducer repeats it at every angle (:77-79).
    Clearing row 0 unsets angle 0 and no other: the list is expanded to the reduction's count (not
    the table's rows) first, and the line names angle 0 — the reducer refuses an unset entry (:82)."""
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, "useBS": [1] * 4, "method_per_run": ["constantQ"]})
    doc.set_angle_field(0, "method_per_run", None)
    assert doc.get("method_per_run") == [None, "constantQ", "constantQ"]
    lines = [line for line in doc.validate() if "(method_per_run)" in line]
    assert len(lines) == 1 and "angles [0]" in lines[0]


def test_a_refusal_does_not_come_back_after_the_field_is_set_and_cleared():
    """W11: an accepted edit forgets the field's refusal. Otherwise, once the field is derived again,
    the old line about a row the user typed into long ago would return."""
    doc = _aug2026_document()
    doc.set_angle_field(5, "LambdaMin", 3.0)
    doc.set_angle_field(0, "LambdaMin", 2.5)
    doc.set_angle_field(0, "LambdaMin", None)
    assert doc.get("LambdaMin") is None
    assert not any("surplus row" in line for line in doc.validate())


def test_a_broadcast_list_unset_at_every_angle_holds_no_unset_entry_once_written_out(tmp_path):
    """W13, W14: a method list unset at every angle (written [], so meanTheta everywhere) that holds an
    unset and a set surplus entry. An edit or an Add writes it out, and the unset surplus entry gets
    the reducer's default too: the reducer lower-cases every entry (nr_reduction_calc.py:82)."""
    shape = {**_THREE_ANGLES, "useBS": [1] * 6, "method_per_run": [None, None, None, None, "constantTOF"]}
    edited = SettingsDocument.from_dict(shape)
    edited.set_angle_field(0, "method_per_run", "constantQ")
    assert edited.get("method_per_run") == ["constantQ"] + ["meanTheta"] * 3 + ["constantTOF"]
    added = SettingsDocument.from_dict(shape)
    added.add_angle(DBname="d.dat", method_per_run="constantQ")
    assert added.get("method_per_run") == ["meanTheta"] * 3 + ["constantQ", "meanTheta", "constantTOF"]
    for doc in (edited, added):
        assert not any("(method_per_run)" in line for line in doc.validate())
    assert _reducer_reading(edited.save(tmp_path / "out.json"))["method_per_run"] == [
        "'constantq'", "'meantheta'", "'meantheta'"]


def test_a_ragged_broadcast_list_padded_past_its_end_leaves_the_missing_angle_unset_and_named():
    """W16: padding fills a broadcast list's SURPLUS rows with the reducer's default, never an angle.
    Angle 2 had no method (the list was short, which the panel already reported); filling it would
    choose one for the scientist, silently."""
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, "useBS": [1] * 6, "method_per_run": ["constantQ"] * 2})
    doc.set_angle_field(5, "method_per_run", "constantTOF")
    assert doc.get("method_per_run") == ["constantQ", "constantQ", None, "meanTheta", "meanTheta", "constantTOF"]
    assert any("(method_per_run)" in line and "angles [2]" in line for line in doc.validate())


def test_a_broadcast_list_unset_at_every_angle_is_no_problem_for_its_unset_surplus_entry(tmp_path):
    """W18: unset at every angle, it is written [] and the reducer uses meanTheta; no entry of it
    reaches the reducer, so its unset surplus entry is no problem."""
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, "useBS": [1] * 4, "method_per_run": [None] * 4})
    assert not any("(method_per_run)" in line for line in doc.validate())
    assert json.loads(doc.save(tmp_path / "out.json").read_text())["method_per_run"] == []


def test_an_explicit_none_given_to_add_leaves_a_compact_list_compact():
    """W23: None is "no value". Writing a derived λ out as [None, None, None, None] gave a list that
    reports a length and carries no values (web_report.py:547 indexes it), and a line."""
    doc = _surplus_document()
    doc.add_angle(DBname="d.dat", LambdaMin=None, ThetaShift=None)
    assert doc.get("LambdaMin") is None
    assert doc.get("ThetaShift") == []
    assert not any("(LambdaMin)" in line or "(ThetaShift)" in line for line in doc.validate())


# --------------------------------------------------------------------------
# editor-combos — what a direct-beam cell offers (C4, C5: the *.txt/*.dat names in the folder the
# reducer joins DBname to), and what a compact list implies at an angle it does not hold (C7)
# --------------------------------------------------------------------------


def test_only_the_direct_beam_column_offers_candidates_and_from_the_folder_the_reducer_reads():
    """nr_reduction_calc.py:402 calls tools.load_db_file(config.DBpath, config.DBname[i]): a DBname entry
    is a file name in DBpath, so that is the folder its candidates come from."""
    assert {f.name: f.candidates_folder for f in fs.FIELD_SPEC if f.candidates_folder} == {"DBname": "DBpath"}
    assert SettingsDocument().candidates("ThetaShift") == ([], 0)  # no folder declared, nothing offered


def _direct_beam_folder(tmp_path, names):
    folder = tmp_path / "transmission"
    folder.mkdir()
    for name in names:
        (folder / name).write_text("")
    return folder


def test_the_direct_beam_candidates_are_the_folders_txt_and_dat_files_sorted(tmp_path):
    """Names are stored verbatim (spaces, non-ASCII); another suffix and a sub-folder are not offered."""
    folder = _direct_beam_folder(tmp_path, ["db_b.dat", "db_a.txt", "notes.md", "db 1 é.dat", "DB_Z.DAT"])
    (folder / "sub.dat").mkdir()
    doc = SettingsDocument.from_dict({"_DBpath_override": str(folder)})
    assert doc.candidates("DBname") == (["DB_Z.DAT", "db 1 é.dat", "db_a.txt", "db_b.dat"], 4)


@pytest.mark.parametrize(
    "state", ["missing", "a-file", "unreadable", "scandir-raises", "unresolvable-path", "null-byte"])
def test_the_direct_beam_candidates_are_empty_whenever_the_folder_cannot_be_listed(tmp_path, monkeypatch, state):
    """Every failure is an empty list, so a cell never raises into a Qt slot. The folder is on a facility
    mount (F6), where any of these is ordinary; "unresolvable-path" is experiment_id None, for which the
    config's DBpath property itself raises TypeError (Path / None)."""
    folder = _direct_beam_folder(tmp_path, ["db_a.dat"])
    values = {"_DBpath_override": str(folder)}
    if state == "missing":
        values["_DBpath_override"] = str(tmp_path / "absent")
    elif state == "a-file":
        values["_DBpath_override"] = str(folder / "db_a.dat")
    elif state == "unreadable":
        if os.geteuid() == 0:
            pytest.skip("root lists a mode-000 folder")
        folder.chmod(0)
    elif state == "scandir-raises":
        def scandir(_path):
            raise OSError("stale file handle")

        monkeypatch.setattr(os, "scandir", scandir)
    elif state == "null-byte":
        values["_DBpath_override"] = str(folder) + "\x00"  # os.scandir raises ValueError
    else:
        values = {"experiment_id": None}
    try:
        assert SettingsDocument.from_dict(values).candidates("DBname") == ([], 0)
    finally:
        folder.chmod(0o755)


def test_the_direct_beam_candidates_are_capped_and_say_how_many_there_were(tmp_path):
    folder = _direct_beam_folder(tmp_path, [f"db_{k:02d}.dat" for k in range(7)])
    doc = SettingsDocument.from_dict({"_DBpath_override": str(folder)})
    assert doc.candidates("DBname", limit=5) == ([f"db_{k:02d}.dat" for k in range(5)], 7)


def test_a_compact_list_implies_the_reductions_value_at_each_angle_it_does_not_hold():
    """C7's model half: where a compact list gives nothing, the reduction uses the broadcast entry
    (nr_reduction_calc.py:77-79) or its own default (:99-110), held in the document's spelling (True for
    useBS's 1)."""
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, "method_per_run": ["constantQ"], "useBS": []})
    assert [doc.implied_entry(k, "method_per_run") for k in range(3)] == [None, "constantQ", "constantQ"]
    assert repr([doc.implied_entry(k, "useBS") for k in range(3)]) == "[True, True, True]"


def test_nothing_is_implied_where_the_reduction_has_no_value_of_its_own():
    """A surplus row (never read); a held entry; an entry left unset in a list the reducer reads as held
    (it reads None there, which is a problem line, not a default); a derived λ; an angle-defining list."""
    surplus = _surplus_document()  # useBS x4, method_per_run []
    assert surplus.implied_entry(2, "method_per_run") == "meanTheta"
    assert surplus.implied_entry(3, "method_per_run") is None
    ragged = SettingsDocument.from_dict({**_THREE_ANGLES, "method_per_run": ["constantQ", None, "constantTOF"]})
    assert ragged.implied_entry(0, "method_per_run") is None
    assert ragged.implied_entry(1, "method_per_run") is None
    assert surplus.implied_entry(0, "LambdaMin") is None
    assert surplus.implied_entry(0, "DBname") is None


def test_an_entry_that_cannot_be_examined_is_left_out_and_the_rest_are_listed(tmp_path, monkeypatch):
    """Added at GREEN, for a branch the RED set did not construct: on a facility mount one entry's stat can
    fail (a stale handle) while the folder lists. That entry is not offered; the listing goes on."""
    class Entry:
        def __init__(self, name, broken=False):
            self.name, self._broken = name, broken

        def is_file(self):
            if self._broken:
                raise OSError("stale file handle")
            return True

    class Listing:
        def __enter__(self):
            return iter([Entry("db_b.dat"), Entry("db_x.dat", broken=True), Entry("db_a.txt")])

        def __exit__(self, *exc):
            return False

    monkeypatch.setattr(os, "scandir", lambda _path: Listing())
    doc = SettingsDocument.from_dict({"_DBpath_override": str(tmp_path)})
    assert doc.candidates("DBname") == (["db_a.txt", "db_b.dat"], 2)


def test_a_list_unset_at_every_angle_implies_the_reductions_default():
    """Added when frame row F34 survived: written [] (G7), so the reducer uses its default there; the
    tests above held only [] and one broadcast entry, never unset entries."""
    doc = SettingsDocument.from_dict({**_THREE_ANGLES, "method_per_run": [None, None, None]})
    assert [doc.implied_entry(k, "method_per_run") for k in range(3)] == ["meanTheta"] * 3


# --------------------------------------------------------------------------
# editor-defaults-and-theta — a new file starts at gaussian / 1.0 (the library keeps rectangular / 0.8); the
# theta enumeration is held canonical when the reducer accepts the value, and kept as loaded when it does not
# --------------------------------------------------------------------------


def test_a_document_the_editor_creates_starts_at_gaussian_and_1():
    """U1, D1: the starting values are part of the seed, so nothing shows as changed in a new file."""
    doc = SettingsDocument.for_new_file()
    assert doc.get("DetResFn") == "gaussian"
    assert doc.get("DetSigma") == 1.0 and type(doc.get("DetSigma")) is float
    assert doc.changed_vs_seed() == {}


def test_the_library_and_a_bare_document_keep_rectangular_and_0_8():
    """U2, D2: only a document the editor creates starts at the new values."""
    config = NRReductionConfig()
    assert (config.DetResFn, config.DetSigma) == ("rectangular", 0.8)
    assert (fs.get("DetResFn").default, fs.get("DetSigma").default) == ("rectangular", 0.8)
    doc = SettingsDocument()
    assert (doc.get("DetResFn"), doc.get("DetSigma")) == ("rectangular", 0.8)


def test_a_loaded_file_without_the_resolution_keys_holds_what_its_reduction_uses(tmp_path):
    """U3, D2, A2: a file that omits them is reduced with the library's rectangular / 0.8 (json_to_config starts
    from NRReductionConfig()), so that is what the editor shows, not the new file's starting values."""
    seed = tmp_path / "s.json"
    seed.write_text(json.dumps({"Sname": "old"}))
    doc = SettingsDocument.from_file(seed)
    assert (doc.get("DetResFn"), doc.get("DetSigma")) == ("rectangular", 0.8)


def test_the_fields_with_a_starting_value_are_exactly_the_resolution_pair():
    """U4: declared on the Field and derived into EDITOR_START_NAMES; every other field starts at its default."""
    assert {name: fs.get(name).starting_value() for name in fs.EDITOR_START_NAMES} == {
        "DetResFn": "gaussian", "DetSigma": 1.0,
    }
    assert all(f.starting_value() == f.default for f in fs.FIELD_SPEC if f.name not in fs.EDITOR_START_NAMES)


_ABSENT = object()

# (id, the value in the file, the value held after load, reported?) — plan §3's table, one row per spelling.
_THETA_LOADS = [
    ("absent", _ABSENT, False, False),
    ("false", False, False, False),
    ("null", None, False, False),
    ("zero", 0, False, False),
    ("empty", "", False, False),
    ("true", True, "detector_angle", False),
    ("detector_angle", "detector_angle", "detector_angle", False),
    ("Detector_Angle", "Detector_Angle", "detector_angle", False),
    ("DETECTOR_ANGLE", "DETECTOR_ANGLE", "detector_angle", False),
    ("sample_angle", "sample_angle", "sample_angle", False),
    ("Sample_Angle", "Sample_Angle", "sample_angle", False),
    ("SAMPLE_ANGLE", "SAMPLE_ANGLE", "sample_angle", False),
    ("string-true", "true", "true", True),
    ("string-TRUE", "TRUE", "TRUE", True),
    ("string-yes", "yes", "yes", True),
    ("one", 1, 1, True),
    ("string-True", "True", "True", True),
    ("detector", "detector", "detector", True),
    ("sample", "sample", "sample", True),
    ("a-list", ["detector_angle"], ["detector_angle"], True),
]


def _theta_document(tmp_path, loaded):
    seed = tmp_path / "seed.json"
    seed.write_text(json.dumps({} if loaded is _ABSENT else {"useCalcTheta": loaded}))
    return SettingsDocument.from_file(seed)


@pytest.mark.parametrize("loaded, held, reported", [pytest.param(*row[1:], id=row[0]) for row in _THETA_LOADS])
def test_loading_holds_a_theta_value_the_reducer_accepts_canonically_and_keeps_the_rest(tmp_path, loaded, held,
                                                                                        reported):
    """U5, D5, D6. What the reducer accepts is held as the value it acts on: True and any case of a name as the
    lower-case name, any falsy value as False. What it rejects (it raises on "true") is kept exactly as loaded
    and reported with the three forms it accepts. Either way nothing shows as changed, and a save writes what
    is held."""
    doc = _theta_document(tmp_path, loaded)
    value = doc.get("useCalcTheta")
    assert value == held and type(value) is type(held)
    lines = [message for message in doc.validate() if "useCalcTheta" in message]
    assert bool(lines) is reported, lines
    if reported:
        assert all(form in lines[0] for form in ("False", "True", "trust sample angle")), lines[0]
    assert doc.changed_vs_seed() == {}
    saved = json.loads(doc.save(tmp_path / "out.json").read_text())["useCalcTheta"]
    assert saved == held and type(saved) is type(held)


@pytest.mark.parametrize("loaded", [pytest.param(row[1], id=row[0]) for row in _THETA_LOADS
                                    if row[1] is not _ABSENT])
def test_the_reduction_reads_the_saved_theta_value_as_it_read_the_loaded_one(tmp_path, loaded):
    """U7, D7: the reducer's own normalisation (NR_Reduction._validate_config, nr_reduction_calc.py:90-97) does
    the same with the saved value as with the loaded one, for every spelling. It reads both as the same name,
    or both as off, or rejects both: a canonicalisation it does not perform itself ("sample" -> "sample_angle")
    would turn a file it rejects into one it reduces. Off is compared by truthiness, which is all its readers
    test (:92, :543, :574; :550 compares a name), so no reader sees which falsy value it was."""

    def outcome(value):
        reduction = _bare_reduction(RBnum=[1], DBname=["db.dat"], RB_Ymin=[1], RB_Ymax=[2], useCalcTheta=value)
        try:
            reduction._validate_config()
        except (ValueError, AttributeError) as exc:  # AttributeError: .lower() on a non-string
            return ("rejected", type(exc).__name__)
        return ("reads", reduction.config.useCalcTheta or False)

    doc = _theta_document(tmp_path, loaded)
    saved = json.loads(doc.save(tmp_path / "out.json").read_text())["useCalcTheta"]
    assert outcome(saved) == outcome(loaded)


def test_the_readers_of_use_calc_theta_are_the_ones_d7_cleared():
    """U7's pin, D7. Holding a falsy value as False and True as "detector_angle" is safe only while no reader of
    useCalcTheta tells them apart. These are its mentions in the library at dispatch (reader, writer or carrier,
    each read): nr_reduction_calc (the rule and three reads), web_report (prints it), new_reduction_from_template
    (copies it into a template), nr_reduction_config (the default) and the example scripts (writers). A new mention
    changes this count: read it, and remove the canonicalisation it can see from D5, before updating the count."""
    library = pathlib.Path(fs.__file__).parent
    declarations = {"field_spec.py", "settings_document.py", "reduction_domains.py"}
    mentions = {
        path.name: count
        for path in sorted(library.glob("*.py"))
        if path.name not in declarations and (count := path.read_text().count("useCalcTheta"))
    }
    assert mentions == {
        "EBW_nr_reduction_test.py": 2,
        "example_nr_from_template.py": 2,
        "example_nr_reduction.py": 5,
        "new_reduction_from_template.py": 3,
        "new_reduction_template_reader.py": 1,
        "nr_reduction_calc.py": 11,
        "nr_reduction_config.py": 1,
        "web_report.py": 1,
    }


def test_every_theta_choice_has_one_label_and_every_label_maps_back():
    """U6, D3, D4: one mapping, declared on the Field, read in both directions; False is an entry of its own."""
    field = fs.get("useCalcTheta")
    assert field.label == "Apply theta calculation"
    assert [label for _, label in field.choice_labels] == ["False", "True", "trust sample angle"]
    assert {stored for stored, _ in field.choice_labels} == {False, *fs.CALC_THETA_CHOICES}
    for stored, label in field.choice_labels:
        assert field.label_for(stored) == label
        back = field.value_for(label)
        assert back == stored and type(back) is type(stored)
    assert field.label_for(True) == "True" and field.label_for("Sample_Angle") == "trust sample angle"
    assert field.label_for(None) == "False" and field.label_for("true") is None


def test_a_stored_choice_without_a_label_fails_the_check_run_at_import():
    """§5, pathological: a choice added to CALC_THETA_CHOICES without a label fails at import (the check runs
    on every field there), never an unlabelled entry. Two choices with one label fail too."""
    import dataclasses

    field = fs.get("useCalcTheta")
    fs._check_choice_labels(field)
    with pytest.raises(ValueError, match="useCalcTheta"):
        fs._check_choice_labels(dataclasses.replace(field, allowed=(*field.allowed, "fitted_angle")))
    with pytest.raises(ValueError, match="share a label"):
        fs._check_choice_labels(dataclasses.replace(
            field, choice_labels=((False, "False"), ("detector_angle", "True"), ("sample_angle", "True"))))
