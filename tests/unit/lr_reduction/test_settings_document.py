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
import subprocess
import sys
import textwrap

import pytest

from lr_reduction import field_spec as fs
from lr_reduction.new_reduction_from_file import json_to_config
from lr_reduction.nr_reduction_config import NRReductionConfig
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
    """All 13 move together, or the arrays silently desynchronise.

    Growing only the obvious eight is the defect this slug exists to prevent:
    `RBnum` and `ScaleFactor` are easy to forget, and a short array shifts every
    subsequent angle's settings by one.
    """
    doc = SettingsDocument()
    doc.add_angle()
    doc.add_angle()
    lengths = {name: len(doc.get(name)) for name in fs.PER_ANGLE_NAMES if doc.get(name) is not None}
    assert set(lengths.values()) == {2}
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
