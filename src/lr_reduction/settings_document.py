"""The editable settings model behind the settings-editor tab (T2).

Wraps one :class:`~lr_reduction.nr_reduction_config.NRReductionConfig` and the
operations an editor needs: seed it from a JSON settings file, a pre-reduced
``.dat`` header, or defaults; add and remove angles; validate against
:mod:`lr_reduction.field_spec`; report what changed; and normalize for handing
back to the reduction.

**Qt-free on purpose** — see the note in :mod:`lr_reduction.field_spec`. The
view is a thin layer over this class, so nearly all of the editor's behaviour
is testable without a display.

Two design points worth stating, because both differ from the obvious reading:

*One index per added angle.* ``add_angle`` puts the new angle at the
reduction's next index (``reduction_angles``) in every list an angle is read
from, so a short array cannot shift every subsequent angle's settings by one —
which would happen silently, since nothing in the config class enforces equal
lengths. Lists the reducer fills, broadcasts or derives stay compact when no
value is given (empty, one broadcast entry, ``None``): the reducer expands them
itself. An edit (``set_angle_field``) changes what the reducer reads at the
angle edited and at no other: in a list it reads as held, exactly the entry
edited; in a compact one, the other angles are written out with what the
reducer read there before (``Field.reducer_default``, or the broadcast entry),
except a derived λ, whose other angles are left unset and reported — and a λ
typed into a row the reduction does not use is refused.

*Validation is not the same as the equal-length invariant.* The reducer
deliberately broadcasts a single ``method_per_run`` entry across all angles
(``nr_reduction_calc.py:77-79``) and defaults an empty one to ``meanTheta``
(``:42-43``). A validator that demanded strict equal lengths everywhere would
reject configurations the reducer accepts, so the broadcastable cases are
exempted here rather than "fixed".
"""

import copy
import json
import os
import tempfile
from pathlib import Path

from lr_reduction import field_spec as fs
from lr_reduction.new_reduction_from_file import json_to_config, load_from_file
from lr_reduction.nr_reduction_config import NRReductionConfig
from lr_reduction.save_reduced_data import make_json_safe

#: Cap on the problems validate() returns. A pathological file (a million
#: angles) otherwise builds a multi-megabyte list on every keystroke, on the GUI
#: thread. Truncation is announced, never silent.
MAX_REPORTED_PROBLEMS = 200

#: Cap on the file names a cell offers (``SettingsDocument.candidates``). A
#: direct-beam folder holds tens of files; a folder of tens of thousands is
#: listed once, sorted, and cut here, and the caller says so.
MAX_CANDIDATES = 2000

#: The file types a direct-beam entry names: the reducer loads ``.txt`` and
#: ``.dat`` direct-beam files. Matched case-insensitively.
CANDIDATE_SUFFIXES = (".txt", ".dat")

#: Returned by ``SettingsDocument._compact_reading`` for a list the reducer reads
#: as it is held. ``None`` cannot mark that: it is what a derived λ reads.
_NOT_COMPACT = object()


def _is_file(entry):
    """``entry.is_file()``, with a failed stat counted as "not a file" rather than ending the listing."""
    try:
        return entry.is_file()
    except OSError:
        return False


def _held_as_bool(value):
    """``value`` as a ``bool`` if it is a boolean spelling, otherwise unchanged."""
    boolean = fs.as_boolean(value)
    return value if boolean is None else boolean


class SettingsDocument:
    """One editable reduction configuration."""

    def __init__(self, config=None):
        self._config = config if config is not None else NRReductionConfig()
        # The state this document was seeded from, for changed_vs_seed(). Copied
        # so later edits cannot reach back and rewrite the baseline.
        self._seed = copy.deepcopy(self._config.__dict__)
        # {field name: row} of the λ edits set_angle_field refused, which
        # validate() reports while they still apply. Not part of the settings.
        self._refused = {}

    # -- construction ------------------------------------------------------

    @classmethod
    def for_new_file(cls):
        """A document for a settings file the editor starts from nothing.

        It holds each field's starting value (``Field.starting_value``): the
        library's default, except where the instrument's current operation
        differs (``fs.EDITOR_START_NAMES``: ``DetResFn`` ``gaussian``,
        ``DetSigma`` ``1.0``). They go in before the seed is taken, so a new
        file shows nothing as changed. Every other document keeps the library's
        values: ``SettingsDocument()``, and a loaded file that omits a field,
        which holds what reducing that file would use.
        """
        config = NRReductionConfig()
        for name in fs.EDITOR_START_NAMES:
            setattr(config, name, fs.get(name).starting_value())
        return cls(config)

    @classmethod
    def from_dict(cls, values):
        """Build from a settings mapping, reporting an unknown key by name.

        ``json_to_config`` raises a bare ``AttributeError`` naming the key; it
        is re-raised as ``ValueError`` because from the editor's point of view
        this is a bad *file*, not a programming error, and the message has to
        reach the scientist.
        """
        try:
            config = json_to_config(values)
        except AttributeError as exc:
            raise ValueError(f"Not a valid reduction setting: {exc}") from exc
        cls._migrate_legacy(config)
        # Before cls(config): the seed is taken there, and it must hold the
        # same spelling as the document, or "Changed from the seed" shows
        # [1, 1, 0] -> [True, False, False] for a one-cell edit.
        cls._canonicalize_booleans(config)
        return cls(config)

    @staticmethod
    def _canonicalize_booleans(config):
        """Hold the entries of every integer-encoded list as ``bool``, whatever spelling the file used.

        The reducer writes ``useBS`` as ``1``/``0`` and reads it by truthiness or
        ``== 1`` (see :func:`~lr_reduction.field_spec.as_boolean`), so a loaded
        integer there is a boolean in all but spelling. Converting it on load
        gives the panel one spelling to show and gives the view a real ``bool``
        to bind. Anything that is not a boolean spelling is left exactly as
        loaded, for ``validate()`` to report. A ``useBS`` that is not a list at
        all is left alone too.

        Scalar booleans are not touched. The reducer reads ``useGravity`` with
        ``is True`` (``nr_reduction_calc.py:1079``), so turning a hand-written
        ``1`` into ``True`` here switched gravity correction on (review
        8b62952). Which fields are canonicalized is ``Field.int_encoded``.
        """
        for field in fs.FIELD_SPEC:
            if not field.int_encoded:
                continue
            value = getattr(config, field.name)
            if isinstance(value, list):
                setattr(config, field.name, [_held_as_bool(entry) for entry in value])

    @staticmethod
    def _migrate_legacy(config):
        """Hold a tri-state choice as the value the reducer acts on, when the reducer accepts it.

        ``useCalcTheta = True`` is the first case: the reducer maps it to
        ``detector_angle`` and then works normally
        (``nr_reduction_calc``, ``NRReduction.__init__``). Reporting it as a
        problem would cry wolf on a file that reduces perfectly well; leaving it
        alone would keep re-saving the deprecated spelling. Migrating it on load
        does what the reducer would have done, so the panel stays quiet and the
        file the scientist saves is explicit.

        The reducer also lower-cases a name and reads any falsy value as off
        (``nr_reduction_calc.py:92-97``), so a name in any case is held as its
        declared spelling and a falsy value as ``False`` (``Field.canonical_choice``).
        No reader of the field tells those apart: its readers test truthiness
        or compare a name, pinned by
        ``test_the_readers_of_use_calc_theta_are_the_ones_d7_cleared``. A value the
        reducer rejects is left as loaded, for ``validate()`` to report.
        """
        for field in fs.FIELD_SPEC:
            if field.falsy_means_off and field.allowed:
                setattr(config, field.name, field.canonical_choice(getattr(config, field.name, None)))

    @classmethod
    def from_file(cls, path):
        """Seed from a ``.json`` settings file or a pre-reduced ``.dat`` header.

        Both are read through the existing ``load_from_file`` rather than
        reimplemented here: the ``.dat`` seed is the ``# Config:`` header line
        the reduction itself writes, and duplicating that parser is how the two
        would drift apart.
        """
        loaded = load_from_file(Path(path))
        values = loaded.get("config")
        if values is None:
            raise ValueError(f"No reduction settings found in {path}")
        return cls.from_dict(values)

    # -- scalar access -----------------------------------------------------

    @property
    def config(self):
        """The wrapped config. The reduction takes this object."""
        return self._config

    def get(self, name):
        return getattr(self._config, fs.get(name).name)

    def set(self, name, value):
        """Set a field. Unknown names raise rather than being silently stored."""
        setattr(self._config, fs.get(name).name, value)

    def to_dict(self):
        """The document's in-memory state, as a plain dict."""
        return dict(self._config.__dict__)

    # -- angles ------------------------------------------------------------

    @property
    def n_angles(self):
        """Rows the table shows: the longest per-angle field.

        This counts every index at which any list holds an entry, so nothing a
        file carries is hidden. ``reduction_angles`` is the reduction's count;
        rows beyond it are surplus.
        """
        lengths = [
            len(self.get(name))
            for name in fs.PER_ANGLE_NAMES
            if isinstance(self.get(name), list)
        ]
        return max(lengths) if lengths else 0

    def _defining_length(self):
        """The longest angle-defining list, or 0 when none holds an entry."""
        lengths = [
            len(self.get(name))
            for name in fs.ANGLE_DEFINING_NAMES
            if isinstance(self.get(name), list)
        ]
        return max(lengths, default=0)

    @property
    def reduction_angles(self):
        """The number of angles a reduction of this document will use.

        The reducer counts by ``RBnum`` and requires every list it indexes with
        no fallback not to be shorter (``nr_reduction_calc.py:61-75``), so the
        count is the longest angle-defining list
        (``field_spec.ANGLE_DEFINING_NAMES``). Entries beyond it in any other
        list are never read: they are surplus, noted rather than reported. When
        no angle-defining list holds an entry the count falls back to
        ``n_angles``, so a document holding only optional columns still has
        angles to show.
        """
        return self._defining_length() or self.n_angles

    @staticmethod
    def _compact_reading(field, current, count):
        """What the reducer reads at every angle while ``current`` is compact, or ``_NOT_COMPACT``.

        Compact is a state the file holds as the reducer's own spelling of "fill
        this in", so the reducer, not the list, decides every angle:

        * an optional list left ``None``: derived per angle
          (``nr_reduction_calc.py:381-383``). Reads ``None`` here;
        * a ``default_if_empty`` or ``broadcast_ok`` list that is ``None``,
          empty, or unset at every angle below ``count``: written ``null`` or
          ``[]`` (``_encode_for_file``), so the reducer's own value applies
          (``Field.reducer_default``; ``:42-43``, ``:99-110``), held in the
          document's spelling (``True`` for ``useBS``'s ``1``);
        * a single broadcast entry: the reducer repeats it (``:77-79``).
        """
        if field.optional_list:
            return None if current is None else _NOT_COMPACT
        if not (field.default_if_empty or field.broadcast_ok):
            return _NOT_COMPACT
        holds = isinstance(current, (list, tuple))
        if current is None or (holds and all(entry is None for entry in current[:count])):
            return _held_as_bool(field.reducer_default) if field.int_encoded else field.reducer_default
        if field.broadcast_ok and holds and len(current) == 1:
            return current[0]
        return _NOT_COMPACT

    def add_angle(self, **values):
        """Add one angle at index ``reduction_angles`` — the reduction's next one — in every list it touches.

        Appending to each list at its own end misaligned a file whose lists
        differ in length: the new values landed in different rows, and the
        author's ``useBS`` sat beside an old surplus entry. Appending after the
        surplus rows instead (v1) turned them into angles with no required
        values. The new angle goes directly after the last real angle, and any
        surplus entries move down one row and stay surplus. By the list's state
        before the gesture:

        * a compact list the reduction accepts as it is is left alone when no
          value is supplied. Compact means an empty ``default_if_empty`` list
          (the reducer fills it, ``nr_reduction_calc.py:99-110``), a single
          ``broadcast_ok`` entry (the reducer repeats it, ``:77-79``), or a
          ``None`` optional list ("derive it from the chopper ranges",
          ``:381-383``; materialising it as ``[None, None]`` would give a list
          that reports a length and carries no values, and ``web_report.py:547``
          indexes it);
        * a compact list given a value (``_compact_reading``: also one unset at
          every angle, or ``None``) is written out the way an edit writes it
          (G9): the existing angles hold what the reducer read there before —
          its own value, or the broadcast entry; a derived λ's are left unset,
          for ``validate()`` to name — and the value sits at index
          ``reduction_angles``, with any surplus entries after it;
        * any other list gets the new entry at index ``reduction_angles``: a
          shorter list is first padded with unset entries to that index, and a
          longer one's surplus entries follow the new entry;
        * a per-angle value that is not a list is left as it is, for
          ``validate()`` to report.
        """
        m = self.reduction_angles
        self._refused.clear()  # rows move; a refusal names a row
        for name in fs.PER_ANGLE_NAMES:
            field = fs.get(name)
            current = self.get(name)
            supplied = values.get(name) is not None
            reading = self._compact_reading(field, current, m)
            if supplied and reading is not _NOT_COMPACT:
                surplus = list(current[m:]) if isinstance(current, (list, tuple)) else []
                if field.broadcast_ok:  # no unset entry anywhere in it (:82)
                    surplus = [reading if entry is None else entry for entry in surplus]
                self.set(name, [reading] * m + [values[name]] + surplus)
                continue
            if not isinstance(current, (list, tuple)):
                continue
            current = list(current)
            # Nothing supplied (a compact list given a value was written out above):
            # a list the reducer expands as it is (empty, or one broadcast entry)
            # stays so. One unset at every angle still takes the new unset entry
            # at m, so that its surplus entries move down with the rest.
            compact = (field.default_if_empty and not current) or (field.broadcast_ok and len(current) <= 1)
            if compact:
                continue
            angles = current[:m] + [None] * (m - len(current[:m]))
            self.set(name, angles + [values.get(name)] + current[m:])

    def remove_angle(self, index):
        """Remove one angle from every per-angle field."""
        if not 0 <= index < self.n_angles:
            raise IndexError(f"No angle at index {index} (have {self.n_angles})")
        self._refused.clear()  # rows move; a refusal names a row
        for name in fs.PER_ANGLE_NAMES:
            current = self.get(name)
            if isinstance(current, list) and index < len(current):
                self.set(name, current[:index] + current[index + 1 :])

    def set_angle_field(self, index, name, value):
        """Set one angle's value for one field.

        The index is explicit and mandatory. The editor's table must pass the
        row it is acting on; there is deliberately no notion of a "current row"
        here to fall back on, which is the shape the active-row-as-hidden-input
        bug takes in reduction GUIs.

        What the reducer reads changes at the angle edited and at no other:

        * in a list it reads as held, exactly the entry edited changes;
        * in a compact list (``_compact_reading``) the angles the edit does not
          touch are written out with what the reducer read there before — its
          own value (``Field.reducer_default``) or the broadcast entry — so they
          reduce as before. A derived λ has no spelling for "derive this one":
          its other angles are left unset, and ``validate()`` names them. A λ
          typed into a row the reduction does not use is refused, and
          ``validate()`` says why;
        * clearing a cell that holds no entry is not an edit, and changes nothing.
        """
        field = fs.get(name)
        if not field.per_angle:
            raise KeyError(f"{name} is not a per-angle field")
        current = self.get(name)
        if not 0 <= index < self.n_angles:
            raise IndexError(f"No angle at index {index} (have {self.n_angles})")
        # isinstance, not len(): a bare string has a length, and list("abc") would
        # explode it into ['a', 'b', 'c'] instead of replacing a wrong type.
        holds = isinstance(current, (list, tuple))
        if value is None and not (holds and index < len(current) and current[index] is not None):
            # A one-entry broadcast shows on row 0 only, so the rows below it look
            # empty. "Clearing" one used to expand the list around an unset entry: a
            # file the reducer refuses, since it lower-cases every entry (:82).
            self._refused.pop(name, None)
            return
        count = self.reduction_angles
        reading = self._compact_reading(field, current, count)
        if reading is not _NOT_COMPACT and value is not None:
            if field.optional_list and index >= count:
                # None is "derive it" at every angle. A value only in a row the
                # reducer never reads would make it a list with nothing at the real
                # angles, which it then reads (:452 raised; review c286e9b, C-3).
                self._refused[name] = index
                return
            # Padded with v2's unset entries instead, one edit of an empty list
            # changed every other angle: [] -> ['constantQ'] was broadcast to all of
            # them, silently, and [] -> [0.01] raised IndexError at :413 (C-1, C-2).
            updated = [reading] * count + (list(current[count:]) if holds else [])
            gap = reading if field.broadcast_ok else None
            if field.broadcast_ok:  # no unset entry anywhere in it, surplus included (:82)
                updated = [reading if entry is None else entry for entry in updated]
            updated += [gap] * (index + 1 - len(updated))
        else:
            # A single broadcast entry stands for the reduction's angles (:77-79):
            # clearing one of them expands it to that count first, not to the table's
            # rows, which are surplus past it.
            updated = list(current) if holds else []
            if field.broadcast_ok and len(updated) == 1:
                updated *= max(count, 1)
            # Padded only far enough to reach the edited index, never to the table's
            # row count: v1 padded to n_angles, and on a file with surplus rows one
            # cell edit grew every edited list into them, which raised the reduction's
            # count (3 -> 7 on a real file) and reported problems that did not exist
            # (review 1568397). Padding a SHORT list rather than indexing past it
            # matters too: an IndexError in a Qt slot reaches qFatal() and kills the
            # launcher. A broadcast list's surplus rows are no angle's, but the reducer
            # lower-cases every entry, so padding there holds its default.
            updated += [
                field.reducer_default if field.broadcast_ok and position >= count else None
                for position in range(len(updated), index + 1)
            ]
        updated[index] = value
        # An optional list that is emptied of every value goes back to None —
        # "derive it from the chopper ranges". Without this, touching one Lambda
        # cell is a one-way door out of that state for the life of the document.
        if field.optional_list and all(entry is None for entry in updated):
            updated = None
        self._refused.pop(name, None)
        self.set(name, updated)

    def implied_entry(self, index, name):
        """The value the reduction uses at angle ``index`` where the list holds none, else ``None``.

        A compact list leaves the angles it does not give to the reducer
        (``_compact_reading``): the broadcast entry, or the reducer's own default.
        The editor shows that value, marked as implied, so a drop-down never says
        "unset" where the reduction has a value (``editor-combos`` C7). ``None``
        when the list holds an entry there; for a row past the reduction's count,
        which it never reads; for a list it reads as held, where an unset entry is
        a problem rather than a default; and for a derived λ, which is computed
        per run and is not known here.
        """
        field = fs.get(name)
        count = self.reduction_angles
        if not field.per_angle or not 0 <= index < count:
            return None
        current = self.get(name)
        if isinstance(current, (list, tuple)) and index < len(current) and current[index] is not None:
            return None
        reading = self._compact_reading(field, current, count)
        return None if reading is _NOT_COMPACT else reading

    def candidates(self, name, limit=MAX_CANDIDATES):
        """The file names a per-angle cell offers, as ``(names, total)``.

        The ``*.txt`` and ``*.dat`` files in the folder ``Field.candidates_folder``
        names (``DBname``: ``DBpath``, where the reducer reads each entry,
        ``nr_reduction_calc.py:402``), sorted, at most ``limit`` of them, with how
        many there were so a caller can say the list was cut. A sub-folder is not
        offered.

        Every failure is ``([], 0)``. The folder is on a facility mount, where a
        missing IPTS folder, a permission, a stale handle or an unresolvable path
        (an ``experiment_id`` of ``None`` makes the property itself raise) are
        ordinary, and the caller is a Qt slot. Nothing is cached: each call lists
        the folder the document resolves now, so the names follow the path. The
        caller asks when a cell offers candidates, never per keystroke or refresh.
        """
        folder_property = fs.get(name).candidates_folder
        if folder_property is None:
            return [], 0
        try:
            with os.scandir(getattr(self._config, folder_property)) as entries:
                names = sorted(
                    entry.name for entry in entries
                    if entry.name.lower().endswith(CANDIDATE_SUFFIXES) and _is_file(entry)
                )
        except (OSError, TypeError, ValueError):
            return [], 0
        return names[:limit], len(names)

    def angle_row(self, index):
        """Every per-angle value for one angle, as a dict."""
        if not 0 <= index < self.n_angles:
            raise IndexError(f"No angle at index {index} (have {self.n_angles})")
        row = {}
        for name in fs.PER_ANGLE_NAMES:
            current = self.get(name)
            row[name] = current[index] if isinstance(current, list) and index < len(current) else None
        return row

    # -- validation --------------------------------------------------------

    def validate(self):
        """Return a list of human-readable problems; empty means clean.

        Reports rather than raises: an editor has to show every problem at
        once, and a partly-filled document is a normal intermediate state, not
        an error. Unset (``None``) entries are therefore not flagged — except
        in an optional list, where a half-specified field is genuinely broken.
        """
        messages = []
        # Lengths are measured against the reduction's count, not the table's
        # rows: a list longer than that is surplus (notes()), and only a list
        # shorter than it is a problem (nr_reduction_calc.py:61-81).
        count = self.reduction_angles

        for field in fs.FIELD_SPEC:
            value = self.get(field.name)

            if field.per_angle:
                if value is None:
                    # A λ edit set_angle_field refused, reported while it still
                    # applies: the field is still derived and the row still surplus.
                    # Add and Remove move rows, and forget refusals; an accepted
                    # edit of the field forgets its own.
                    row = self._refused.get(field.name)
                    if row is not None and count <= row < self.n_angles:
                        messages.append(
                            f"{field.label} ({field.name}) is derived from the chopper ranges, so the "
                            f"value typed into surplus row {row + 1} was not kept: the reduction does "
                            f"not use that row, and a value there alone would stop it deriving the "
                            f"field at every angle; set it for each angle first, or leave it derived"
                        )
                    continue
                # A wrong TYPE is a problem to report, not a thing to iterate.
                # validate() used to walk straight into len()/enumerate() on
                # whatever a settings file happened to contain, so {"tof_min": 5}
                # raised TypeError out of a Qt slot and aborted the process.
                if not isinstance(value, (list, tuple)):
                    messages.append(field.check(value))
                    continue
                # "Unset at some angles" looks at the angles the reduction uses,
                # entries below the count. Past it is surplus that the reducer never
                # reads, and naming a surplus row as unset was a false problem
                # (review 1568397, B-1).
                angles = list(value[:count])
                missing = [i for i, entry in enumerate(angles) if entry is None]
                if field.optional_list and missing:
                    messages.append(
                        f"{field.label} ({field.name}) is set for some angles but not "
                        f"angles {missing}: either give every angle a value or clear "
                        f"the field to derive it from the chopper ranges"
                    )
                # A list the reducer fills or broadcasts: unset at EVERY angle is
                # written [] (its default applies, _encode_for_file); unset at
                # SOME angles cannot be written as a default and reaches the
                # reducer as None (:413 adds it; :82 calls .lower() on it).
                fills_itself = field.default_if_empty or field.broadcast_ok
                if fills_itself and missing and len(missing) < len(angles):
                    messages.append(
                        f"{field.label} ({field.name}) is set for some angles but not "
                        f"angles {missing}: either give every angle a value or clear "
                        f"the field to use the reduction's default"
                    )
                # The reducer lower-cases EVERY entry of a broadcast list (:82),
                # surplus rows included, and fails on an unset one there too —
                # unless the list is unset at every angle, and so written [].
                if field.broadcast_ok and len(missing) < len(angles):
                    unset = [i + 1 for i in range(count, len(value)) if value[i] is None]
                    if unset:
                        where = (f"row {unset[0]}" if len(unset) == 1
                                 else f"rows {', '.join(str(row) for row in unset)}")
                        messages.append(
                            f"{field.label} ({field.name}) has no value in surplus {where}: the "
                            f"reduction never uses a surplus row, but it reads every entry of this "
                            f"list and fails on an unset one; fill in a value or remove the row"
                        )
                length = 0 if fills_itself and len(missing) == len(angles) else len(value)
                if length < count and not self._length_is_allowed(field, length):
                    messages.append(
                        f"{field.label} ({field.name}) has {len(value)} entries "
                        f"for {count} angles"
                    )
                messages.extend(
                    field.check_element(entry, f" at angle {i}")
                    for i, entry in enumerate(value)
                    if entry is not None
                )
            elif field.runtime_owned:
                # The runtime record (LambdaMinUse/LambdaMaxUse) is not an
                # input. The reduction assigns it on every run
                # (nr_reduction_calc.py:385-391, in _load_and_extract_lambda,
                # before that method reads it at :452), and the editor shows it
                # read-only. A problem line here could not be acted on. The
                # declared list[float] disagrees with the writer's one scalar
                # per call; header-scale-factors-per-position decides the
                # record's shape, and until then this check stays out of it.
                continue
            else:
                messages.append(field.check(value))

        found = [m for m in messages if m]
        if len(found) > MAX_REPORTED_PROBLEMS:
            extra = len(found) - MAX_REPORTED_PROBLEMS
            found = found[:MAX_REPORTED_PROBLEMS]
            found.append(f"... and {extra} more problems not listed")
        return found

    def notes(self):
        """What the panel shows beside the problems: true of a file that reduces, worth knowing.

        ``validate()`` does not return these, so a reducible file still reads
        "No problems found."

        * A per-angle list longer than the reduction's count: the reducer
          indexes ``[i]`` with ``i < len(RBnum)`` (``nr_reduction_calc.py:61``)
          and never reads the extra entries. One line per field, never per
          entry; Remove angle on a surplus row drops them.
        * A boolean default list unset at every angle the reduction uses
          (``useBS`` — pinned by a test on the derivation): it is written ``[]``, which the reducer
          fills with 1, on, for every angle (``:102-103``). Shown only when
          angle-defining entries exist, because otherwise there is no reduction
          to describe.
        """
        count = self.reduction_angles
        defining = self._defining_length()
        lines = []
        for field in fs.FIELD_SPEC:
            value = self.get(field.name)
            if not field.per_angle or not isinstance(value, list):
                continue
            if len(value) > count:
                extra = len(value) - count
                lines.append(
                    f"{field.label} ({field.name}) has {extra} extra "
                    f"{'entry' if extra == 1 else 'entries'} beyond the {count} angles the "
                    f"reduction uses; it ignores them, and removing the surplus angle drops them"
                )
            if (field.default_if_empty and field.element_type == "bool" and defining
                    and all(entry is None for entry in value[:count])):
                lines.append(
                    f"{field.label} ({field.name}) is unset, so the reduction uses its default: "
                    f"on (1) at every angle (nr_reduction_calc.py:102-103)"
                )
        return lines

    @staticmethod
    def _length_is_allowed(field, length):
        """Is a per-angle length that differs from n_angles still legitimate?

        Three ways it can be, all from the reducer rather than from taste:
        a length-1 ``method_per_run`` is broadcast to every angle; the five
        ``default_if_empty`` arrays are filled in when left empty; and the
        runtime-owned ``RBnum`` comes from the runs being reduced, not from the
        author. Reporting these as problems is not harmless — a panel that
        cries wolf on a valid file teaches the scientist to ignore it, which is
        how a real problem goes unread.
        """
        if field.runtime_owned:
            return True
        if field.default_if_empty and length == 0:
            return True
        if field.broadcast_ok and length in (0, 1):
            return True
        return False

    # -- output ------------------------------------------------------------

    def normalize(self):
        """JSON-safe settings with the runtime-owned fields dropped.

        Those fields (``RBnum`` and the ``Lambda*Use`` record) are filled in by
        the reduction from the runs it is given; carrying an authored value for
        them would silently override the run. Integer-encoded lists are written
        as ``save()`` writes them (``_encode_for_file``).
        """
        return {
            key: value
            for key, value in self._encode_for_file(make_json_safe(self.to_dict()), self.reduction_angles).items()
            if key not in fs.RUNTIME_OWNED_NAMES
        }

    @staticmethod
    def _encode_for_file(values, count):
        """Write lists the way the reduction reads them.

        * A list the reducer fills or broadcasts (``default_if_empty``,
          ``broadcast_ok``) that is unset at every angle below ``count`` (the
          reduction's; entries past it are surplus it never reads) is written ``[]``. That
          is the only spelling of "use your default" the reducer has: it tests
          ``if not self.config.<name>`` (``nr_reduction_calc.py:99-110``), and
          ``[None, None]`` is a list, so the default is skipped and the ``None``
          used. Unset at some angles only is written as held; ``validate()``
          reports it. An angle-defining list is never collapsed, because unset
          there is not a default.
        * Integer-encoded lists are written as the reducer writes them, ``1``/``0``
          (``:103``). Which fields is declared (``Field.int_encoded``), not
          decided here by name. Only ``bool`` entries change; an unset entry
          stays ``None`` (``null``) and anything else is written as held.

        ``values`` is the fresh mapping ``make_json_safe`` built, never the
        document's own state.
        """
        for field in fs.FIELD_SPEC:
            if not (field.default_if_empty or field.broadcast_ok):
                continue
            entries = values.get(field.name)
            if isinstance(entries, list) and entries and all(entry is None for entry in entries[:count]):
                values[field.name] = []
        for name in fs.INT_ENCODED_NAMES:
            entries = values.get(name)
            if isinstance(entries, list):
                values[name] = [int(entry) if isinstance(entry, bool) else entry for entry in entries]
        return values

    def save(self, path):
        """Write the full document as a JSON settings file, atomically.

        The whole document, not ``normalize()``: this is the scientist's file
        and round-tripping it must not quietly drop fields. Use ``normalize()``
        when handing settings to a reduction.

        Written to a temporary file in the same directory, fsynced, then
        ``os.replace``d over the target. A plain ``open(path, "w")`` truncates
        before a single byte is produced, so anything that interrupts the dump —
        a full IPTS quota, a stalled ``/SNS`` mount, the process dying — leaves
        the previous good settings destroyed and a partial file in their place.
        The temp file shares the target's directory so the rename is atomic on
        that filesystem, and the fsync is not redundant: on NFS and FUSE mounts
        ``close()`` does not imply durability.

        Refuses to write through a symbolic link. The save dialog's overwrite
        confirmation names the link, not its target, so following one would
        overwrite a file the user never saw named.
        """
        path = Path(path)
        if path.is_symlink():
            raise ValueError(
                f"{path} is a symbolic link to {os.path.realpath(path)}; "
                f"refusing to write through it — save to the target directly if that is the intent"
            )
        payload = json.dumps(self._encode_for_file(make_json_safe(self.to_dict()), self.reduction_angles), indent=2)
        # mkstemp opens O_CREAT|O_EXCL on a fresh name, so there is no link to
        # follow and no pre-existing file to clobber.
        handle_fd, temporary = tempfile.mkstemp(
            dir=str(path.parent), prefix=path.name + ".", suffix=".tmp"
        )
        try:
            with os.fdopen(handle_fd, "w") as handle:
                handle.write(payload)
                handle.flush()
                os.fsync(handle.fileno())
            # Explicit, not umask: a settings file is meant to be readable by
            # collaborators on a shared IPTS directory. Deliberately not 0600.
            os.chmod(temporary, 0o644)
            os.replace(temporary, path)
        except BaseException:
            try:
                os.unlink(temporary)
            except OSError:
                pass
            raise
        return path

    def overrides(self):
        """The fields that differ from a fresh config — what this layer contributes.

        A resolution stack needs each layer's contribution, which is neither
        ``to_dict()`` (everything, so every layer would override every other)
        nor ``normalize()`` (everything minus the runtime-owned fields).
        """
        defaults = NRReductionConfig().__dict__
        current = self.to_dict()
        return {
            key: value
            for key, value in current.items()
            if key not in defaults or value != defaults[key]
        }

    def changed_vs_seed(self):
        """``{name: (seed_value, current_value)}`` for every field that moved."""
        current = self.to_dict()
        return {
            key: (self._seed[key], current[key])
            for key in current
            if key in self._seed and current[key] != self._seed[key]
        }
