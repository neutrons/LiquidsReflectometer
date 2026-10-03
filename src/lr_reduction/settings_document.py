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
itself. An edit (``set_angle_field``) changes exactly the entry edited.

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

    # -- construction ------------------------------------------------------

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
        """Rewrite values the reducer accepts only as legacy aliases.

        ``useCalcTheta = True`` is the case: the reducer maps it to
        ``detector_angle`` and then works normally
        (``nr_reduction_calc``, ``NRReduction.__init__``). Reporting it as a
        problem would cry wolf on a file that reduces perfectly well; leaving it
        alone would keep re-saving the deprecated spelling. Migrating it on load
        does what the reducer would have done, so the panel stays quiet and the
        file the scientist saves is explicit.
        """
        for field in fs.FIELD_SPEC:
            if not (field.falsy_means_off and field.allowed):
                continue
            if getattr(config, field.name, None) is True:
                setattr(config, field.name, field.allowed[0])

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
        * a compact list given a value is expanded: a broadcast list by
          repeating its single entry for the existing angles, the others with
          unset entries;
        * any other list gets the new entry at index ``reduction_angles``: a
          shorter list is first padded with unset entries to that index, and a
          longer one's surplus entries follow the new entry;
        * a per-angle value that is not a list is left as it is, for
          ``validate()`` to report.
        """
        m = self.reduction_angles
        for name in fs.PER_ANGLE_NAMES:
            field = fs.get(name)
            current = self.get(name)
            supplied = name in values
            if current is None:
                if supplied:
                    self.set(name, [None] * m + [values[name]])
                continue
            if not isinstance(current, (list, tuple)):
                continue
            current = list(current)
            compact = (field.default_if_empty and not current) or (field.broadcast_ok and len(current) <= 1)
            if compact:
                if supplied:
                    head = current * m if (field.broadcast_ok and len(current) == 1) else [None] * m
                    self.set(name, head + [values[name]])
                continue
            angles = current[:m] + [None] * (m - len(current[:m]))
            self.set(name, angles + [values.get(name)] + current[m:])

    def remove_angle(self, index):
        """Remove one angle from every per-angle field."""
        if not 0 <= index < self.n_angles:
            raise IndexError(f"No angle at index {index} (have {self.n_angles})")
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
        """
        field = fs.get(name)
        if not field.per_angle:
            raise KeyError(f"{name} is not a per-angle field")
        current = self.get(name)
        if not 0 <= index < self.n_angles:
            raise IndexError(f"No angle at index {index} (have {self.n_angles})")
        # A single broadcast entry stands for the reduction's angles
        # (nr_reduction_calc.py:77-79). Setting one of them expands it to the
        # reduction's count first, so the others keep that value instead of
        # becoming unset entries the reducer would choke on. It expands to that
        # count, not to the table's rows, because rows past it are surplus.
        if field.broadcast_ok and isinstance(current, (list, tuple)) and len(current) == 1:
            current = list(current) * max(self.reduction_angles, 1)
        # An edit changes exactly the entry edited. A list is padded only far
        # enough to reach the edited index, never to the table's row count: v1
        # padded to n_angles, and on a file with surplus rows one cell edit grew
        # every edited list into them, which raised the reduction's count (3 -> 7
        # on a real file) and reported problems that did not exist (review 1568397).
        # Padding a SHORT list rather than indexing past it matters too: an
        # IndexError in a Qt slot reaches qFatal() and kills the launcher, and short
        # lists are ordinary (a length-1 method_per_run; normalize() drops RBnum).
        # isinstance, not len(): a bare string has a length, and list("abc") would
        # explode it into ['a', 'b', 'c'] instead of replacing a wrong type.
        updated = list(current) if isinstance(current, (list, tuple)) else []
        updated += [None] * (index + 1 - len(updated))
        updated[index] = value
        field = fs.get(name)
        # An optional list that is emptied of every value goes back to None —
        # "derive it from the chopper ranges". Without this, touching one Lambda
        # cell is a one-way door out of that state for the life of the document.
        if field.optional_list and all(entry is None for entry in updated):
            updated = None
        self.set(name, updated)

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
