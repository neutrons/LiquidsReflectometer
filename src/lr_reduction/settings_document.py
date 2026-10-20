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

*Angles grow together.* ``add_angle`` mutates **every** per-angle field in one
operation. Growing only the obvious ones leaves the others short, and a short
array shifts every subsequent angle's settings by one — silently, since nothing
in the config class enforces equal lengths.

*Validation is not the same as the equal-length invariant.* The reducer
deliberately broadcasts a single ``method_per_run`` entry across all angles
(``nr_reduction_calc.py:76-78``) and defaults an empty one to ``meanTheta``
(``:41-42``). A validator that demanded strict equal lengths everywhere would
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
        return cls(config)

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
        """Number of angles, taken as the longest per-angle field."""
        lengths = [
            len(self.get(name))
            for name in fs.PER_ANGLE_NAMES
            if isinstance(self.get(name), list)
        ]
        return max(lengths) if lengths else 0

    def add_angle(self, **values):
        """Append one angle, growing every per-angle field together.

        Unsupplied entries are ``None`` — "not filled in yet", which an editor
        must be able to represent. The exception is the optional lists
        (``LambdaMin``/``LambdaMax``): while they are ``None`` the whole field
        means "derive it from the chopper ranges"
        (``nr_reduction_calc.py:381-383``), which is a valid configuration, not
        a missing one. Materialising them into ``[None, None]`` on the first add
        would turn that into a list that reports a length but carries no values
        — and ``web_report.py:547`` indexes it. So they are left alone until a
        value is actually supplied, at which point the list is created at full
        length and ``validate()`` reports the angles still lacking a value.
        """
        n = self.n_angles
        for name in fs.PER_ANGLE_NAMES:
            current = self.get(name)
            if current is None:
                if name not in values:
                    continue
                self.set(name, [None] * n + [values[name]])
            else:
                self.set(name, list(current) + [values.get(name)])

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
        # Pad a None field AND a SHORT one. Only the None case was handled
        # before, so any per-angle column shorter than n_angles raised
        # IndexError here — and an unhandled exception in a Qt slot calls
        # qFatal(), killing the whole launcher. Short columns arise from
        # ordinary files: the reducer sanctions a length-1 method_per_run, and
        # normalize() itself drops the runtime-owned RBnum, so the editor's own
        # save/reload round trip produces one.
        # isinstance, not len(): a per-angle field holding a bare string has a
        # length, so a len() guard let list("abc") explode it into
        # ['a','b','c'] instead of treating it as the wrong type it is.
        if not isinstance(current, (list, tuple)) or len(current) < self.n_angles:
            padded = [None] * self.n_angles
            if isinstance(current, (list, tuple)):
                padded[: len(current)] = list(current)
            current = padded
        if not 0 <= index < self.n_angles:
            raise IndexError(f"No angle at index {index} (have {self.n_angles})")
        updated = list(current)
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
        n = self.n_angles

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
                if field.optional_list and any(entry is None for entry in value):
                    missing = [i for i, entry in enumerate(value) if entry is None]
                    messages.append(
                        f"{field.label} ({field.name}) is set for some angles but not "
                        f"angles {missing}: either give every angle a value or clear "
                        f"the field to derive it from the chopper ranges"
                    )
                if len(value) != n and not self._length_is_allowed(field, len(value)):
                    messages.append(
                        f"{field.label} ({field.name}) has {len(value)} entries "
                        f"for {n} angles"
                    )
                messages.extend(
                    field.check_element(entry, f" at angle {i}")
                    for i, entry in enumerate(value)
                    if entry is not None
                )
            else:
                messages.append(field.check(value))

        found = [m for m in messages if m]
        if len(found) > MAX_REPORTED_PROBLEMS:
            extra = len(found) - MAX_REPORTED_PROBLEMS
            found = found[:MAX_REPORTED_PROBLEMS]
            found.append(f"... and {extra} more problems not listed")
        return found

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
        them would silently override the run.
        """
        return {
            key: value
            for key, value in make_json_safe(self.to_dict()).items()
            if key not in fs.RUNTIME_OWNED_NAMES
        }

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
        payload = json.dumps(make_json_safe(self.to_dict()), indent=2)
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
