"""Convenient, normalized access to a Mantid workspace's sample logs."""

from collections.abc import Callable, Collection, Mapping, Sequence
from typing import Any

import numpy as np
from mantid.api import Run, Workspace
from mantid.kernel import Property, StringTimeSeriesProperty

from lr_reduction.exceptions import (
    AmbiguousLogError,
    LogNotFoundError,
    LogTypeError,
    LogUnitError,
)
from lr_reduction.types import MantidWorkspace
from lr_reduction.utils.workspace import workspace_handle

#: Python types a Mantid sample log may hold as a plain scalar, as opposed to a series.
_SCALAR_TYPES = (bool, int, float, str)


class SampleLogs:
    """Convenient, normalized access to a workspace's sample logs.

    Sample logs are the single source of truth for run metadata, so this is the intended
    way to read any of it. Wraps either a workspace object or a workspace name::

        logs = SampleLogs(workspace)
        sequence_id = logs["sequence_id"]
        incidence_angle = logs.single_value("ths", operation=np.mean)

    Reading with ``[...]``
    ---------------------
    ``logs[name]`` returns a plain Python value, but only when the log has an
    unambiguous one: a scalar, or a time series whose entries are all equal. A series
    that genuinely varies raises `AmbiguousLogError` rather than silently picking an
    entry, because there is no single correct choice — the first and last values of a
    motor angle that settles during a run are different numbers, and `ths` feeds the
    incidence angle and therefore Q.

    So reserve ``[...]`` for logs that are constant across a run (identity,
    configuration, mode flags), and state the reduction explicitly for anything
    physics-bearing, via `single_value` or `time_average`.

    Note that this makes ``[...]`` deliberately non-total: ``name in logs`` being True
    does **not** guarantee ``logs[name]`` succeeds, since the log may be present but
    varying. Use `single_value` when you do not know a log's shape ahead of time.

    Values are normalized to Python natives — a log read never leaks a `numpy` scalar,
    which matters because `numpy.int32` is not an `int` and is not JSON serializable.
    The write side normalizes in the same direction, since Mantid rejects numpy scalars.

    Resolution is deferred, not cached
    ----------------------------------
    The constructor argument is kept as given and re-resolved on every access. A wrapper
    built from a *name* must behave this way: any algorithm writing its output back to
    the same name replaces the entry in the analysis data service with a new workspace
    object, and a wrapper holding the old object would silently read stale logs — or
    report a present log as missing.
    """

    def __init__(self, workspace: MantidWorkspace) -> None:
        """
        Parameters
        ----------
        workspace
            The workspace to read logs from, either the object or its name in the
            analysis data service.
        """
        self._workspace_ref = workspace

    def __contains__(self, name: str) -> bool:
        """Whether a log of this name exists. See the class docstring: this does not
        imply ``self[name]`` will succeed."""
        return self._run().hasProperty(name)

    def __getitem__(self, name: str) -> bool | int | float | str:
        """The log's single unambiguous value.

        Raises
        ------
        LogNotFoundError
            No log of this name exists.
        AmbiguousLogError
            The log is a series whose entries are not all equal, or records no values.
        """
        return self._normalize(name, self._property(name).value)

    # NOTE: this method shadows the builtin `property` for the remainder of the class
    # body, so no `@property` can be declared below it. The name is deliberate: it is
    # carried over from `legacy/mantid_utils.py::SampleLogValues`.
    def property(self, name: str) -> Property:
        """The raw Mantid `Property`, for the full time series, its units, or statistics.

        The escape hatch from this class's normalization; prefer the other accessors.
        """
        return self._property(name)

    def mean(self, name: str) -> float:
        """The **arithmetic** mean of a numeric log's values.

        Not time-weighted — each recorded entry counts equally regardless of how long it
        held. Use `time_average` for the time-weighted value, which is usually what a
        motor or beam log means.

        Equivalent to ``single_value(name, np.mean)``, but reduces the values this class
        has already validated rather than reading the log a second time. Deliberately not
        Mantid's statistics, which report `nan` for a vector-valued log — one written from
        a plain sequence, as `insert` does — instead of averaging its entries.

        Raises
        ------
        LogNotFoundError
            No log of this name exists.
        LogTypeError
            The log holds strings, so there is nothing to average.
        AmbiguousLogError
            The log records no values, or contains NaN entries, so no meaningful mean
            exists. Reduce it explicitly with ``single_value(name, np.nanmean)`` to skip
            the NaN entries.
        """
        _, values = self._numeric_values(name)
        return self._to_python(np.mean(values))

    def time_average(self, name: str) -> float:
        """The **time-weighted** mean of a numeric log's values.

        Weights each entry by how long it held, so a value recorded briefly contributes
        proportionally less. Contrast `mean`, which is the unweighted arithmetic mean.

        Defined only for a scalar or a time series. A vector-valued log records its
        entries without times, so there is no duration to weight them by; Mantid reports
        `nan` for such a log rather than saying so, so this rejects it outright.

        Raises
        ------
        LogNotFoundError
            No log of this name exists.
        LogTypeError
            The log holds strings, or is vector-valued and so carries no times to weight
            by.
        AmbiguousLogError
            The log records no values, or contains NaN entries, so no meaningful mean
            exists. Reduce it explicitly with ``single_value(name, np.nanmean)`` to skip
            the NaN entries.
        """
        prop, values = self._numeric_values(name)
        if self._is_vector(prop, values):
            raise LogTypeError(
                f"Sample log {name!r} is vector-valued and records no times, so it has no time-weighted "
                f"mean; read it with mean({name!r}) or single_value({name!r}, operation=...) instead"
            )
        return self._run().getStatistics(name).time_mean

    def single_value(self, name: str, operation: Callable[[Sequence], float] = np.mean) -> Any:
        """Reduce a log to one value, stating the reduction explicitly.

        Handles any log shape, so the caller need not know it ahead of time: a scalar is
        returned as-is, a numeric series is passed through `operation`, and a
        non-numeric series must be unambiguous (as for ``[...]``).

        The reduction is **unweighted** — `operation` sees the raw values with no regard
        for how long each held. Use `time_average` for the time-weighted mean.

        `operation` also owns the NaN policy, and is the only accessor here that does:
        the default `np.mean` propagates NaN, so a series with a NaN entry reduces to NaN.
        Pass `np.nanmean` (or another NaN-aware reduction) to skip those entries. The
        fixed-reduction accessors — `mean`, `time_average`, `find_log_with_units` — reject
        a NaN-bearing log instead, since they have no operation through which a caller
        could state the intent.

        Parameters
        ----------
        name
            Name of the sample log.
        operation
            Reduction applied to a numeric series; ignored for scalars.

        Raises
        ------
        LogNotFoundError
            No log of this name exists.
        AmbiguousLogError
            The log records no values, or is a non-numeric series whose entries are not
            all equal — no reduction function is meaningful over those.
        """
        value = self._property(name).value
        if isinstance(value, _SCALAR_TYPES):
            return value
        values = np.asarray(value)
        if values.size > 0 and np.issubdtype(values.dtype, np.number):
            return self._to_python(operation(values))
        return self._normalize(name, value)

    def find_log_with_units(
        self,
        name: str,
        unit: str | Collection[str] | None = None,
        operation: Callable[[Sequence], float] | None = None,
    ) -> float:
        """A numeric log's value, asserting its recorded unit is one the caller expects.

        Guards against a silently-wrong unit corrupting a geometry-derived quantity — a
        distance logged in the wrong unit would mis-scale Q with no other symptom.

        Parameters
        ----------
        name
            Name of the sample log.
        unit
            A unit, or several acceptable spellings of one. Pass a collection when the
            instrument is inconsistent: a single file records `ths` in ``"degree"`` but
            `thi` and `tthd` in ``"deg"``, so ``("deg", "degree")`` accepts both. `None`
            skips the check.
        operation
            Reduction for a varying series. When omitted the log must have a single
            unambiguous value, consistent with ``[...]``.

        Raises
        ------
        LogUnitError
            The recorded unit is not among those accepted.
        LogNotFoundError
            No log of this name exists.
        LogTypeError
            The log holds strings.
        AmbiguousLogError
            `operation` was omitted and the log has no single unambiguous value, or the
            log contains NaN entries. Note that the NaN check applies even when
            `operation` is NaN-aware, because the unit assertion runs first; read such a
            log with ``single_value(name, np.nanmean)`` and check the unit separately.
        """
        prop, _ = self._numeric_values(name)
        if unit is not None:
            accepted = (unit,) if isinstance(unit, str) else tuple(unit)
            if prop.units not in accepted:
                expected = " or ".join(repr(candidate) for candidate in accepted)
                raise LogUnitError(f"Sample log {name!r} is recorded in {prop.units!r}, not {expected}")
        if operation is None:
            return self._normalize(name, prop.value)
        return self.single_value(name, operation=operation)

    def insert(self, name: str, value: Any, unit: str | None = None) -> None:
        """Write a sample log onto the wrapped workspace, replacing any log of that name.

        The write goes straight onto the workspace's `Run` rather than through
        `AddSampleLog`. The trade-off is deliberate and costs provenance: an
        `AddSampleLog` call would show up in the workspace's algorithm history and this
        does not. It is made because `AddSampleLog` cannot express what this must write.

        - No vector form. `LogText` is a string, so a sequence can only go in as the
          *text* ``"[11111, 22222]"``; `Run.addProperty` records a real
          `VectorIntPropertyWithValue` that `mean` and `single_value` can reduce.
        - No boolean form, only ``LogText="1"``, which reads back as the integer 1 rather
          than `True`.
        - It fails outright on a workspace that is not in the analysis data service — as
          the intermediate workspaces inside an algorithm are not — and fails *after*
          adding the property, so the log is written and an exception raised.

        Writing to an intermediate workspace is not wasted effort even though the
        workspace itself is discarded: Mantid copies the `Run` from an algorithm's input
        to its output, so a log written early travels down the chain onto the result that
        does get saved.

        Numpy scalars and numpy array elements are converted to Python natives first:
        Mantid registers no converter for them and rejects them outright.

        Parameters
        ----------
        name
            Name of the sample log to write.
        value
            A scalar, or a list, tuple or array to record as a vector-valued log.
        unit
            Unit to record alongside the value.

        Raises
        ------
        LogTypeError
            `value` is not a shape Mantid can record: an empty sequence, a nested
            sequence, or anything that is neither a scalar nor a flat sequence.
        """
        if isinstance(value, np.generic):
            value = value.item()
        elif not isinstance(value, _SCALAR_TYPES):
            value = self._insertable_sequence(name, value)
        try:
            self._run().addProperty(name, value, unit or "", True)
        except (RuntimeError, TypeError, ValueError) as exc:
            # Mantid rejects a shape its converters do not cover — a nested sequence, for
            # one — with a bare RuntimeError or ValueError from deep inside its mapping
            # layer. Translate it so a caller can catch `SampleLogsError` and mean it.
            # The type, not the value: a long sequence would put a wall of numbers in the
            # message, and Mantid's own text already says what it could not convert.
            raise LogTypeError(f"Sample log {name!r} cannot be written from {type(value).__name__}: {exc}") from exc

    def _insertable_sequence(self, name: str, value: Any) -> list:
        """`value` as a list of Python natives, or a raise if it is not a flat sequence.

        Validates the shape *before* iterating. Iteration is what previously decided this,
        which meant an unsupported value failed as a bare `TypeError` ("'NoneType' object
        is not iterable") rather than the documented `LogTypeError` — `None` from an unset
        optional being the likely one. A mapping was worse than an error: iterating it
        yields its keys, so the log was silently written from those.
        """
        if isinstance(value, Mapping):
            raise LogTypeError(
                f"Sample log {name!r} cannot be written from a mapping; iterating one would "
                f"record its keys and silently discard the values"
            )
        if not isinstance(value, (Sequence, np.ndarray)):
            raise LogTypeError(
                f"Sample log {name!r} cannot be written from {type(value).__name__}; "
                f"pass a scalar, or a list, tuple or array of them"
            )
        if isinstance(value, np.ndarray) and value.ndim != 1:
            raise LogTypeError(
                f"Sample log {name!r} cannot be written from a {value.ndim}-dimensional array; "
                f"Mantid records a log as a scalar or a flat vector"
            )
        entries = [self._to_python(entry) for entry in value]
        if not entries:
            # Mantid raises a bare RuntimeError from deep inside its mapping layer
            # ("Cannot have a sequence type of length zero"); translate it like the
            # rest of this class translates Mantid's errors.
            raise LogTypeError(f"Sample log {name!r} cannot be written from an empty sequence")
        return entries

    def _workspace(self) -> Workspace:
        """The wrapped workspace, resolved fresh — see the class docstring."""
        return workspace_handle(self._workspace_ref)

    def _run(self) -> Run:
        """The wrapped workspace's `Run`, resolved fresh — see the class docstring."""
        return self._workspace().getRun()

    def _property(self, name: str) -> Property:
        """The named `Property`, or `LogNotFoundError`.

        Mantid raises a bare `RuntimeError` for an unknown log; every read here goes
        through this so callers get the package's own exception family instead.
        """
        run = self._run()
        if not run.hasProperty(name):
            raise LogNotFoundError(f"No sample log named {name!r} on this workspace")
        return run.getProperty(name)

    def _numeric_values(self, name: str) -> tuple[Property, np.ndarray]:
        """The named `Property` and its values, or a raise if there is no number to reduce.

        Every fixed-reduction accessor goes through this, so all three reject the same
        shapes. Mantid returns `nan` rather than raising for the statistics of a string
        log, an empty series or a NaN-bearing one, each of which would otherwise propagate
        a silently meaningless number into Q.

        Returns the values alongside the property because `prop.value` builds a fresh
        array on every access: a long run's `proton_charge` is large, and callers need both
        the array and the checks derived from it. Reading it once here is the difference
        between one copy and three per call.
        """
        prop = self._property(name)
        values = np.asarray(prop.value)

        # Type comes from the property, not from `values`: `np.asarray([])` is float64
        # whatever the log holds, so an empty string series would pass a dtype check and
        # be misreported as merely empty.
        if self._is_string(prop):
            raise LogTypeError(f"Sample log {name!r} holds strings, so it cannot be reduced to a number")
        if values.size == 0:
            raise AmbiguousLogError(f"Sample log {name!r} records no values, so it cannot be reduced to a number")
        # `np.bool_` is not a subtype of `np.number`, so booleans need naming explicitly.
        # They are reducible: the mean of a flag is the fraction of entries it held, and
        # `time_average` gives its duty cycle. Rejecting them would refuse a log `insert`
        # itself writes.
        if not (np.issubdtype(values.dtype, np.number) or np.issubdtype(values.dtype, np.bool_)):
            raise LogTypeError(f"Sample log {name!r} does not hold numeric values, so it cannot be reduced")
        if np.issubdtype(values.dtype, np.floating) and bool(np.isnan(values).any()):
            raise AmbiguousLogError(
                f"Sample log {name!r} contains NaN entries, so it has no meaningful mean; "
                f"reduce it explicitly with single_value({name!r}, operation=np.nanmean) "
                f"to skip them"
            )
        return prop, values

    @staticmethod
    def _is_string(prop: Property) -> bool:
        """Whether the property holds strings, determined without reading its values.

        A series' `type` is a mangled C++ template name, so the series case goes by class
        instead; the scalar and vector cases Mantid names plainly.
        """
        return prop.type in ("string", "str list") or isinstance(prop, StringTimeSeriesProperty)

    @staticmethod
    def _is_vector(prop: Property, values: np.ndarray) -> bool:
        """Whether the property holds a bare vector of values — neither scalar nor series.

        This pair of tests separates the three shapes a log can take: a scalar reduces to
        a zero-dimensional array, a time series carries `times`, and a vector log — what
        `insert` writes from any sequence — carries neither. Takes the already-materialized
        `values` so `prop.value` is not built a second time.
        """
        return values.ndim > 0 and not hasattr(prop, "times")

    def _normalize(self, name: str, value: Any) -> Any:
        """Reduce a raw log value to one Python scalar, or raise if it is ambiguous."""
        if isinstance(value, _SCALAR_TYPES):
            return value

        # np.asarray matters: a string series arrives as a Python list, for which
        # `value == value[0]` is False outright rather than an elementwise comparison.
        values = np.asarray(value)
        if values.size == 0:
            raise AmbiguousLogError(f"Sample log {name!r} records no values, so it has no single value to read")

        first = values.reshape(-1)[0]
        if bool(np.all(values == first)):
            return self._to_python(first)

        numeric = np.issubdtype(values.dtype, np.number)
        if numeric and bool(np.isnan(values).any()):
            # Deliberately this class and not a kind of its own: such a log is not
            # strictly ambiguous, but NaN != NaN means the entries cannot be compared at
            # all, and returning NaN from `[...]` would mis-scale Q with no other symptom.
            raise AmbiguousLogError(
                f"Sample log {name!r} contains NaN entries, so its values cannot be compared; "
                f"read it with single_value({name!r}, operation=...) to state the reduction"
            )
        if numeric:
            spread = f"{values.size} values, {values.min()} to {values.max()}"
        else:
            spread = f"{values.size} values, {len(set(values.tolist()))} distinct"
        raise AmbiguousLogError(
            f"Sample log {name!r} varies across the run ({spread}), so it has no single value; "
            f"read it with single_value({name!r}, operation=...) or time_average({name!r}) "
            f"to state the reduction"
        )

    @staticmethod
    def _to_python(value: Any) -> Any:
        """Convert a numpy scalar to its Python equivalent, leaving other values alone."""
        return value.item() if isinstance(value, np.generic) else value
