"""Declarative description of every :class:`NRReductionConfig` field.

One table drives widget construction, the prompts the scientist reads,
validation messages, and the runtime-owned list used when normalizing a
document for saving. Without it each of those grows its own copy of "what
fields exist and what may they hold", and they drift.

**Qt-free on purpose.** This module and :mod:`lr_reduction.settings_document`
import nothing from ``qtpy``/``PyQt``. That seam is what lets the settings
model be tested in milliseconds without a display, and T3's resolution layer
builds on it.

**Names are storage names, not display names.** Every ``Field.name`` is exactly
a key of ``NRReductionConfig().__dict__`` — including the four private
``_*_override`` path fields. That is deliberate:

* ``json_to_config`` (``new_reduction_from_file.py:441``) gates on
  ``hasattr(config, key)`` and raises ``AttributeError`` on anything else, so a
  name that is not a real attribute is a hard failure at load, not a warning;
* the saved JSON is written from ``config.__dict__``
  (``new_reduction_from_file.py:145-151``), so ``__dict__`` keys are what a
  settings file actually contains;
* ``Spath``/``NEXUSpathRB``/``DBpath``/``BINpath`` are properties backed by
  those private names. Naming the public property instead would mean carrying a
  public-to-private mapping that can drift from the class — the exact failure
  this table exists to prevent.

The user-facing name lives in ``Field.label``.

``base_path`` is deliberately absent: it is a property with **no setter**, so
``hasattr`` passes and ``setattr`` raises. Any "mirror every attribute" loop
must exclude it.
"""

from dataclasses import dataclass
from typing import Any, Optional, Tuple

from lr_reduction.reduction_domains import (
    CALC_THETA_CHOICES,
    DET_RES_CHOICES,
    DET_RES_NOTES,
    METHOD_CHOICES,
    PEAK_TYPE_CHOICES,
)

# Re-exported from the single source, not mirrored here. These used to be
# hand-copies of bare local lists inside nr_reduction_calc — a copy waiting to
# drift, where the editor would go on offering a value the reducer had stopped
# accepting.
#
# What is and is not guaranteed, stated exactly, because an earlier version of
# this comment claimed "drift is structurally impossible" and that was an
# overclaim:
#
#   * the enforcing VALIDATORS derive from reduction_domains — nr_reduction_calc
#     builds valid_methods/valid_calc_theta from it, and nr_tools builds its
#     error messages from it. Those cannot drift.
#   * the DISPATCH does not. `if method == 'meantheta' ... elif 'constantq'` in
#     nr_reduction_calc, and the same shape in nr_tools, are literal by
#     necessity: each branch computes something different, so there is nothing
#     to derive them from. Adding a value to a domain here does NOT make the
#     reducer able to compute it.
#   * what closes that gap is the contents-equality pins in
#     tests/unit/lr_reduction/test_settings_document.py, which assert each
#     domain equals the exact set the dispatch handles, plus positive drivers
#     that call the consumers with every offered value.
#
# So: adding a value to a domain without teaching the dispatch is caught by a
# test, not prevented by construction. That is a weaker guarantee than the
# earlier text claimed, and it is the true one.
__all__ = [
    "CALC_THETA_CHOICES", "DET_RES_CHOICES", "METHOD_CHOICES", "PEAK_TYPE_CHOICES",
    "Field", "FIELD_SPEC", "BY_NAME", "PER_ANGLE_NAMES", "OPTIONAL_LIST_NAMES",
    "RUNTIME_OWNED_NAMES", "DEFAULT_IF_EMPTY_NAMES", "GROUPS", "TYPES",
    "get", "fields_in",
]


@dataclass(frozen=True)
class Field:
    """One configuration field.

    Attributes
    ----------
    name
        Exact ``NRReductionConfig.__dict__`` key. Pinned by a guard test.
    label
        User-facing name, shown in the editor.
    group
        Section the editor groups this field under.
    type
        Vocabulary term describing the value shape (see ``TYPES``).
    default
        The value a fresh ``NRReductionConfig`` carries.
    help
        The prompt shown to the scientist. Written to say what the field
        *does*, not to restate its name.
    allowed
        Permitted values, for enumerated fields. Empty means unconstrained.
    minimum, maximum
        Inclusive numeric bounds, where a bound is physically meaningful.
    per_angle
        True when the value is one entry per angle. ``add_angle`` grows every
        such field together.
    broadcast_ok
        Per-angle fields the reducer will broadcast from a single entry, so a
        length of 1 is valid rather than a mismatch.
    optional_list
        Per-angle fields initialised to ``None`` rather than ``[]``, where
        ``None`` means "derive it" and is a first-class state, not a missing
        value.
    runtime_owned
        Filled in by the reduction run, not authored by the scientist. Dropped
        by ``SettingsDocument.normalize()``.
    """

    name: str
    label: str
    group: str
    type: str
    default: Any
    help: str
    allowed: Tuple[Any, ...] = ()
    minimum: Optional[float] = None
    maximum: Optional[float] = None
    per_angle: bool = False
    broadcast_ok: bool = False
    optional_list: bool = False
    runtime_owned: bool = False
    default_if_empty: bool = False
    no_separators: bool = False
    falsy_means_off: bool = False
    case_sensitive: bool = False
    value_notes: Tuple[Tuple[str, str], ...] = ()

    # -- type vocabulary ------------------------------------------------

    @property
    def element_type(self):
        """Type of ONE entry of a list field; the field's own type otherwise.

        Strips a single level, so ``list[list[int]]`` yields ``list[int]`` —
        one ``BkgROI`` entry really is a list of pixel bounds.
        """
        if self.type.startswith("list[") and self.type.endswith("]"):
            return self.type[len("list[") : -1]
        return self.type

    @property
    def is_list(self):
        return self.type.startswith("list[")

    def default_value(self):
        """A copy of the default, never the shared object.

        ``frozen=True`` freezes the binding, not the list behind it. A caller
        that mutated ``field.default`` in place would corrupt this table for
        every consumer in the process — and a resolution stack *starts* from
        defaults, so the first such caller is likely rather than hypothetical.
        """
        if isinstance(self.default, list):
            return [list(e) if isinstance(e, list) else e for e in self.default]
        return self.default

    # -- text -> value ---------------------------------------------------

    def canonical(self, value):
        """Return the declared spelling of ``value``, if one matches case-insensitively."""
        if not self.allowed or not isinstance(value, str):
            return value
        for choice in self.allowed:
            if isinstance(choice, str) and choice.lower() == value.lower():
                return choice
        return value

    def coerce_element(self, text):
        """Coerce the text of ONE entry (a table cell) to this field's element type.

        An enumerated value is normalised to its declared spelling, so a user
        who types "Gaussian" stores "gaussian" — the spelling the reduction
        compares against — rather than a value that validates and then matches
        no branch.
        """
        return self.canonical(_coerce_typed(text, self.element_type))

    def coerce(self, text):
        """Coerce the text of a whole field value to its declared type.

        A list field parses a comma- or whitespace-separated sequence, which is
        what a single-line editor for ``data_x_range`` or
        ``emission_coefficients`` actually receives.
        """
        # One splitter, and one normalisation point. For a non-list field this
        # IS coerce_element, so it delegates rather than repeating the
        # canonical() call — a second copy is how the earlier coerce/check pair
        # drifted apart, and a mutation that removed one of them left the other
        # covering for it.
        if not self.is_list:
            return self.coerce_element(text)
        return self.canonical(_coerce_typed(text, self.type))

    # -- value -> problem ------------------------------------------------

    def check(self, value, where=""):
        """Return a problem with ``value``, or ``""``.

        ``None`` is "not filled in yet", which an editor has to be able to
        hold, so it is never reported here.
        """
        if value is None:
            return ""
        if self.is_list:
            if not isinstance(value, (list, tuple)):
                return (
                    f"{self.label} ({self.name}){where}: expected a list of "
                    f"{self.element_type}, got {type(value).__name__} {value!r}"
                )
            # Recurse. Checking only the container let a corrupted render round
            # trip silently: "[50, 200]" re-parsed as ['[50', '200]'] is a list,
            # so a container-only check called it clean while the reduction
            # received strings where it expected pixels.
            for index, entry in enumerate(value):
                problem = self.check_element(entry, f"{where}[{index}]")
                if problem:
                    return problem
            return ""
        return self.check_element(value, where)

    def check_element(self, value, where=""):
        """Return a problem with one entry, or with a scalar; or ``""``."""
        if value is None:
            return ""
        # Tri-state: a falsy value means "off", and only a truthy one has to be
        # in the allowed set. useCalcTheta is the case — the reducer skips the
        # whole block when it is falsy, and False is the class default, so a
        # fresh document must not report a problem.
        if self.falsy_means_off and not value:
            return ""
        if self.allowed:
            # Case matters where the CONSUMER compares exactly. nr_reduction_calc
            # lower-cases method_per_run and useCalcTheta before checking, so
            # those are genuinely case-insensitive; nr_tools compares DetResFn
            # and peak_type with ==, so "Gaussian" passes validation here and
            # then falls through every branch, and the reduction dies partway
            # with the settings file looking correct.
            if self.case_sensitive and value not in self.allowed:
                canonical = self.canonical(value)
                if canonical != value:
                    return (
                        f"{self.label} ({self.name}){where}: {value!r} differs in case from "
                        f"{canonical!r}, and this field is matched exactly by the reduction"
                    )
            if str(value).lower() not in {str(a).lower() for a in self.allowed}:
                for noted, reason in self.value_notes:
                    if str(value).lower() == noted.lower():
                        return f"{self.label} ({self.name}){where}: {value!r} — {reason}"
                return (
                    f"{self.label} ({self.name}){where}: {value!r} is not one of "
                    f"{', '.join(str(a) for a in self.allowed)}"
                )
            return ""
        expected = self.element_type
        problem = _type_problem(value, expected)
        if problem:
            return f"{self.label} ({self.name}){where}: {problem}"
        if expected == "path":
            problem = _path_problem(value)
            if problem:
                return f"{self.label} ({self.name}){where}: {problem}"
        if self.no_separators and isinstance(value, str):
            if "/" in value or "\\" in value:
                return (
                    f"{self.label} ({self.name}){where}: {value!r} contains a path separator; "
                    f"it names a single component, not a location"
                )
            if value.strip() in ("..", "."):
                return (
                    f"{self.label} ({self.name}){where}: {value!r} is a directory traversal, "
                    f"not a name"
                )
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            return ""
        if self.minimum is not None and value < self.minimum:
            return f"{self.label} ({self.name}){where}: {value} is below {self.minimum}"
        if self.maximum is not None and value > self.maximum:
            return f"{self.label} ({self.name}){where}: {value} is above {self.maximum}"
        return ""


#: Every type string a Field may declare. Asserted below, so a typo like
#: "flaot" fails at import rather than silently producing an unvalidated,
#: uncoerced text box — which is the premise the v1 coercion bug rested on.
TYPES = (
    "str", "int", "float", "bool", "path",
    "list[str]", "list[int]", "list[float]", "list[bool]", "list[list[int]]",
)

# Accepted spellings for a boolean in text. Never bool(text): bool("False") is
# True, which is how a scientist who turned background subtraction OFF got it
# applied anyway.
_TRUE = {"true", "1", "yes", "on", "t", "y"}
_FALSE = {"false", "0", "no", "off", "f", "n"}


def _coerce_typed(text, type_name):
    """Coerce one piece of text to ``type_name``. Empty means unset (``None``).

    A value that cannot be coerced is returned unchanged rather than raised on:
    the editor must be able to hold what the user typed, and ``check()`` is what
    reports it. Raising here would abort a Qt slot.
    """
    if not isinstance(text, str):
        return text
    stripped = text.strip()
    if stripped == "":
        return None
    if type_name == "bool":
        lowered = stripped.lower()
        if lowered in _TRUE:
            return True
        if lowered in _FALSE:
            return False
        return stripped
    if type_name == "int":
        try:
            return int(stripped)
        except ValueError:
            return stripped
    if type_name == "float":
        try:
            return float(stripped)
        except ValueError:
            return stripped
    if type_name.startswith("list["):
        inner = type_name[len("list[") : -1]
        parts = [p for p in stripped.replace(",", " ").split() if p]
        return [_coerce_typed(p, inner) for p in parts]
    return stripped


def _type_problem(value, type_name):
    """Describe how ``value`` contradicts ``type_name``, or return ``""``."""
    if type_name == "bool":
        return "" if isinstance(value, bool) else f"expected true/false, got {value!r}"
    if type_name == "int":
        # bool is an int subclass; a checkbox value in an int field is a bug.
        if isinstance(value, bool) or not isinstance(value, int):
            return f"expected a whole number, got {type(value).__name__} {value!r}"
        return ""
    if type_name == "float":
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            return f"expected a number, got {type(value).__name__} {value!r}"
        return ""
    if type_name in ("str", "path"):
        return "" if isinstance(value, str) else f"expected text, got {type(value).__name__} {value!r}"
    if type_name.startswith("list["):
        if not isinstance(value, (list, tuple)):
            return f"expected a list, got {type(value).__name__} {value!r}"
        inner = type_name[len("list[") : -1]
        for index, entry in enumerate(value):
            problem = _type_problem(entry, inner)
            if problem:
                return f"entry {index}: {problem}"
        return ""
    return ""


def _path_problem(value):
    """Reject a traversal in a field that genuinely holds a whole path.

    An earlier version also rejected absolute paths, which mis-modelled these
    fields: the four ``_*_override`` values are not joined to a base — they
    **are** the location (``nr_reduction_config``'s path properties). An
    absolute path is their normal shape, and it is exactly what
    ``QFileDialog.getExistingDirectory`` returns, so rejecting it flagged the
    only legitimate value while a relative one silently resolved against the
    process working directory instead.

    What stays rejected is ``..``, which is a traversal whether the path is
    absolute or relative. Fields that name a *component* rather than a path —
    ``experiment_id``, ``Sname`` and the subname siblings — carry
    ``no_separators`` instead, which is the stricter rule they need.
    """
    if ".." in value.replace("\\", "/").split("/"):
        return f"{value!r} contains '..', which points outside the intended directory"
    return ""

RUNS = "Runs and angles"
PATHS = "Paths"
NAMING = "Output naming"
PROCESSING = "Processing"
BACKGROUND = "Background"
QSPACE = "Q-space"
WAVELENGTH = "Wavelength and TOF"
THETA = "Theta and scaling"
GEOMETRY = "Instrument geometry"
DEADTIME = "Dead time"
RESOLUTION = "Detector resolution"
PEAK = "Peak fitting"
RUNTIME = "Runtime record"


FIELD_SPEC = (
    # ---- per-angle -------------------------------------------------------
    Field("method_per_run", "Q method", RUNS, "list[str]", [],
          "Lambda-to-Q conversion used for each angle. One entry per angle; a "
          "single entry is broadcast to all angles, and an empty list defaults "
          "to meanTheta.",
          allowed=METHOD_CHOICES, per_angle=True, broadcast_ok=True),
    Field("DBname", "Direct-beam file", RUNS, "list[str]", [],
          "Pre-processed direct-beam file backing each angle.", per_angle=True),
    Field("RBnum", "Run numbers", RUNS, "list[int]", [],
          "Run numbers reduced at each angle. Supplied by the reduction run, "
          "not authored here.",
          per_angle=True, runtime_owned=True),
    Field("RB_Ymin", "Peak Y min (pixel)", RUNS, "list[int]", [],
          "Lower edge of the specular peak window, in detector pixels.",
          minimum=0, per_angle=True),
    Field("RB_Ymax", "Peak Y max (pixel)", RUNS, "list[int]", [],
          "Upper edge of the specular peak window, in detector pixels.",
          minimum=0, per_angle=True),
    Field("BkgROI", "Background ROI", BACKGROUND, "list[list[int]]", [],
          "Background region per angle, as pixel bounds.", per_angle=True),
    Field("useBS", "Subtract background", BACKGROUND, "list[bool]", [],
          "Whether to subtract background at each angle.", per_angle=True, default_if_empty=True),
    Field("tof_min", "TOF min", WAVELENGTH, "list[float]", [],
          "Lower time-of-flight bound per angle.",
          minimum=0.0, per_angle=True, default_if_empty=True),
    Field("tof_max", "TOF max", WAVELENGTH, "list[float]", [],
          "Upper time-of-flight bound per angle.",
          minimum=0.0, per_angle=True, default_if_empty=True),
    Field("LambdaMin", "Lambda min", WAVELENGTH, "list[float]", None,
          "Lower wavelength bound per angle. Leave unset to derive it from the "
          "chopper ranges; if set, every angle needs a value.",
          per_angle=True, optional_list=True),
    Field("LambdaMax", "Lambda max", WAVELENGTH, "list[float]", None,
          "Upper wavelength bound per angle. Leave unset to derive it from the "
          "chopper ranges; if set, every angle needs a value.",
          per_angle=True, optional_list=True),
    Field("ThetaShift", "Theta shift (deg)", THETA, "list[float]", [],
          "Correction added to the measured theta at each angle.", per_angle=True, default_if_empty=True),
    Field("ScaleFactor", "Scale factor", THETA, "list[float]", [],
          "Multiplier applied to each angle's reflectivity before stitching.",
          per_angle=True, default_if_empty=True),

    # ---- scalars ---------------------------------------------------------
    Field("Sname", "Output name", NAMING, "str", "reduction_output",
          "Base name for the reduced output files.", no_separators=True),
    # A directory NAME, not a path: it is joined verbatim onto /SNS/REF_L, so
    # an absolute value replaces the base entirely and a '..' walks out of it.
    Field("experiment_id", "IPTS", NAMING, "str", "",
          "IPTS identifier. Also the root of every default path.",
          no_separators=True),
    Field("subname", "Output subtitle", NAMING, "str", None,
          "Optional subtitle appended to saved file names.", no_separators=True),
    Field("DTCsubname", "Dead-time-corrected suffix", NAMING, "str", "_DTC",
          "Suffix for dead-time-corrected outputs.", no_separators=True),
    Field("BINsubname", "Binned suffix", NAMING, "str", "_DTC",
          "Suffix for binned outputs.", no_separators=True),
    Field("errBINsubname", "Binned-error suffix", NAMING, "str", "_err_DTC",
          "Suffix for binned uncertainty outputs.", no_separators=True),
    Field("data_x_range", "Detector X range", RUNS, "list[int]", [50, 200],
          "Detector pixel range integrated over in X. Two values, not per angle."),

    Field("_Spath_override", "Output path", PATHS, "path", None,
          "Where reduced data is written. Unset uses <IPTS>/shared/reduced."),
    Field("_NEXUSpathRB_override", "NeXus path", PATHS, "path", None,
          "Where run NeXus files are read from. Unset uses <IPTS>/nexus."),
    Field("_DBpath_override", "Direct-beam path", PATHS, "path", None,
          "Where direct-beam files are read from. Unset uses "
          "<IPTS>/shared/transmission."),
    Field("_BINpath_override", "Binned-output path", PATHS, "path", None,
          "Where binned output is written. Unset uses <IPTS>/shared/reduced."),

    Field("Normalize", "Normalize to critical edge", PROCESSING, "bool", False,
          "Scale reflectivity to 1 over the critical-edge region set by Qnorm."),
    Field("AutoScale", "Auto-scale between angles", PROCESSING, "bool", False,
          "Scale each angle to its neighbour using the overlap region."),
    # NOT a bool, despite the name and the False default. The reducer accepts
    # 'detector_angle'/'sample_angle' and treats a legacy True as an alias for
    # the former (nr_reduction_calc, NRReduction.__init__). Declaring it bool
    # rendered a checkbox that could not express 'sample_angle' at all and
    # silently downgraded a loaded one on any toggle.
    Field("useCalcTheta", "Theta source", PROCESSING, "str", False,
          "Where theta comes from: the detector angle, or the fitted sample "
          "angle. Leave blank to keep the THS/THI log values.",
          allowed=CALC_THETA_CHOICES, falsy_means_off=True),
    Field("plotON", "Show plots", PROCESSING, "bool", True,
          "Display plots during reduction. Turn off for batch processing."),
    Field("plotQ4", "Plot as R*Q^4", PROCESSING, "bool", False,
          "Plot R*Q^4 instead of R."),
    Field("save8col", "Save 8-column output", PROCESSING, "bool", False,
          "Also write the 8-column form, adding L, dL, T and dT."),
    Field("useGravity", "Gravity correction", PROCESSING, "bool", True,
          "Apply the gravity correction to the neutron trajectory."),
    Field("use_emission_time", "Emission-time correction", PROCESSING, "bool", True,
          "Apply the moderator emission-time correction."),

    Field("qmin", "Q min", QSPACE, "float", 0.001,
          "Lower edge of the output Q range.", minimum=0.0),
    Field("qmax", "Q max", QSPACE, "float", 0.5,
          "Upper edge of the output Q range.", minimum=0.0),
    Field("dqbin", "Q bin width", QSPACE, "float", 0.005,
          "Width of the output Q bins.", minimum=0.0),
    Field("Qline_threshold", "Q-line threshold", QSPACE, "float", 1.0,
          "Fraction of a Q-line that must fall inside a bin for it to count, "
          "outside constantTOF mode.", minimum=0.0, maximum=1.0),
    Field("Qnorm", "Normalization Q", QSPACE, "float", 0.015,
          "Q below which data is treated as the critical-edge plateau when "
          "normalizing.", minimum=0.0),
    Field("tof_bin", "TOF bin width", WAVELENGTH, "float", 50,
          "Width of the time-of-flight bins.", minimum=0.0),

    Field("mmpix", "Pixel size (mm)", GEOMETRY, "float", None,
          "Detector pixel size. Unset reads it from the instrument settings."),
    Field("dSampDet", "Sample-detector distance", GEOMETRY, "float", None,
          "Unset reads it from the instrument settings."),
    Field("ny", "Vertical pixels", GEOMETRY, "int", None,
          "Number of pixels in Y. Unset reads it from the instrument settings."),
    # The source comments both ny and nx as "number of vertical pixels"; nx is
    # the horizontal count. Described correctly here rather than copying the
    # slip into the scientist-facing prompt.
    Field("nx", "Horizontal pixels", GEOMETRY, "int", None,
          "Number of pixels in X. Unset reads it from the instrument settings."),
    Field("dMod", "Moderator-detector distance", GEOMETRY, "float", None,
          "Unset reads it from the instrument settings."),
    Field("xi_ref", "xi reference distance", GEOMETRY, "float", None,
          "Distance defining xi = 0. Unset reads it from the instrument settings."),
    Field("dS1Samp", "S1-sample distance", GEOMETRY, "float", None,
          "Unset reads it from the instrument settings."),
    Field("IncidentTheta", "Incident theta (deg)", GEOMETRY, "float", None,
          "Beamline angle relative to earth, positive downwards. Unset reads "
          "the PV, falling back to 4.0 for older runs."),
    Field("emission_coefficients", "Emission-time coefficients", GEOMETRY,
          "list[float]", None,
          "Coefficients of the TOF emission-time correction."),

    Field("dead_time", "Dead time (us)", DEADTIME, "float", 4.2,
          "Detector dead time.", minimum=0.0),
    Field("dead_time_tof_step", "Dead-time TOF step", DEADTIME, "float", 50,
          "TOF bin width used when computing the dead-time correction.",
          minimum=0.0),

    Field("DetResFn", "Resolution function", RESOLUTION, "str", "rectangular",
          "Shape of the detector resolution function.", allowed=DET_RES_CHOICES,
          value_notes=DET_RES_NOTES, case_sensitive=True),
    Field("DetSigma", "Resolution sigma", RESOLUTION, "float", 0.8,
          "Width of the detector resolution function.", minimum=0.0),

    Field("peak_pad", "Peak fit padding (pixels)", PEAK, "int", 1,
          "Extra pixels included outside the background range when fitting the "
          "peak.", minimum=0),
    Field("peak_type", "Peak shape", PEAK, "str", "supergauss",
          "Function fitted to the specular peak.", allowed=PEAK_TYPE_CHOICES, case_sensitive=True),

    Field("LambdaMinUse", "Lambda min used", RUNTIME, "list[float]", None,
          "Wavelength bound the run actually used. Recorded by the reduction.",
          runtime_owned=True),
    Field("LambdaMaxUse", "Lambda max used", RUNTIME, "list[float]", None,
          "Wavelength bound the run actually used. Recorded by the reduction.",
          runtime_owned=True),
)


BY_NAME = {f.name: f for f in FIELD_SPEC}

# Asserted at import. The three FIELD_SPEC<->NRReductionConfig guards live in
# the tests; these two cannot wait for a test run, because an unknown type
# string silently degrades a field to an uncoerced, unvalidated text box.
_unknown_types = sorted({f.type for f in FIELD_SPEC} - set(TYPES))
assert not _unknown_types, f"FIELD_SPEC declares unknown type(s): {_unknown_types}"
assert len(BY_NAME) == len(FIELD_SPEC), "duplicate field name in FIELD_SPEC"

#: Storage names of every per-angle field, in FIELD_SPEC order. ``add_angle``
#: grows all of them together; a field missing from here is a field that
#: silently ends up a different length from its siblings.
PER_ANGLE_NAMES = tuple(f.name for f in FIELD_SPEC if f.per_angle)

#: Per-angle fields whose default is ``None`` rather than ``[]``.
OPTIONAL_LIST_NAMES = tuple(f.name for f in FIELD_SPEC if f.optional_list)

#: Fields the reduction run fills in, dropped when normalizing for save.
RUNTIME_OWNED_NAMES = tuple(f.name for f in FIELD_SPEC if f.runtime_owned)

#: Per-angle fields the reducer fills in itself when left empty
#: (nr_reduction_calc, "Set defaults for optional arrays"). An empty one is a
#: deliberate "use the default", not a length mismatch to report — reporting it
#: trains the scientist to ignore the panel, which is how a real problem hides.
DEFAULT_IF_EMPTY_NAMES = tuple(f.name for f in FIELD_SPEC if f.default_if_empty)

#: Groups in the order the editor should present them.
GROUPS = tuple(dict.fromkeys(f.group for f in FIELD_SPEC))


def get(name):
    """Return the :class:`Field` named ``name``."""
    return BY_NAME[name]


def fields_in(group):
    """Return the fields belonging to ``group``, in table order."""
    return tuple(f for f in FIELD_SPEC if f.group == group)
