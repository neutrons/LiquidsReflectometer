"""Warning hygiene — keep third-party noise out so first-party signal shows.

matplotlib 3.9.4's `_mathtext` / `_fontconfig_pattern` still call pyparsing's
deprecated camelCase API, and pyparsing 3.3.2 warns about it. That is
third-party and there is no matplotlib release that removes it, so the suite
filters the pyparsing category — narrowly, because a bare category ignore would
also hide our own deprecations.

These two tests are the guard: the first fails if the filter stops covering the
noise, the second fails if it is ever widened enough to hide us.

The last two guard the helper that lets already-emitted warnings fire again:
it resets the registries without running any module's lookup hook.
"""

import sys
import types
import warnings

import pytest


def _reset_warning_registries():
    """Let already-emitted warnings fire again.

    Python records emitted warnings per module and suppresses repeats, so a
    warning triggered by an earlier test would not reappear here and the check
    would pass vacuously. Note this deliberately does NOT use
    `simplefilter("always")`: that replaces the filter list, which would discard
    the very ini configuration under test and make the assertion unfalsifiable.

    The registry is read from each module's own namespace, never with getattr:
    a getattr runs a module-level __getattr__ (PEP 562) or a lazy loader when
    the name is missing, and mpmath 1.4's deprecated `rational` and `math2`
    modules warn from theirs.
    """
    for module in list(sys.modules.values()):
        try:
            registry = object.__getattribute__(module, "__dict__").get("__warningregistry__")
        except (AttributeError, TypeError):
            # sys.modules may hold None (a blocked import) or an object with no namespace.
            continue
        if registry:
            registry.clear()


def _render_mathtext():
    """Render the kind of label the Overplot tab draws (`R·Q⁴`)."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    figure = Figure()
    FigureCanvasAgg(figure)
    axes = figure.add_subplot(111)
    axes.set_ylabel(r"$R \cdot Q^4$")
    figure.canvas.draw()


def test_no_pyparsing_warning_from_mathtext():
    """The suite's filterwarnings must actually cover the mathtext noise."""
    pyparsing_warnings = pytest.importorskip("pyparsing.warnings")
    category = pyparsing_warnings.PyparsingDeprecationWarning

    _reset_warning_registries()
    with warnings.catch_warnings(record=True) as caught:
        _render_mathtext()

    leaked = [w for w in caught if issubclass(w.category, category)]
    assert not leaked, (
        f"{len(leaked)} PyparsingDeprecationWarning(s) escaped the ini filter, "
        f"first from {leaked[0].filename}:{leaked[0].lineno}"
    )


def test_first_party_deprecation_still_visible():
    """The filter must not be widened into a bare category ignore.

    If someone replaces the pyparsing-scoped entry with
    `ignore::DeprecationWarning`, our own deprecations vanish with the noise —
    this goes red instead.
    """
    _reset_warning_registries()
    with warnings.catch_warnings(record=True) as caught:
        warnings.warn("first-party probe", DeprecationWarning, stacklevel=1)

    assert any(
        issubclass(w.category, DeprecationWarning) and "first-party probe" in str(w.message) for w in caught
    ), "a first-party DeprecationWarning is being suppressed; the warning filter is too broad"


def _warn_on_lookup(name):
    warnings.warn(f"a module-level __getattr__ ran for {name!r}", DeprecationWarning, stacklevel=2)
    raise AttributeError(name)


def test_resetting_the_registries_runs_no_module_hook(monkeypatch):
    """W4: the reset reads each module's own namespace, so a module-level __getattr__ (PEP 562) never runs.

    mpmath 1.4's deprecated `rational` and `math2` modules warn from theirs. sympy imports mpmath, and the suite
    imports sympy at collection (`lr_reduction.user_defined_function`), so this file's own reset put their two
    DeprecationWarnings in the suite's summary. The reset still clears a registry it finds.
    """
    hooked = types.ModuleType("hygiene_probe_hooked")
    hooked.__getattr__ = _warn_on_lookup
    registered = types.ModuleType("hygiene_probe_registered")
    registered.__warningregistry__ = {("probe", UserWarning, 1): True}
    # Out of sys.modules again before pytest reports: its traceback code reads __file__ from every module.
    with monkeypatch.context() as patch:
        patch.setitem(sys.modules, hooked.__name__, hooked)
        patch.setitem(sys.modules, registered.__name__, registered)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            _reset_warning_registries()
    assert registered.__warningregistry__ == {}


def test_resetting_the_registries_with_mpmath_loaded_warns_nothing():
    """W4, the reported case: with mpmath imported, the reset raises no DeprecationWarning."""
    pytest.importorskip("mpmath")
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        _reset_warning_registries()
