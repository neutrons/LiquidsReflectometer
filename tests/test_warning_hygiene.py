"""Warning hygiene — keep third-party noise out so first-party signal shows.

matplotlib 3.9.4's `_mathtext` / `_fontconfig_pattern` still call pyparsing's
deprecated camelCase API, and pyparsing 3.3.2 warns about it. That is
third-party and there is no matplotlib release that removes it, so the suite
filters the pyparsing category — narrowly, because a bare category ignore would
also hide our own deprecations.

These two tests are the guard: the first fails if the filter stops covering the
noise, the second fails if it is ever widened enough to hide us.
"""

import warnings

import pytest


def _reset_warning_registries():
    """Let already-emitted warnings fire again.

    Python records emitted warnings per module and suppresses repeats, so a
    warning triggered by an earlier test would not reappear here and the check
    would pass vacuously. Note this deliberately does NOT use
    `simplefilter("always")`: that replaces the filter list, which would discard
    the very ini configuration under test and make the assertion unfalsifiable.
    """
    import sys

    for module in list(sys.modules.values()):
        registry = getattr(module, "__warningregistry__", None)
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
