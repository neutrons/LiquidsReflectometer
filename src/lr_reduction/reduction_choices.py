"""Canonical value domains for the reduction configuration.

One definition per domain, imported by both the code that *enforces* it and the
code that *offers* it to a user. Before this module the settings editor's
choice lists hand-mirrored bare local lists inside
:mod:`lr_reduction.nr_reduction_calc`, which is a copy waiting to drift: the
editor would go on offering a value the reducer had stopped accepting, and
nothing would notice until a reduction failed.

Deliberately dependency-free — stdlib only, no numpy, matplotlib or Qt — so the
settings model can import it without dragging the reduction stack in.

Spellings here are the canonical, human-facing ones. The reducer lower-cases
before comparing (``nr_reduction_calc.py`` in ``NRReduction.__init__``), so
``lowered()`` produces the form it matches against.
"""


def lowered(choices):
    """Return ``choices`` lower-cased, the form the reducer compares against."""
    return [choice.lower() for choice in choices]


#: Lambda-to-Q conversion, per angle (``NRReductionConfig.method_per_run``).
METHOD_CHOICES = ("meanTheta", "constantQ", "constantTOF")

#: The subset of :data:`METHOD_CHOICES` that
#: ``nr_reduction_calc._calculate_theta_and_bins`` dispatches on.
#:
#: ``constantTOF`` is validated by ``_validate_config`` but has no branch there —
#: it takes a different route through the reduction. Recorded as its own tuple
#: so the theta dispatch's error message can derive from what that function
#: actually handles, rather than from the full domain, which would name a method
#: it cannot compute.
THETA_DISPATCH_CHOICES = ("constantQ", "meanTheta")

#: Source of the theta value (``NRReductionConfig.useCalcTheta``).
#:
#: NOT a boolean, despite its name and its ``False`` default: the reducer
#: accepts these two strings and treats a legacy ``True`` as an alias for
#: ``detector_angle``. A settings editor that renders it as a checkbox cannot
#: express ``sample_angle`` at all, and silently downgrades a loaded one.
CALC_THETA_CHOICES = ("detector_angle", "sample_angle")

#: Detector resolution function (``NRReductionConfig.DetResFn``).
#:
#: Only these two are safe. The two consumers **disagree** about a third value,
#: and this module records that rather than silently picking a side:
#:
#: * ``nr_tools.calc_beam_on_detector`` accepts ``'none'``/``None`` and simply
#:   skips the convolution;
#: * ``nr_reduction_calc._calc_detector_convolution`` binds ``pad`` only under
#:   ``rectangular`` and ``gaussian``, so ``'none'`` reaches
#:   ``max(verts[:,1]) + pad`` with ``pad`` unbound and raises
#:   ``UnboundLocalError``.
#:
#: So ``'none'`` is not offered by the editor, and a settings file carrying it
#: is reported with the reason above rather than a bare "not one of" — it is a
#: real hazard in one code path, not merely an unlisted spelling. Whichever way
#: the inconsistency is resolved belongs in one of the two consumers, not here.
DET_RES_CHOICES = ("rectangular", "gaussian")

#: Values a consumer tolerates that this domain deliberately does not offer,
#: with the reason. Keyed by field-declared value.
#: What nr_tools tolerates in addition to DET_RES_CHOICES. Kept separate so
#: the editor's offer and the reducer's acceptance are not conflated.
DET_RES_TOLERATED = ("none",)

DET_RES_NOTES = (
    (
        "none",
        "nr_tools skips the convolution for 'none', but "
        "nr_reduction_calc._calc_detector_convolution raises UnboundLocalError "
        "on it (pad is only bound for 'rectangular'/'gaussian') — the two "
        "consumers disagree, so this value is unsafe",
    ),
)

#: Specular peak shape (``NRReductionConfig.peak_type``), dispatched in
#: ``nr_tools.fit_peak``.
PEAK_TYPE_CHOICES = ("gauss", "supergauss")
