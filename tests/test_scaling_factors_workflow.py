import os
import pathlib

import mantid.simpleapi as mtd_api
import numpy as np
import pytest

mtd_api.config["default.facility"] = "SNS"
mtd_api.config["default.instrument"] = "REF_L"

from lr_reduction.scaling_factors import workflow as sf_workflow
from lr_reduction.utils import amend_config

# Per-field comparison bars, `|calc - ref| <= atol + rtol * |ref|`.
#
# Metadata fields are copied verbatim from the run logs and must round-trip;
# 1e-12 only absorbs a one-ulp difference between Mantid builds. Fitted fields
# reproduce bit-exactly on one build and drift ~1e-11 across builds; the bars
# are >= ~2000x that drift and ~6 orders inside the fit's own uncertainty
# (error_a/a is 1.7-4.4%). Widen only from a new measurement, never to make a
# red suite green.
#
# `b` gets 1e-7: a converged `b` moves 8.7e-11 when only Mantid's minimizer
# changes (Levenberg-MarquardtMD -> Levenberg-Marquardt); see
# scripts/test/measure_fit_path_dependence.py. The atol covers rows whose
# slope sits within 1 sigma of zero.
#
# Ceiling: the _46_200/_46_300 references differ only by TOF binning (`b` by
# 6.7e-04, `a` by 5.6e-05). Widening BOTH `b` past 6.7e-04 AND _FITTED past
# 5.6e-05 makes test_reference_files_are_pairwise_distinguishable vacuous.
_METADATA = (1e-12, 0.0)
_FITTED = (1e-8, 0.0)
_TOL = {
    "LambdaRequested": _METADATA,
    "S1H": _METADATA,
    "S2iH": _METADATA,
    "S1W": _METADATA,
    "S2iW": _METADATA,
    "a": _FITTED,
    "error_a": _FITTED,
    "error_b": _FITTED,
    "b": (1e-7, 1e-12),
}


def _parse_cfg(path):
    """One {field: value-string} dict per data row; duplicate keys are an error."""
    rows = []
    with open(path, "r") as fd:
        for lineno, line in enumerate(fd, 1):
            if line.startswith("#") or not line.strip():
                continue
            pairs = [tok.split("=", 1) for tok in line.split() if "=" in tok]
            keys = [k for k, _ in pairs]
            duplicates = sorted({k for k in keys if keys.count(k) > 1})
            assert not duplicates, f"{path} line {lineno}: duplicate field(s) {duplicates}"
            rows.append(dict(pairs))
    return rows


def check_results(data_file, reference):
    """
    Compare every field of a scaling-factor file against its reference.

    Numeric fields use the per-field bar in `_TOL`; non-numeric fields must
    match exactly. Row count and field order are checked too — the consumer
    (`template.py`) reads this format positionally. Returns the number of
    comparisons made, so a caller can assert the helper did the work.
    """
    cfg_data = _parse_cfg(data_file)
    cfg_ref = _parse_cfg(reference)

    assert len(cfg_data) == len(cfg_ref), (
        f"{data_file} has {len(cfg_data)} data rows, reference {reference} has {len(cfg_ref)}"
    )

    comparisons = 0
    for i, (row, ref) in enumerate(zip(cfg_data, cfg_ref)):
        assert list(row) == list(ref), (
            f"row {i}: field order differs from the reference; "
            f"got {list(row)} expected {list(ref)}"
        )
        for key, ref_str in ref.items():
            value_str = row[key]
            try:
                v_ref = float(ref_str)
            except ValueError:
                assert value_str == ref_str, f"row {i} {key}: {value_str!r} != reference {ref_str!r}"
                comparisons += 1
                continue
            try:
                v_calc = float(value_str)
            except ValueError:
                raise AssertionError(
                    f"row {i} {key}: {value_str!r} is not numeric but the reference {ref_str!r} is"
                )
            rtol, atol = _TOL[key]
            delta = np.fabs(v_calc - v_ref)
            bar = atol + rtol * np.fabs(v_ref)
            assert delta <= bar, (
                f"row {i} {key}: {v_calc!r} vs reference {v_ref!r} "
                f"(|delta| {delta:.3e} > {bar:.3e} = {atol:.0e} + {rtol:.0e}*|ref|)"
            )
            comparisons += 1
    return comparisons


# Discovered, not listed, so a new reference cannot be left out of the
# pairwise guard. sf_201043_Si.cfg is a reduction input and is excluded.
_REFERENCE_CFGS = tuple(
    sorted(p.name for p in (pathlib.Path(__file__).parent / "data").glob("sf_197912_Si*.cfg"))
)
# An empty parameter set would be reported as "skipped", not as a failure.
assert len(_REFERENCE_CFGS) >= 4, (
    f"reference glob found {_REFERENCE_CFGS}; expected at least the four committed references"
)

_REF_ROWS = [
    "IncidentMedium=Si LambdaRequested=9.74 S1H=0.391 S2iH=0.25 S1W=20.005 S2iW=20.0 "
    "a=1.1051209892255538 b=-5.536686634970115e-07 error_a=0.046673303933524965 error_b=1.3017622676214164e-06",
    "IncidentMedium=Si LambdaRequested=7.043 S1H=0.39 S2iH=0.25 S1W=19.952 S2iW=19.95 "
    "a=6.624303827034047 b=1.9231042128909928e-05 error_a=0.23125978138549352 error_b=9.279372359246127e-06",
]


def _write_cfg(path, rows):
    with open(path, "w") as fd:
        fd.write("# y=a+bx\n#\n")
        fd.writelines(row + "\n" for row in rows)
    return str(path)


def _mutated(row, key, value):
    """Return `row` with `key` set to `value`, preserving field order."""
    toks = []
    for tok in row.split():
        k, _, v = tok.partition("=")
        toks.append(f"{k}={value}" if k == key else f"{k}={v}")
    return " ".join(toks)


# 2 rows x 10 fields, IncidentMedium included.
_EXPECTED_COMPARISONS = 20


def test_check_results_accepts_an_identical_file(tmp_path):
    """Positive control, asserting the comparison count rather than 'did not raise'."""
    ref = _write_cfg(tmp_path / "ref.cfg", _REF_ROWS)
    data = _write_cfg(tmp_path / "data.cfg", _REF_ROWS)
    assert check_results(data, ref) == _EXPECTED_COMPARISONS


def test_tol_covers_every_numeric_reference_field():
    """A field missing from _TOL is a loud error, not a loose default."""
    fields = list(_parse_cfg(os.path.join(os.path.dirname(__file__), "data", "sf_197912_Si_auto.cfg"))[0])
    string_fields = {"IncidentMedium"}
    assert set(_TOL) == set(fields) - string_fields
    assert set(_TOL).isdisjoint(string_fields)


@pytest.mark.parametrize(
    "left, right",
    [
        pytest.param(a, b, id=f"{a}--vs--{b}")
        for i, a in enumerate(_REFERENCE_CFGS)
        for b in _REFERENCE_CFGS[i + 1 :]
    ],
)
def test_reference_files_are_pairwise_distinguishable(template_dir, left, right):
    """No two references may be interchangeable (the _46_300 file was once a copy of _46_200)."""
    with pytest.raises(AssertionError):
        check_results(os.path.join(template_dir, left), os.path.join(template_dir, right))


def test_check_results_detects_a_field_reorder(tmp_path):
    """Field order is part of the format; template.py reads it positionally."""
    ref = _write_cfg(tmp_path / "ref.cfg", _REF_ROWS)
    toks = _REF_ROWS[1].split()
    swapped = " ".join([toks[0]] + [toks[2], toks[1]] + toks[3:])
    data = _write_cfg(tmp_path / "data.cfg", [_REF_ROWS[0], swapped])
    with pytest.raises(AssertionError):
        check_results(data, ref)


def test_check_results_rejects_a_duplicate_field(tmp_path):
    """A repeated key must not silently last-wins into the dict."""
    ref = _write_cfg(tmp_path / "ref.cfg", _REF_ROWS)
    # The duplicate carries the reference value, so only the duplicate guard can raise.
    data = _write_cfg(tmp_path / "data.cfg", [_REF_ROWS[0], _REF_ROWS[1] + " a=6.624303827034047"])
    with pytest.raises(AssertionError):
        check_results(data, ref)


@pytest.mark.parametrize(
    "key, value",
    [
        # Gross corruption; none of these was detected before the fix.
        pytest.param("a", "999999.0", id="a-the-scaling-factor-itself"),
        pytest.param("b", "999999.0", id="b-the-slope"),
        pytest.param("error_a", "999999.0", id="error_a-uncertainty-on-a"),
        pytest.param("LambdaRequested", "999.0", id="LambdaRequested-the-row-s-wavelength"),
        pytest.param("S1W", "999.0", id="S1W-slit-width"),
        # Non-numeric on either side of the comparison.
        pytest.param("IncidentMedium", "Air", id="IncidentMedium-non-numeric"),
        pytest.param("a", "n/a", id="a-non-numeric-under-a-numeric-reference"),
        # 1% on `a` passed the old 0.02 bar.
        pytest.param("a", "6.690546865304387", id="a-1pc-high-under-the-old-0.02-bar"),
        # 0.04 mm on a 20 mm slit: the stale-reference delta this fix exposed.
        pytest.param("S1W", "19.992", id="S1W-0.04mm-the-stale-reference-delta"),
    ],
)
def test_check_results_detects_mutation(tmp_path, key, value):
    ref = _write_cfg(tmp_path / "ref.cfg", _REF_ROWS)
    mutated = [_REF_ROWS[0], _mutated(_REF_ROWS[1], key, value)]
    data = _write_cfg(tmp_path / "data.cfg", mutated)
    with pytest.raises(AssertionError):
        check_results(data, ref)


def test_check_results_detects_a_missing_field(tmp_path):
    ref = _write_cfg(tmp_path / "ref.cfg", _REF_ROWS)
    dropped = " ".join(t for t in _REF_ROWS[1].split() if not t.startswith("error_a="))
    data = _write_cfg(tmp_path / "data.cfg", [_REF_ROWS[0], dropped])
    with pytest.raises(AssertionError):
        check_results(data, ref)


def test_check_results_detects_a_truncated_file(tmp_path):
    ref = _write_cfg(tmp_path / "ref.cfg", _REF_ROWS)
    data = _write_cfg(tmp_path / "data.cfg", _REF_ROWS[:1])
    with pytest.raises(AssertionError):
        check_results(data, ref)


def test_compute_sf(nexus_dir, template_dir, tmp_path):
    """
    Test the computation of scaling factors
    """
    with amend_config(data_dir=nexus_dir):
        ws = mtd_api.Load("REF_L_197912")

    output_dir = str(tmp_path)

    # We are passing the first run of the set. For the autoreduction,
    # we would be missing runs from the complete set so we will want to
    # wait for the whole set to be acquired.
    output = sf_workflow.process_scaling_factors(ws, output_dir, use_deadtime=False, wait=True, postfix="_test")
    assert output is False

    output_cfg = os.path.join(output_dir, "sf_197912_Si_test.cfg")
    if os.path.isfile(output_cfg):
        os.remove(output_cfg)

    output = sf_workflow.process_scaling_factors(ws, output_dir, use_deadtime=False, wait=False, postfix="_test")
    assert output is True

    check_results(output_cfg, os.path.join(template_dir, "sf_197912_Si_auto.cfg"))


def test_compute_sf_with_deadtime(nexus_dir, template_dir, tmp_path):
    """
    Test the computation of scaling factors
    """
    with amend_config(data_dir=nexus_dir):
        ws = mtd_api.Load("REF_L_197912")

    output_dir = str(tmp_path)

    output_cfg = os.path.join(output_dir, "sf_197912_Si_test_dt.cfg")
    if os.path.isfile(output_cfg):
        os.remove(output_cfg)

    output = sf_workflow.process_scaling_factors(ws, output_dir, use_deadtime=True, wait=False, postfix="_test_dt")
    assert output is True

    check_results(output_cfg, os.path.join(template_dir, "sf_197912_Si_dt_par_42_200.cfg"))


def test_compute_sf_with_deadtime_tof_300(nexus_dir, template_dir, tmp_path):
    """
    Test the computation of scaling factors
    """
    with amend_config(data_dir=nexus_dir):
        ws = mtd_api.Load("REF_L_197912")

    output_dir = str(tmp_path)

    output_cfg = os.path.join(output_dir, "sf_197912_Si_test_dt.cfg")
    if os.path.isfile(output_cfg):
        os.remove(output_cfg)

    output = sf_workflow.process_scaling_factors(
        ws,
        output_dir,
        use_deadtime=True,
        deadtime=4.6,
        deadtime_tof_step=300,
        paralyzable=False,
        wait=False,
        postfix="_test_dt",
    )
    assert output is True

    check_results(output_cfg, os.path.join(template_dir, "sf_197912_Si_dt_par_46_300.cfg"))


def test_compute_sf_with_deadtime_tof_200(nexus_dir, template_dir, tmp_path):
    """
    Test the computation of scaling factors
    """
    with amend_config(data_dir=nexus_dir):
        ws = mtd_api.Load("REF_L_197912")

    output_dir = str(tmp_path)

    output_cfg = os.path.join(output_dir, "sf_197912_Si_test_dt.cfg")
    if os.path.isfile(output_cfg):
        os.remove(output_cfg)

    output = sf_workflow.process_scaling_factors(
        ws,
        output_dir,
        use_deadtime=True,
        deadtime=4.6,
        deadtime_tof_step=200,
        paralyzable=False,
        wait=False,
        postfix="_test_dt",
    )
    assert output is True

    check_results(output_cfg, os.path.join(template_dir, "sf_197912_Si_dt_par_46_200.cfg"))


def test_compute_sf_with_deadtime_tof_200_sort(nexus_dir, template_dir, tmp_path):
    """
    Test the computation of scaling factors with order_by_runs=False.

    Compares against the same reference as its order_by_runs=True sibling
    because the two orderings coincide for this run set, so the flag is not
    actually exercised (pre-existing).
    """
    with amend_config(data_dir=nexus_dir):
        ws = mtd_api.Load("REF_L_197912")

    output_dir = str(tmp_path)

    output_cfg = os.path.join(output_dir, "sf_197912_Si_test_dt.cfg")
    if os.path.isfile(output_cfg):
        os.remove(output_cfg)

    output = sf_workflow.process_scaling_factors(
        ws,
        output_dir,
        order_by_runs=False,
        use_deadtime=True,
        deadtime=4.6,
        deadtime_tof_step=200,
        paralyzable=False,
        wait=False,
        postfix="_test_dt",
    )
    assert output is True

    check_results(output_cfg, os.path.join(template_dir, "sf_197912_Si_dt_par_46_200.cfg"))
