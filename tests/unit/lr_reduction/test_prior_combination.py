"""Re-reducing runs of a sequence must not change the new-workflow output headers.

Autoreduction calls ``reduce_from_file([run], ..., check_for_prior=True)`` once per
run; each call merges the new run with the files already in the output folder
(the "priors") and rewrites every file of the sequence.  The ``# Run Title:``,
``# Angles:`` and ``NR_runs`` headers are lists indexed by sequence position
(``seq_num - 1``), padded with ``null`` up to the highest position present, and
each position's entry comes from the same source as that position's data -- so
any number of reductions, in any order, gives the same per-run log headers
(``NR_runs``, ``Run Title``, ``Angles``) and data as reducing the sequence once.
``Scaling factors`` and ``Lambda Range`` are per position too (F17,
header-scale-factors-per-position): the factor applied to each position's data
and the range its reduction used; ``Config.ScaleFactor`` stays the authored list.

The inputs here are **synthetic**: ``NR_Reduction._reduce_single_run`` (the
physics) is replaced by a stub returning a made-up R(Q) = 1e-6 Q^-4 -- not a
measurement.  Everything else is the real code: sequence grouping from the NeXus
logs (tiny HDF5 files written per test), ``NR_Reduction.reduce()`` bookkeeping,
prior discovery and merge, and header writing.
"""

import ast
import json
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

import lr_reduction.new_reduction_from_file as nrff
from lr_reduction.nr_reduction_calc import NR_Reduction

EXPERIMENT = "IPTS-00000"
SEQ_ID = 221472
THI = -0.002
THCEN_OFFSET = 0.1  # the stub's ThCen differs from ths, so a ThCen read back from NeXus is detectable
RUNS = {  # run: (sequence number, ths, title) -- modelled on IPTS-36119 221472-221474
    221472: (1, -0.45, "Si Ir Air-221472-1."),
    221473: (2, -1.201, "Si Ir Air-221472-2."),
    221474: (3, -3.5, "Si Ir Air-221472-3."),
}
R1, R2, R3 = RUNS
N_POINTS = 40


def title(run):
    return RUNS[run][2]


def ths(run):
    return RUNS[run][1]


def thcen(run):
    return round(ths(run) + THCEN_OFFSET, 3)


def write_nexus(nexus_dir, run):
    seq, theta, name = RUNS[run]
    with h5py.File(nexus_dir / f"REF_L_{run}.nxs.h5", "w") as f:
        f["entry/title"] = [name.encode()]
        logs = f.create_group("entry/DASlogs")
        logs["BL4B:CS:Autoreduce:Sequence:Num/value"] = np.array([seq])
        logs["BL4B:CS:Autoreduce:Sequence:Id/value"] = np.array([SEQ_ID])
        logs["BL4B:Mot:ths.RBV/value"] = np.array([theta])
        logs["BL4B:Mot:thi.RBV/value"] = np.array([THI])


def synthetic_curve(position, reverse_q):
    """Overlapping Q bands, one per sequence position (reversed: position 1 has the highest Q).

    Each position is off by its own factor, so autoscaling has real work to do.
    """
    band = (2 - position) if reverse_q else position
    q = np.geomspace(0.008, 0.03, N_POINTS) * 1.8**band
    r = 1e-6 * q**-4 * 1.5**position
    return q, r, 0.05 * r, 0.02 * q


@pytest.fixture
def env(tmp_path, monkeypatch):
    nexus = tmp_path / "nexus"
    out = tmp_path / "autoreduce"
    nexus.mkdir()
    out.mkdir()
    for run in RUNS:
        write_nexus(nexus, run)
    settings = {
        "method_per_run": ["meanTheta"] * 3,
        "DBname": ["DB_A1.dat", "DB_A2.dat", "DB_A3.dat"],
        "RB_Ymin": [146, 142, 130],
        "RB_Ymax": [155, 159, 167],
        "Sname": "blank",
        "experiment_id": EXPERIMENT,
        "AutoScale": True,
        "LambdaMin": [2.7, 2.6, 2.5],
        "LambdaMax": [9.5, 9.5, 9.5],
        "ThetaShift": [0, 0, 0],
        "ScaleFactor": [1, 1, 1],
        "LambdaMinUse": 2.5,
        "LambdaMaxUse": 9.5,
    }
    settings_file = tmp_path / "reduce_settings.json"
    settings_file.write_text(json.dumps(settings))
    state = SimpleNamespace(nexus=nexus, out=out, settings=settings_file, reverse_q=False)

    def fake_reduce_single_run(self, i, rb_num, save=True, start_times=None, end_times=None):  # noqa: ARG001 -- the real method's signature; save is not used
        seq = RUNS[rb_num][0]
        assert i == seq - 1, "reduce() must pass the run's sequence position"
        # The real one's signature (time-slicing-reconcile): these reductions have no time window, so none arrives
        assert start_times is None and end_times is None, (start_times, end_times)
        q, r, dr, dq = synthetic_curve(i, state.reverse_q)
        zeros = np.zeros_like(q)
        # reduce() takes the title from self.log_values and the angles from the returned log_vals
        self.log_values = {"title": title(rb_num), "ths": ths(rb_num), "thi": THI, "ThCen": ths(rb_num) + THCEN_OFFSET}
        # As the real _reduce_single_run does (header-scale-factors-per-position): the call's wavelength range
        # (nr_reduction_calc.py:391-397) and the authored scale factor applied to the data (:1139-1140).
        self.config.LambdaMinUse = self.config.LambdaMin[i]
        self.config.LambdaMaxUse = self.config.LambdaMax[i]
        r, dr = r * self.config.ScaleFactor[i], dr * self.config.ScaleFactor[i]
        result = {"q": q, "r": r, "dr": dr, "dq": dq, "t": zeros, "l": zeros, "dt": zeros, "dl": zeros}
        return result, self.config, self.log_values

    monkeypatch.setattr(NR_Reduction, "_reduce_single_run", fake_reduce_single_run)
    return state


def reduce_runs(env, runs):
    """Reduce each run by itself, as autoreduction does, into the same output folder."""
    for run in runs:
        nrff.reduce_from_file([run], env.settings, EXPERIMENT, datapath=env.nexus, plot=False,
                              override_params={"Spath": env.out, "subname": "autoreduction"},
                              check_for_prior=True)


def read_outputs(out):
    """{file name: (header dict, data array)} for every .dat file in the output folder."""
    outputs = {}
    for path in sorted(Path(out).glob("*.dat")):
        header, records = {}, {}
        with open(path) as f:
            for line in f:
                if not line.startswith("#"):
                    break
                line = line[1:].strip()
                if line.startswith("NR_runs = "):
                    header["NR_runs"] = ast.literal_eval(line[len("NR_runs = "):])
                elif line.startswith("Run Title: "):
                    header["title"] = json.loads(line[len("Run Title: "):])["title"]
                elif line.startswith("Angles: "):
                    header.update(json.loads(line[len("Angles: "):]))
                elif line.startswith("Header format: "):
                    header["format"] = line
                elif line.startswith("Scaling factors = "):
                    records["scale"] = json.loads(line[len("Scaling factors = "):])["scale_factor"]
                elif line.startswith("Lambda Range = "):
                    try:
                        lam = json.loads(line[len("Lambda Range = "):])
                        records["lambda_min"], records["lambda_max"] = lam["lambda_min"], lam["lambda_max"]
                    except (ValueError, TypeError, KeyError):
                        records["lambda_min"] = records["lambda_max"] = line[len("Lambda Range = "):]
                elif line.startswith("Config: "):
                    config = json.loads(line[len("Config: "):])
                    records["config_scale"] = config.get("ScaleFactor")
                    records["config_lambda_min"] = config.get("LambdaMinUse")
                    records["config_lambda_max"] = config.get("LambdaMaxUse")
        header["records"] = records
        outputs[path.name] = (header, np.loadtxt(path, unpack=True))
    return outputs


def canonical_names():
    return {f"REFL_{SEQ_ID}_{RUNS[run][0]}_{run}_autoreduction.dat" for run in RUNS} | {
        f"REFL_{SEQ_ID}_combined_autoreduction.dat"}


def expected_header(present):
    """The header every file of the sequence carries once `present` runs have been reduced."""
    top = max(RUNS[run][0] for run in present)
    slot = {RUNS[run][0]: run for run in present}
    runs = [slot.get(p) for p in range(1, top + 1)]
    return {
        "NR_runs": runs,
        "title": [title(r) if r else None for r in runs],
        "THS": [ths(r) if r else None for r in runs],
        "THI": [THI if r else None for r in runs],
        "ThCen": [thcen(r) if r else None for r in runs],
    }


def headers_only(outputs):
    """The per-run log headers (M1's scope); the F17 records are compared by their own tests."""
    return {name: {k: v for k, v in header.items() if k not in ("format", "records")}
            for name, (header, _) in outputs.items()}


# ---------------------------------------------------------------------------
# One pass, in order: the reference every other scenario must reproduce


def test_in_order_headers_describe_each_position(env):
    reduce_runs(env, [R1, R2, R3])
    outputs = read_outputs(env.out)
    assert set(outputs) == canonical_names()
    for name, header in headers_only(outputs).items():
        assert header == expected_header(RUNS), name


def test_header_format_marker_precedes_config(env):
    reduce_runs(env, [R1])
    path = next(env.out.glob("*_1_*.dat"))
    lines = [line for line in path.read_text().splitlines() if line.startswith("#")]
    marker = [i for i, line in enumerate(lines) if line.startswith("# Header format: 3")]
    config = [i for i, line in enumerate(lines) if line.startswith("# Config: ")]
    assert marker and config and marker[0] < config[0]
    # load_from_file() stops at "# Config:" and must still find the settings
    assert nrff.load_from_file(path)["config"]["DBname"] == ["DB_A1.dat", "DB_A2.dat", "DB_A3.dat"]


# ---------------------------------------------------------------------------
# Re-reduction, order and gaps: identical headers and data to one clean pass

SCENARIOS = {
    "re-reduce all twice": [R1, R2, R3, R1, R2, R3],
    "re-reduce all three times": [R1, R2, R3] * 3,
    "re-reduce one run": [R1, R2, R3, R2],
    "out of order": [R3, R1, R2],
    "gap, then filled": [R1, R3, R2],
    "last run first, then re-reduce": [R3, R2, R1, R3],
}


@pytest.fixture
def canonical(env, tmp_path_factory):
    """Outputs of one in-order pass, reduced into a separate folder."""
    reference = SimpleNamespace(**{**vars(env), "out": tmp_path_factory.mktemp("canonical")})
    reduce_runs(reference, [R1, R2, R3])
    return read_outputs(reference.out)


@pytest.mark.parametrize("runs", SCENARIOS.values(), ids=SCENARIOS.keys())
def test_rereduction_reproduces_single_pass(env, canonical, runs):
    reduce_runs(env, runs)
    outputs = read_outputs(env.out)
    assert set(outputs) == set(canonical)
    assert headers_only(outputs) == headers_only(canonical)
    for name, (_, data) in outputs.items():
        np.testing.assert_allclose(data, canonical[name][1], rtol=1e-9, err_msg=name)


def test_batch_of_the_whole_sequence_matches_single_pass(env, canonical):
    # the GUI batch path: one call for all runs, so no priors are involved
    nrff.reduce_from_file([R1, R2, R3], env.settings, EXPERIMENT, datapath=env.nexus, plot=False,
                          override_params={"Spath": env.out, "subname": "autoreduction"},
                          check_for_prior=True)
    outputs = read_outputs(env.out)
    assert headers_only(outputs) == headers_only(canonical)
    for name, (_, data) in outputs.items():
        np.testing.assert_allclose(data, canonical[name][1], rtol=1e-9, err_msg=name)


def test_gap_pads_with_null_and_does_not_crash(env):
    reduce_runs(env, [R1, R3])
    outputs = read_outputs(env.out)
    assert set(outputs) == canonical_names() - {f"REFL_{SEQ_ID}_2_{R2}_autoreduction.dat"}
    for name, header in headers_only(outputs).items():
        assert header == expected_header([R1, R3]), name
    combined = outputs[f"REFL_{SEQ_ID}_combined_autoreduction.dat"][1]
    assert combined.shape[1] == 2 * N_POINTS  # both runs' data, not just the last one


def test_later_position_reduced_first_is_padded(env):
    reduce_runs(env, [R3])
    for name, header in headers_only(read_outputs(env.out)).items():
        assert header == expected_header([R3]), name


def test_file_names_and_headers_follow_sequence_not_q_order(env):
    env.reverse_q = True  # position 1 has the highest Q
    reduce_runs(env, [R1, R2, R3, R2])
    outputs = read_outputs(env.out)
    assert set(outputs) == canonical_names()
    for name, header in headers_only(outputs).items():
        assert header == expected_header(RUNS), name


@pytest.mark.parametrize("runs", [[R1, R2, R3], [R3, R1, R2, R1]], ids=["in order", "out of order"])
def test_eight_column_files_carry_the_same_headers(env, runs):
    for run in runs:
        nrff.reduce_from_file([run], env.settings, EXPERIMENT, datapath=env.nexus, plot=False,
                              override_params={"Spath": env.out, "subname": "autoreduction", "save8col": True},
                              check_for_prior=True)
    outputs = read_outputs(env.out)
    assert {name for name in outputs if name.endswith("_8col.dat")} == {
        name.replace(".dat", "_8col.dat") for name in canonical_names()}
    for name, header in headers_only(outputs).items():
        assert header == expected_header(RUNS), name


R2B = 221475  # a second run at sequence position 2: the step was measured again
# R2B's stub curve sits on a stretched Q grid, so a file shows WHOSE data it holds. A factor on R
# would not: AutoScale is on in these settings and absorbs a multiplicative difference.
R2B_Q_STRETCH = 1.01


@pytest.fixture
def remeasured(env, monkeypatch):
    monkeypatch.setitem(RUNS, R2B, (2, -1.3, "Si Ir Air-221472-2. again"))
    write_nexus(env.nexus, R2B)
    stub = NR_Reduction._reduce_single_run  # env's physics stub

    def remeasured_stub(self, i, rb_num, save=True, start_times=None, end_times=None):
        result, config, logs = stub(self, i, rb_num, save, start_times=start_times, end_times=end_times)
        if rb_num == R2B:
            result = {**result, "q": result["q"] * R2B_Q_STRETCH}
        return result, config, logs

    monkeypatch.setattr(NR_Reduction, "_reduce_single_run", remeasured_stub)
    return env


def q_grid(run):
    """The stub's Q grid for `run` (reverse_q off)."""
    q = synthetic_curve(RUNS[run][0] - 1, reverse_q=False)[0]
    return q * R2B_Q_STRETCH if run == R2B else q


def test_current_run_wins_at_its_position(remeasured):
    reduce_runs(remeasured, [R1, R2, R3, R2B])
    header = read_outputs(remeasured.out)[f"REFL_{SEQ_ID}_combined_autoreduction.dat"][0]
    assert header["NR_runs"] == [R1, R2B, R3]
    assert header["THS"][1] == ths(R2B)


def test_two_runs_at_one_position_are_reported(remeasured, capsys):
    reduce_runs(remeasured, [R1, R2, R3, R2B])
    capsys.readouterr()
    reduce_runs(remeasured, [R1])
    assert f"sequence position 2 has files for runs [{R2}, {R2B}]" in capsys.readouterr().out


def test_remeasured_position_survives_rereduction_of_another_run(remeasured):
    reduce_runs(remeasured, [R1, R2, R3, R2B, R1])
    header = read_outputs(remeasured.out)[f"REFL_{SEQ_ID}_combined_autoreduction.dat"][0]
    assert header["NR_runs"] == [R1, R2B, R3]


# ---------------------------------------------------------------------------
# Prior files: take a file's own entry only when its header vouches for it


def header(nr_runs, entries, marked):
    """A parsed prior header with `entries` = list of runs (or None) whose values fill the lists."""
    return {
        "format": 2 if marked else None,
        "NR_runs": nr_runs,
        "title": [title(r) if r else None for r in entries],
        "ths": [ths(r) if r else None for r in entries],
        "thi": [THI if r else None for r in entries],
        "ThCen": [thcen(r) if r else None for r in entries],
    }


def own(run):
    return {"title": title(run), "ths": ths(run), "thi": THI, "ThCen": thcen(run)}


@pytest.mark.parametrize(("prior", "seq", "run", "expected"), [
    pytest.param(header([R1, R2, R3], [R1, R2, R3], True), 2, R2, own(R2), id="marked, positional"),
    pytest.param(header([R1, None, R3], [R1, None, R3], True), 3, R3, own(R3), id="marked, with a gap"),
    pytest.param(header([R1, R2, R3], [R1, R2, R3], True), 2, R3, None, id="marked, wrong run at position"),
    pytest.param(header([R1, R2], [R1, R2, R3], True), 2, R2, None, id="marked, length mismatch"),
    # legacy: written by a single-run call, which appended its own entry last
    pytest.param(header([None, R2], [R1, R2], False), 2, R2, own(R2), id="legacy single-run, in order"),
    pytest.param(header([None, R2], [R3, R2], False), 2, R2, own(R2), id="legacy single-run, out of order"),
    pytest.param(header([R1], [R1], False), 1, R1, own(R1), id="legacy, first run"),
    pytest.param(header([None, None, R3], [R3], False), 3, R3, own(R3), id="legacy, first run at position 3"),
    pytest.param(header([None, None, R3], [R1, R2, R3] * 2, False), 3, R3, own(R3), id="legacy, duplicated"),
    # legacy files that cannot vouch for their own entry
    pytest.param(header([None, None, R3], [R1, R2, R3], False), 1, R1, None, id="legacy, rewritten by another run"),
    pytest.param(header([R1, None, R3], [R2, R1, R3], False), 1, R1, None, id="legacy multi-run batch"),
    pytest.param(header([None, R2], [R1, R2], False) | {"ths": [ths(R2)]}, 2, R2, None,
                 id="legacy, lists of different lengths"),
    pytest.param({"format": None, "NR_runs": None, "title": None, "ths": None, "thi": None, "ThCen": None},
                 1, R1, None, id="no header"),
])
def test_own_logs_from_header(prior, seq, run, expected):
    assert nrff.own_logs_from_header(prior, seq, run) == expected


def test_read_logs_from_nexus(env):
    assert nrff.read_logs_from_nexus(R2, env.nexus) == {"title": title(R2), "ths": ths(R2), "thi": THI, "ThCen": None}


def write_legacy_file(out, run, nr_runs, entries):
    """A partial file as the pre-fix code wrote it: dense, append-order lists, no format marker."""
    seq = RUNS[run][0]
    q, r, dr, dq = synthetic_curve(seq - 1, reverse_q=False)
    angles = {"THS": [ths(e) for e in entries], "THI": [THI for _ in entries], "ThCen": [thcen(e) for e in entries]}
    head = "\n".join([
        f"NR_runs = {nr_runs}",
        f"Run Title: {json.dumps({'title': [title(e) for e in entries]})}",
        f"Angles: {json.dumps(angles)}",
        "---" * 20,
        f"Config: {json.dumps({'RBnum': nr_runs})}",
        "---" * 20,
        "columns = Q, R, dR, dQ (sigma)",
        "---" * 20,
    ])
    np.savetxt(out / f"REFL_{SEQ_ID}_{seq}_{run}_autoreduction.dat", np.column_stack((q, r, dr, dq)),
               header=head, delimiter="\t")


def assert_headers_after_fallback(out, thcen_known):
    """Every header is the canonical one, except ThCen is null where it could not be recovered."""
    expected = expected_header(RUNS)
    expected["ThCen"] = [thcen(r) if r in thcen_known else None for r in (R1, R2, R3)]
    for name, header in headers_only(read_outputs(out)).items():
        assert header == expected, name


# The legacy folders below are what the pre-fix code (exp 9aaaefa) leaves behind: its merge
# rewrote every file of the sequence with the header of the call that reduced the last run.

def test_legacy_in_order_pass_then_rereduce(env, canonical):
    # pre-fix, in order: all three files carry the header written while reducing R3
    for run in (R1, R2, R3):
        write_legacy_file(env.out, run, [None, None, R3], [R1, R2, R3])

    reduce_runs(env, [R1])  # R3's file vouches for its own entry; R2's does not (falls back to NeXus)
    assert_headers_after_fallback(env.out, thcen_known={R1, R3})

    reduce_runs(env, [R2])
    assert headers_only(read_outputs(env.out)) == headers_only(canonical)


def test_legacy_multi_run_batch_is_not_trusted(env):
    # pre-fix: R2 alone, then a batch of R1 and R3 -> [t2] + [t1, t3] under NR_runs [R1, None, R3]
    for run in (R1, R2, R3):
        write_legacy_file(env.out, run, [R1, None, R3], [R2, R1, R3])

    reduce_runs(env, [R2])  # position 1 must not take t2, the first list entry
    assert_headers_after_fallback(env.out, thcen_known={R2})


def test_corrupted_legacy_priors_fall_back_to_nexus_then_heal(env, canonical):
    # a folder left by the pre-fix code after two passes: every list doubled
    for run in (R1, R3):
        write_legacy_file(env.out, run, [None, None, R3], [R1, R2, R3] * 2)

    reduce_runs(env, [R2])  # R3's file wrote its own entry last; R1's file was written by R3's call
    assert_headers_after_fallback(env.out, thcen_known={R2, R3})

    reduce_runs(env, [R1])  # re-reducing the stale position heals the sequence
    assert headers_only(read_outputs(env.out)) == headers_only(canonical)


# ---------------------------------------------------------------------------
# D-5: at a re-measured sequence position the highest run number wins, whatever the order of
# reduction, and the superseded run's file is left in place (charter D-5, [human, 2026-10-02])

COMBINED = f"REFL_{SEQ_ID}_combined_autoreduction.dat"


def combined_header(out):
    return read_outputs(out)[COMBINED][0]


def partial_name(seq, run):
    return f"REFL_{SEQ_ID}_{seq}_{run}_autoreduction.dat"


@pytest.mark.parametrize("runs", [
    pytest.param([R1, R2, R3, R2B, R2], id="superseded run re-reduced last"),
    pytest.param([R2B, R1, R2, R3], id="newer run reduced first"),
    pytest.param([R1, R2B, R3, R2, R1], id="superseded run reduced after it, then another"),
])
def test_highest_run_wins_in_any_order(remeasured, runs):
    reduce_runs(remeasured, runs)
    outputs = read_outputs(remeasured.out)
    header = outputs[COMBINED][0]
    assert header["NR_runs"] == [R1, R2B, R3]
    assert (header["title"][1], header["THS"][1]) == (title(R2B), ths(R2B))
    # the data of position 2 are R2B's too, not only its header entries
    np.testing.assert_allclose(outputs[partial_name(2, R2B)][1][0], q_grid(R2B), rtol=1e-9)


def test_superseded_file_is_left_in_place(remeasured):
    reduce_runs(remeasured, [R1, R2, R3, R2B, R1])
    assert (remeasured.out / partial_name(2, R2)).exists()


def test_shared_position_report_names_the_run_used(remeasured, capsys):
    reduce_runs(remeasured, [R1, R2, R3, R2B])
    capsys.readouterr()
    reduce_runs(remeasured, [R1])
    assert f"sequence position 2 has files for runs [{R2}, {R2B}]; using run {R2B}" in capsys.readouterr().out


def write_misnamed_copy(out, run, wrong_seq, nr_runs, entries):
    """What the pre-fix merge rewrite left when it named files by Q-order index: run `run`'s file
    under another position's name, beside the true one."""
    write_legacy_file(out, run, nr_runs, entries)
    (out / partial_name(RUNS[run][0], run)).rename(out / partial_name(wrong_seq, run))


def test_misnamed_legacy_copy_does_not_take_a_position(env, capsys):
    # the V10 folder: the true files of R1 and R2, plus each run's copy under the other's position
    write_misnamed_copy(env.out, R2, 1, [R1, R2], [R1, R2])
    write_misnamed_copy(env.out, R1, 2, [R1, R2], [R1, R2])
    for run in (R1, R2):
        write_legacy_file(env.out, run, [R1, R2], [R1, R2])
    reduce_runs(env, [R3])
    out = capsys.readouterr().out
    assert combined_header(env.out)["NR_runs"] == [R1, R2, R3]
    for copy in (partial_name(1, R2), partial_name(2, R1)):
        assert f"{copy} is not used" in out, copy


def test_a_misnamed_copy_does_not_displace_the_current_run(env, capsys):
    # R2's copy under position 1 has the higher run number than the run being reduced there
    write_misnamed_copy(env.out, R2, 1, [R1, R2, R3], [R1, R2, R3])
    for run in (R2, R3):
        write_legacy_file(env.out, run, [R1, R2, R3], [R1, R2, R3])
    reduce_runs(env, [R1])
    out = capsys.readouterr().out
    assert combined_header(env.out)["NR_runs"] == [R1, R2, R3]
    assert f"{partial_name(1, R2)} is not used" in out


def write_nexus_without_sequence(nexus_dir, run):
    seq, theta, name = RUNS[run]
    with h5py.File(nexus_dir / f"REF_L_{run}.nxs.h5", "w") as f:
        f["entry/title"] = [name.encode()]
        logs = f.create_group("entry/DASlogs")
        logs["BL4B:Mot:ths.RBV/value"] = np.array([theta])
        logs["BL4B:Mot:thi.RBV/value"] = np.array([THI])


@pytest.mark.parametrize("damage", ["file removed", "sequence log missing"])
def test_unreadable_nexus_keeps_the_candidate(remeasured, capsys, damage):
    reduce_runs(remeasured, [R1, R2, R3, R2B])
    (remeasured.nexus / f"REF_L_{R2B}.nxs.h5").unlink()
    if damage == "sequence log missing":
        write_nexus_without_sequence(remeasured.nexus, R2B)
    capsys.readouterr()
    reduce_runs(remeasured, [R1])  # position 2 is shared by R2 and R2B: both are checked
    out = capsys.readouterr().out
    assert combined_header(remeasured.out)["NR_runs"] == [R1, R2B, R3]
    assert f"{partial_name(2, R2B)} could not be checked" in out


def test_unshared_positions_do_not_consult_nexus(env, canonical, capsys):
    reduce_runs(env, [R1, R2, R3])
    for run in (R2, R3):
        (env.nexus / f"REF_L_{run}.nxs.h5").unlink()
    capsys.readouterr()
    reduce_runs(env, [R1])
    out = capsys.readouterr().out
    assert headers_only(read_outputs(env.out)) == headers_only(canonical)
    assert "could not be checked" not in out
    assert "is not used" not in out


@pytest.mark.parametrize("runs", [[R2B, R2, R3], [R2, R2B, R3]], ids=["newer run first", "newer run second"])
def test_two_runs_of_one_position_in_one_call(remeasured, capsys, runs):
    nrff.reduce_from_file(runs, remeasured.settings, EXPERIMENT, datapath=remeasured.nexus, plot=False,
                          override_params={"Spath": remeasured.out, "subname": "autoreduction"},
                          check_for_prior=True)
    out = capsys.readouterr().out
    outputs = read_outputs(remeasured.out)
    assert partial_name(2, R2B) in outputs
    assert partial_name(2, R2) not in outputs
    assert outputs[COMBINED][0]["NR_runs"] == [None, R2B, R3]
    # The whole report line: the run left out must be named in it. A bare `str(R2) in out` matched
    # the call's own "Beginning run set [...]" echo, so it held whatever the report said (review 4617053).
    assert f"this call names runs [{R2}, {R2B}] at sequence position 2; reducing run {R2B} only" in out


R2C = 221476  # a third run at sequence position 2


def test_three_runs_at_one_position(remeasured, monkeypatch, capsys):
    monkeypatch.setitem(RUNS, R2C, (2, -1.25, "Si Ir Air-221472-2. third"))
    write_nexus(remeasured.nexus, R2C)
    reduce_runs(remeasured, [R1, R2C, R2, R2B])
    capsys.readouterr()
    reduce_runs(remeasured, [R3])  # position 2 is then three prior files, none of them this call's
    assert combined_header(remeasured.out)["NR_runs"] == [R1, R2C, R3]
    assert (f"sequence position 2 has files for runs [{R2}, {R2B}, {R2C}]; using run {R2C}"
            in capsys.readouterr().out)


def test_run_numbers_compare_as_integers(env, monkeypatch):
    # name order puts "..._2_100000_..." before "..._2_99999_..."; text comparison ranks "99999" highest
    low, high = 99999, 100000
    for run, theta in ((low, -1.21), (high, -1.22)):
        monkeypatch.setitem(RUNS, run, (2, theta, f"position 2, run {run}"))
        write_nexus(env.nexus, run)
    reduce_runs(env, [R1, high, low])
    assert combined_header(env.out)["NR_runs"] == [R1, high]


def test_eight_column_set_follows_the_same_rule(remeasured):
    for run in [R1, R2, R3, R2B, R1]:
        nrff.reduce_from_file([run], remeasured.settings, EXPERIMENT, datapath=remeasured.nexus, plot=False,
                              override_params={"Spath": remeasured.out, "subname": "autoreduction", "save8col": True},
                              check_for_prior=True)
    outputs = read_outputs(remeasured.out)
    for name in (COMBINED, COMBINED.replace(".dat", "_8col.dat")):
        assert outputs[name][0]["NR_runs"] == [R1, R2B, R3], name
    np.testing.assert_allclose(outputs[partial_name(2, R2B).replace(".dat", "_8col.dat")][1][0],
                               q_grid(R2B), rtol=1e-9)


def test_a_position_whose_every_file_is_a_misnamed_copy_stays_a_gap(remeasured, capsys):
    """Added when mutation F7 (the guard for a position left with no candidate) survived the battery:
    nothing put a position where every file belongs elsewhere, and there the highest of nothing raises.
    """
    for run in (R2, R3):
        write_misnamed_copy(remeasured.out, run, 1, [R1, R2, R3], [R1, R2, R3])
        write_legacy_file(remeasured.out, run, [R1, R2, R3], [R1, R2, R3])
    reduce_runs(remeasured, [R2B])
    out = capsys.readouterr().out
    assert combined_header(remeasured.out)["NR_runs"] == [None, R2B, R3]
    for copy in (partial_name(1, R2), partial_name(1, R3)):
        assert f"{copy} is not used" in out, copy
    assert f"sequence position 1 has files for runs [{R2}, {R3}]; none of them belongs to it" in out


# ---------------------------------------------------------------------------
# header-scale-factors-per-position (F17): the scale factor applied to each position's data and the
# wavelength range its reduction used, per sequence position, the same in every file whatever the order

AUTHORED = [1, 2, 1]


def records(header):
    return header["records"]


@pytest.fixture
def batch(env, tmp_path_factory):
    """The whole sequence in one call (the GUI batch): ground truth for the applied factors, no priors."""
    reference = SimpleNamespace(**{**vars(env), "out": tmp_path_factory.mktemp("batch")})
    nrff.reduce_from_file([R1, R2, R3], env.settings, EXPERIMENT, datapath=env.nexus, plot=False,
                          override_params={"Spath": reference.out, "subname": "autoreduction"},
                          check_for_prior=True)
    return read_outputs(reference.out)


def scale_of(outputs, name=None):
    header = outputs[name or COMBINED][0]
    return records(header)["scale"]


@pytest.mark.parametrize("runs", SCENARIOS.values(), ids=SCENARIOS.keys())
def test_scale_record_is_the_same_in_any_order(env, batch, runs):
    """Every file of every scenario carries the batch's per-position scale factors (rel 1e-9: V7, a product of
    per-call factors agrees across orders to rounding). Before: the list of the last call (V1)."""
    reduce_runs(env, runs)
    expected = scale_of(batch)
    assert all(v is not None for v in expected)
    for name, (header, _) in read_outputs(env.out).items():
        assert records(header)["scale"] == pytest.approx(expected, rel=1e-9), name


@pytest.mark.parametrize("reverse_q", [False, True], ids=["sequence order is Q order", "reverse Q"])
@pytest.mark.parametrize("runs", [[R1, R2, R3], [R3, R1, R2, R1]], ids=["in order", "out of order"])
def test_scale_record_is_what_was_applied(env, runs, reverse_q):
    """Each position file's R column divided by the stub's unscaled R is the recorded factor: the record states
    what was applied to that position's data, whichever position the merge scaled (Q order is not sequence
    order with reverse_q)."""
    env.reverse_q = reverse_q
    reduce_runs(env, runs)
    outputs = read_outputs(env.out)
    for run in RUNS:
        seq = RUNS[run][0]
        header, data = outputs[partial_name(seq, run)]
        _, unscaled, _, _ = synthetic_curve(seq - 1, reverse_q)
        applied = data[1] / unscaled
        assert applied == pytest.approx(np.full_like(applied, records(header)["scale"][seq - 1]), rel=1e-9), run


def with_settings(env, **changes):
    settings = json.loads(env.settings.read_text())
    settings.update(changes)
    env.settings.write_text(json.dumps(settings))


@pytest.mark.parametrize("runs", [*SCENARIOS.values(), "batch"], ids=[*SCENARIOS.keys(), "batch"])
def test_config_scale_factor_is_the_authored_list(env, runs):
    """R4 (decision F17-1): Config.ScaleFactor stays the authored input in every file. The applied factor lives in
    the Scaling factors line only; written into the config it would be multiplied into the data again when the
    header is used as a settings file (V5)."""
    with_settings(env, ScaleFactor=list(AUTHORED))
    if runs == "batch":
        nrff.reduce_from_file([R1, R2, R3], env.settings, EXPERIMENT, datapath=env.nexus, plot=False,
                              override_params={"Spath": env.out, "subname": "autoreduction"}, check_for_prior=True)
    else:
        reduce_runs(env, runs)
    for name, (header, _) in read_outputs(env.out).items():
        assert records(header)["config_scale"] == AUTHORED, name


@pytest.mark.parametrize("runs", SCENARIOS.values(), ids=SCENARIOS.keys())
def test_lambda_range_is_recorded_per_position(env, runs):
    """The wavelength range each position's reduction used, in the Lambda Range line and in Config's
    LambdaMinUse/LambdaMaxUse, by sequence position. Before: the last call's scalars."""
    reduce_runs(env, runs)
    for name, (header, _) in read_outputs(env.out).items():
        r = records(header)
        assert (r["lambda_min"], r["lambda_max"]) == ([2.7, 2.6, 2.5], [9.5, 9.5, 9.5]), name
        assert (r["config_lambda_min"], r["config_lambda_max"]) == ([2.7, 2.6, 2.5], [9.5, 9.5, 9.5]), name


def test_a_first_call_for_a_later_position_records_null_at_the_gaps(env):
    """Battery row 10 survived without this: the merge rebuilds gaps from its positions, so reduce()'s own gap
    entries reach a file only when no merge runs, as in the first call of a sequence for a later position."""
    reduce_runs(env, [R3])
    for name, (header, _) in read_outputs(env.out).items():
        r = records(header)
        assert r["scale"][:2] == [None, None] and r["lambda_min"] == [None, None, 2.5], name


def test_gap_records_null(env):
    """A position with no run has no applied factor and no range: null in the three records."""
    reduce_runs(env, [R1, R3])
    for name, (header, _) in read_outputs(env.out).items():
        r = records(header)
        assert r["scale"][1] is None and r["lambda_min"][1] is None and r["lambda_max"][1] is None, name
        assert None not in (r["scale"][0], r["scale"][2]), name


def test_autoscale_off_records_the_authored_factor(env):
    """With autoscale off, the factor applied to each position is the authored one."""
    with_settings(env, AutoScale=False, ScaleFactor=list(AUTHORED))
    reduce_runs(env, [R1, R2, R3, R2])
    for name, (header, _) in read_outputs(env.out).items():
        assert records(header)["scale"] == AUTHORED, name


def write_format2_file(out, run, nr_runs, scale, lam=(2.5, 9.5), list_lines=False):
    """A partial file as M1 (header format 2) wrote it: positional logs, but the Scaling factors line holds the
    writing call's list and the Lambda Range line its scalars. With `list_lines`, the Lambda Range line parses as
    per-position lists too: a format-2 file whose lines look like records, which only the marker can refuse."""
    seq = RUNS[run][0]
    q, r, dr, dq = synthetic_curve(seq - 1, reverse_q=False)
    entries = [e or None for e in nr_runs]
    angles = {"THS": [ths(e) if e else None for e in entries], "THI": [THI if e else None for e in entries],
              "ThCen": [thcen(e) if e else None for e in entries]}
    head = "\n".join([
        f"NR_runs = {nr_runs}",
        f"Run Title: {json.dumps({'title': [title(e) if e else None for e in entries]})}",
        f"Scaling factors = {json.dumps({'scale_factor': scale})}",
        (f"Lambda Range = {json.dumps({'lambda_min': [lam[0]] * len(nr_runs), 'lambda_max': [lam[1]] * len(nr_runs)})}"
         if list_lines else f"Lambda Range = {lam[0]}Å to {lam[1]}Å"),
        f"Angles: {json.dumps(angles)}",
        "Header format: 2 (Run Title, Angles and NR_runs are indexed by sequence position)",
        "---" * 20,
        f"Config: {json.dumps({'RBnum': nr_runs, 'ScaleFactor': scale})}",
        "---" * 20,
        "columns = Q, R, dR, dQ (sigma)",
        "---" * 20,
    ])
    np.savetxt(out / partial_name(seq, run), np.column_stack((q, r, dr, dq)), header=head, delimiter="\t")


@pytest.mark.parametrize("legacy", [False, "list lines", True], ids=["format 2", "format 2, list-shaped lines", "no marker"])
def test_prior_without_format_3_does_not_vouch_then_heals_when_rereduced(env, capsys, legacy):
    """R2: before this slug, every writer put the writing call's list into every file, so no format-2 (or legacy)
    entry is provably the file's own. Such a prior gives null for its position's scale and range, with a notice,
    and the position heals when its run is reduced again."""
    for run in (R1, R2, R3):
        if legacy is True:
            write_legacy_file(env.out, run, [R1, R2, R3], [R1, R2, R3])
        else:
            # "list lines" (battery row 4 survived without it): the records parse, so only the marker refuses them
            write_format2_file(env.out, run, [R1, R2, R3], [1, 0.5, 0.25], list_lines=legacy == "list lines")
    reduce_runs(env, [R2])
    out = capsys.readouterr().out
    r = records(read_outputs(env.out)[COMBINED][0])
    assert r["scale"][0] is None and r["scale"][2] is None and r["scale"][1] is not None
    assert r["lambda_min"] == [None, 2.6, None] and r["lambda_max"] == [None, 9.5, None]
    for run in (R1, R3):
        assert f"{partial_name(RUNS[run][0], run)}: the header does not vouch for its scale factor" in out, run
    reduce_runs(env, [R1, R3])
    healed = records(read_outputs(env.out)[COMBINED][0])
    assert healed["lambda_min"] == [2.7, 2.6, 2.5] and None not in healed["scale"]


def rewrite_header_line(path, prefix, line):
    text = path.read_text().splitlines()
    path.write_text("\n".join(line if t.startswith(prefix) else t for t in text) + "\n")


@pytest.mark.parametrize("damage", ["list shorter than NR_runs", "another run at the position"])
def test_inconsistent_format_3_header_does_not_vouch(env, damage):
    """A format-3 prior whose records do not line up with its NR_runs is not vouched for: null, never a guess."""
    reduce_runs(env, [R1, R2, R3])
    path = env.out / partial_name(1, R1)
    if damage == "list shorter than NR_runs":
        rewrite_header_line(path, "# Scaling factors = ", '# Scaling factors = {"scale_factor": [1.0, 0.5]}')
    else:
        rewrite_header_line(path, "# NR_runs = ", f"# NR_runs = [{R2}, {R2}, {R3}]")
    reduce_runs(env, [R3])
    assert records(read_outputs(env.out)[COMBINED][0])["scale"][0] is None


def test_returned_config_keeps_scalar_lambda_use(env):
    """R5: in memory, LambdaMinUse/LambdaMaxUse stay the per-call scalars web_report formats with %6.4g (V4); the
    per-position lists exist in the records and the written header only."""
    reduce_runs(env, [R1, R2])
    config = nrff.reduce_from_file([R3], env.settings, EXPERIMENT, datapath=env.nexus, plot=False,
                                   override_params={"Spath": env.out, "subname": "autoreduction"},
                                   check_for_prior=True)[3]
    assert isinstance(config.LambdaMinUse, float) and isinstance(config.LambdaMaxUse, float)


def test_template_style_caller_gets_no_format_3_marker(tmp_path):
    """V10/R3: a save_results caller that supplies no records (the template path) gets today's lines and marker 2,
    never a format-3 claim it cannot back."""
    from lr_reduction import save_reduced_data
    from lr_reduction.nr_reduction_config import NRReductionConfig

    config = NRReductionConfig()
    config.RBnum, config.ScaleFactor, config.Spath, config.Sname = [R1], [1.5], tmp_path, "template"
    config.LambdaMinUse, config.LambdaMaxUse = 2.5, 9.5
    logs = {"title": [title(R1)], "ths": [ths(R1)], "thi": [THI], "ThCen": [thcen(R1)]}
    q = np.linspace(0.01, 0.1, 5)
    save_reduced_data.save_results({"Q": q, "R": q, "dR": q, "dQ": q}, config, logs)
    lines = [line for line in (tmp_path / "template.dat").read_text().splitlines() if line.startswith("#")]
    assert any(line.startswith("# Header format: 2") for line in lines)
    assert "# Scaling factors = {\"scale_factor\": [1.5]}" in lines
    assert "# Lambda Range = 2.5Å to 9.5Å" in lines
    # a partial record set is no record set (frame F20): format 2 and today's lines, never a half format 3
    save_reduced_data.save_results({"Q": q, "R": q, "dR": q, "dQ": q}, config, {**logs, "scale": [2.0]},
                                   sname="partial")
    partial = (tmp_path / "partial.dat").read_text().splitlines()
    assert any(line.startswith("# Header format: 2") for line in partial)
    assert "# Scaling factors = {\"scale_factor\": [1.5]}" in partial


@pytest.mark.parametrize("runs", [[R1, R2, R3], [R3, R1, R2, R1], "one call"], ids=["in order", "out of order", "one call"])
def test_eight_column_files_carry_the_same_records(env, runs):
    """The _8col files carry the 4-column files' records. "one call" reduces the whole sequence at once, with no
    priors, so reduce()'s own 8-column files are the ones kept (frame rows F2 and F4 survived without it)."""
    for call in ([[R1, R2, R3]] if runs == "one call" else [[run] for run in runs]):
        nrff.reduce_from_file(call, env.settings, EXPERIMENT, datapath=env.nexus, plot=False,
                              override_params={"Spath": env.out, "subname": "autoreduction", "save8col": True},
                              check_for_prior=True)
    outputs = read_outputs(env.out)
    for name in [n for n in outputs if n.endswith("_8col.dat")]:
        assert records(outputs[name][0]) == records(outputs[name.replace("_8col.dat", ".dat")][0]), name


def test_remeasured_position_records_follow_the_winner(remeasured):
    """D-5: position 2's records are the winning run's: its scale is what was applied to R2B's data."""
    reduce_runs(remeasured, [R1, R2, R3, R2B, R2])
    outputs = read_outputs(remeasured.out)
    header, data = outputs[partial_name(2, R2B)]
    _, unscaled, _, _ = synthetic_curve(1, reverse_q=False)
    applied = data[1] / unscaled
    assert records(outputs[COMBINED][0])["scale"][1] == pytest.approx(float(applied[0]), rel=1e-9)
    assert applied == pytest.approx(np.full_like(applied, applied[0]), rel=1e-9)


def test_nr_runs_header_round_trips_numpy_integers(tmp_path):
    """R6: a run passed as a numpy integer, or as digits, is written as a plain int, so read_prior_header can read
    the line back; a run number that is not an integer is refused, never written."""
    from lr_reduction import save_reduced_data
    from lr_reduction.nr_reduction_config import NRReductionConfig

    def save(runs, name):
        config = NRReductionConfig()
        config.RBnum, config.Spath, config.Sname = runs, tmp_path, name
        logs = {"title": [None] * len(runs), "ths": [None] * len(runs), "thi": [None] * len(runs),
                "ThCen": [None] * len(runs)}
        q = np.linspace(0.01, 0.1, 5)
        save_reduced_data.save_results({"Q": q, "R": q, "dR": q, "dQ": q}, config, logs)
        return tmp_path / f"{name}.dat"

    path = save([np.int64(R1), None, "221474"], "numpy")
    assert nrff.read_prior_header(path)["NR_runs"] == [R1, None, R3]
    path = save([float(R1), np.int32(R2)], "integral")  # an integral float is a whole run number (frame F18)
    assert nrff.read_prior_header(path)["NR_runs"] == [R1, R2]
    for bad in (3.5, "x", True):  # True is an int subclass, not a run number (frame F17)
        with pytest.raises(ValueError, match="run number"):
            save([bad], f"bad-{bad}")
