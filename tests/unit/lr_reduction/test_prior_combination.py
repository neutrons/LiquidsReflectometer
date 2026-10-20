"""Re-reducing runs of a sequence must not change the new-workflow output headers.

Autoreduction calls ``reduce_from_file([run], ..., check_for_prior=True)`` once per
run; each call merges the new run with the files already in the output folder
(the "priors") and rewrites every file of the sequence.  The ``# Run Title:``,
``# Angles:`` and ``NR_runs`` headers are lists indexed by sequence position
(``seq_num - 1``), padded with ``null`` up to the highest position present, and
each position's entry comes from the same source as that position's data -- so
any number of reductions, in any order, gives the same per-run log headers
(``NR_runs``, ``Run Title``, ``Angles``) and data as reducing the sequence once.
(``Scaling factors`` and ``Lambda Range`` still describe only the last call:
tasking finding F17, not covered here.)

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

    def fake_reduce_single_run(self, i, rb_num):
        seq = RUNS[rb_num][0]
        assert i == seq - 1, "reduce() must pass the run's sequence position"
        q, r, dr, dq = synthetic_curve(i, state.reverse_q)
        zeros = np.zeros_like(q)
        # reduce() takes the title from self.log_values and the angles from the returned log_vals
        self.log_values = {"title": title(rb_num), "ths": ths(rb_num), "thi": THI, "ThCen": ths(rb_num) + THCEN_OFFSET}
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
        header = {}
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
    return {name: {k: v for k, v in header.items() if k != "format"} for name, (header, _) in outputs.items()}


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
    marker = [i for i, line in enumerate(lines) if line.startswith("# Header format: 2")]
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

    def remeasured_stub(self, i, rb_num):
        result, config, logs = stub(self, i, rb_num)
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
