#!/usr/bin/env python3
"""Mutation battery for roi-estimate (T1 slug 1). Run from the repo root:

    pixi run python plans/scripts/roi_estimate_mutations.py

Restore safety per the Developer contract: originals held in memory AND a
mode-600 backup, restore in a ``finally``, sha256 compared after every mutation
with an abort on mismatch, and a per-invocation timeout far below the 600 s
harness ceiling.
"""

import hashlib
import os
import signal
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MOD = os.path.join(REPO, "src/lr_reduction/roi_estimate.py")
T = "tests/unit/lr_reduction/test_roi_estimate.py"

MUTATIONS = [
    (1, "motors read [0] (the angle BEFORE the move) instead of [-1]",
     'meta[name] = float(f[dpath][-1])', 'meta[name] = float(f[dpath][0])'),
    (2, "a missing chopper log defaults instead of refusing",
     '''                raise KeyError(
                    f"{path} has no chopper log at {dpath!r} — the wavelength "
                    f"band cannot be derived for this run, and must not be guessed"
                )''',
     '                return [2.4, 5.8]'),
    (3, "inline the band maths (a fifth copy) instead of calling the library",
     '        return nr_tools.get_lam_range(chopper_lam, chopper_speed)',
     '        return [chopper_lam - 1.75 * 60.0 / chopper_speed - 0.15,\n                chopper_lam + 1.75 * 60.0 / chopper_speed - 0.15]'),
    (4, "fork the id unpacking instead of calling get_y_tof",
     '''    _, y_tof, _ = binary_processing.get_y_tof(
        tof_array, event_id, e_offset, list(lowres), pcharge, n_y=n_y, n_x=n_x
    )''',
     '''    y_tof = np.zeros((n_y, len(tof_array)))
    np.add.at(y_tof, (event_id % n_y, np.zeros(len(event_id), dtype=int)), 1)'''),
    (5, "hard-code the detector shape instead of the instrument DB",
     '''    settings = nr_tools.read_settings(start_time)
    return int(settings["num_x_pixels"]), int(settings["num_y_pixels"])''',
     '    return 256, 304'),
    (6, "drop the no-counts guard",
     '''    if counts.size == 0 or not np.any(counts > 0):
        raise ValueError("no counts on the detector — no peak can be estimated")''',
     '    pass'),
    (7, "drop the contrast guard",
     '''    if contrast < min_contrast:
        raise ValueError(
            f"contrast {contrast:.2f} is below {min_contrast} — the detector is "
            f"featureless here, so the largest bin is noise, not a peak"
        )''',
     '    pass'),
    (8, "drop the low-side FIT CHECK (the real protection; the clamps were dead)",
     '    if high_edge - width + 1 >= 0:',
     '    if True:'),
    (9, "drop the no-room refusal",
     '''    raise ValueError(
        f"no room for a {width}-pixel background with a {gap}-pixel gap beside "
        f"peak {peak_range} on a {n_y}-pixel detector"
    )''',
     '    return 0, width - 1'),
    (10, "normalise by the DASlogs time series instead of the run total",
     '        pcharge = np.asarray(f["entry/proton_charge"][:])',
     '        pcharge = np.asarray(f["entry/DASlogs/proton_charge/value"][:])'),
]


def sha(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def verify_baseline_matches_head():
    """Refuse to start unless the target matches its committed blob.

    The first version took `sha(MOD)` of whatever was on disk and called that
    "clean". With zero git references, a leftover from a killed run BECAME the
    baseline: the battery then restored to the mutated text and printed
    `restored: OK` with the mutation still in the file. That is precisely the
    failure `todo-mutation-harness-restore-safety` was written about, and since
    this is the campaign's first committed battery it is the reference
    implementation — so the bug propagates by being copied.

    Detection complete: compare against `git show HEAD:<path>`, not against the
    working tree. A battery that cannot tell dirty from clean cannot honestly
    report anything.
    """
    rel = os.path.relpath(MOD, REPO)
    try:
        blob = subprocess.run(
            ["git", "show", f"HEAD:{rel}"], cwd=REPO,
            capture_output=True, check=True,
        ).stdout
    except (subprocess.CalledProcessError, OSError) as exc:
        raise SystemExit(f"ABORT: cannot read the HEAD blob for {rel}: {exc}") from exc
    if hashlib.sha256(blob).hexdigest() != sha(MOD):
        raise SystemExit(
            f"ABORT: {rel} differs from HEAD — the working tree is dirty, so it "
            f"cannot serve as the mutation baseline. Commit or restore it first. "
            f"(A leftover from a killed run looks exactly like this.)"
        )


def main():
    # Baseline first: never adopt the working tree sight-unseen.
    verify_baseline_matches_head()

    orig = open(MOD, encoding="utf-8").read()
    clean = sha(MOD)
    fd, bak = tempfile.mkstemp(prefix="roi-", suffix=".bak")
    os.close(fd)
    os.chmod(bak, 0o600)
    open(bak, "w", encoding="utf-8").write(orig)

    # A `finally` does not run when the process is SIGTERM'd, which is exactly
    # how the 600 s harness ceiling kills a long battery — demonstrated: killed
    # at the limit, the target was left mutated. Restore from the handler too.
    def _restore_and_die(signum, _frame):
        with open(MOD, "w", encoding="utf-8") as fh:
            fh.write(orig)
        print(f"\nsignal {signum}: restored {MOD} before exiting", file=sys.stderr)
        raise SystemExit(128 + signum)

    for _sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(_sig, _restore_and_die)

    rows = []
    try:
        for row, desc, old, new in MUTATIONS:
            if orig.count(old) != 1:
                rows.append((row, desc, f"ANCHOR x{orig.count(old)}"))
                print(f"[{row}] ANCHOR MISS — {desc}")
                continue
            try:
                open(MOD, "w", encoding="utf-8").write(orig.replace(old, new))
                proc = subprocess.run(
                    [sys.executable, "-m", "pytest", "-q", "--no-header", "-p",
                     "no:cacheprovider", "--timeout=90", "--timeout-method=thread", T],
                    cwd=REPO, capture_output=True, text=True, timeout=240,
                )
                lines = [x for x in proc.stdout.strip().splitlines() if x.strip()]
                observed = lines[-1] if lines else f"exit {proc.returncode}"
            except subprocess.TimeoutExpired:
                observed = "HUNG (harness timeout)"
            finally:
                open(MOD, "w", encoding="utf-8").write(orig)
                if sha(MOD) != clean:
                    raise SystemExit(f"ABORT: {MOD} did not restore cleanly")
            rows.append((row, desc, observed))
            print(f"[{row}] {desc}\n     {observed}")
    finally:
        open(MOD, "w", encoding="utf-8").write(orig)
        print("restored:", "OK" if sha(MOD) == clean else "*** DIRTY ***")
        os.unlink(bak)
    print("\n=== ledger rows ===")
    for row, desc, observed in rows:
        print(f"| {row} | {desc} | {observed} |")


if __name__ == "__main__":
    main()

# Measured 2026-09-20 on feature/roi-estimate (10 rows, 10 red):
#
# | # | mutation | observed |
# |---|---|---|
# | 1 | motors read [0] not [-1]                              | 2 failed |
# | 2 | a missing chopper log defaults instead of refusing    | 2 failed |
# | 3 | inline the band maths (a fifth copy)                  | 1 failed |
# | 4 | fork the id unpacking instead of get_y_tof            | 2 failed |
# | 5 | hard-code the detector shape                          | 1 failed |
# | 6 | drop the no-counts guard                              | 1 failed |
# | 7 | drop the contrast guard                               | 1 failed |
# | 8 | drop the low-side fit check                           | 2 failed |
# | 9 | drop the no-room refusal                              | 1 failed |
# | 10| normalise by the DASlogs series not the run total     | 5 failed |
#
# Rows 5 and 8 SURVIVED the first pass and are the useful part of the record.
#
# Row 5: `settings.json` holds exactly one entry for each pixel count, so the
# database value and the literal 256/304 agree and no value-equality assertion
# can separate them. The guard now asserts PROVENANCE — move the database under
# the function and the answer must move — which is the property the slug needs,
# since the failure being prevented is a future geometry change a literal would
# not follow.
#
# Row 8: the `max(0, ...)` / `min(n_y - 1, ...)` clamps it mutated were DEAD
# CODE. Each sat inside a branch whose own condition already forbade the
# out-of-range case, so deleting them could not change any value. They were
# removed rather than kept as reassurance, the guard became a sweep over every
# peak position, and the mutation was retargeted at the fit check that is the
# actual protection.
