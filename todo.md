# Integrator review gate — `harness-hardening` v2 REJECTED

Supersedes the v1 rejection recorded in this file at `7497b28`
(`git log -- todo.md`).

**Disposition: blocking, per the plan's own declaration.** The Analyst
upgraded `test-reviewer` to blocking for this attempt — *"the slug's
deliverable is self-tested guarantees, so the domain that measures whether
the guarantees hold must gate it"* — and that domain returned three
blocking findings from a fully measured disarm matrix. No Integrator
departure was needed.

**Gate: green.** `test-launcher` **8 passed** (10.56 s), `test-reduction`
**107 passed**, exit 0, 530.14 s, no race. `pixi.lock` restored to
`version: 6` after the run re-stamped it.

**Retry budget: the plan records attempt 2; charter §1 sets N = 3.** A v3
is the last retry before escalation — and `port-overplot-axes-refresh` and
`port-settings-persistence` are *both* already at attempt 3, so all three
active slugs are now at or approaching the cap.

## What v1 was rejected for IS fixed — state this first

- **All 8 self-tests are now falsifiable.** Every mechanism they claim to
  guard was disarmed and each test has at least one mutation that reds it.
  v1's decisive defect is closed.
- **Order-independence is closed.** All 8 pass individually by node id, in
  suite order, and in fully reversed order — including the settings test
  that was green-when-run-alone at v1.
- **F2 (modal cover) is closed and pinned.** The patch targets `QDialog`,
  where `exec_`/`exec`/`open` are actually defined, covering the 6 real
  production sites (`direct_beam.py:113,185,355,363`;
  `roi_selector.py:230,2038`), and returns `QDialog.Accepted` (1). Both the
  wrong-class and wrong-return-value mutations red.
- **The collection hook is closed and fully pinned** — path filter,
  `--timeout` escape, and marker-arming each red under mutation.
- **Sibling slugs are measured clean**: S1 6/6 and S3 32/32 identical under
  the old and new conftest; the `autouse` contract change and the
  base-class patch are safe as written.

## B1 — BLOCKING: the headline F1 fix is deletable in silence

Deleting `conftest.py:20-23` — `_SCRATCH_ROOT`, `setDefaultFormat`, and
the two-format `setPath` loop, i.e. **the entire fix v1 was rejected for
only half-doing** — yields **8 passed**.

The mechanism is correct; it was verified end-to-end. Fixture-less
`QSettings` under the shipped conftest resolves into the scratch tree,
and with `:20-23` deleted resolves to the real config dir:

```
shipped :  /tmp/launcher-tests-settings-*/Unknown Organization.ini
deleted :  $XDG_CONFIG_HOME/Unknown Organization.conf      <- the real ~/.config
```

`test_qsettings_file_lands_in_tmp` cannot see this because it always runs
**under** `isolated_qapp`, whose own per-test `setPath` satisfies the
assertion regardless of whether the import block exists. Confirmed by the
Integrator: `_SCRATCH_ROOT` appears only at `conftest.py:20,23` and no test
references it.

The `NativeFormat` halves (at `:22-23` and at the fixture's `:40-41`) and
the deliberate "never restore `defaultFormat`" decision at `:53-55` are
each separately unpinned — the last of which **is literally the v1 bug**,
restorable with 8 passed.

Failure scenario: a future slug tidies the redundant-looking import block,
all 8 stay green, and every fixture-less test — including all 15 launcher
classes that construct `QSettings()` in `__init__` — resumes writing the
developer's real `~/.config`. F1 re-opened with zero signal.

**Fix (both verified by the reviewer):** a subprocess test with fake
`HOME`/`XDG_CONFIG_HOME` running a fixture-less, import-time `QSettings`
in both formats and asserting nothing escapes; plus
`test_native_format_also_lands_in_tmp` under the fixture. Measured green on
base, red on the import-block deletion, on an `IniFormat`-only import
redirect, and on an `IniFormat`-only fixture narrowing.

## B2 — BLOCKING: the timeout self-test does not test the shipped configuration

`test_harness.py:134` writes its own ini into `tmp_path`:

```python
(tmp_path / "pytest.ini").write_text("[pytest]\ntimeout = 10\ntimeout_method = thread\n")
```

A `pytest.ini` in `tmp_path` makes `tmp_path` the rootdir, so the repo's
`pyproject.toml` is **never read** by the inner run. Measured: **deleting
`timeout_method = "thread"` from `pyproject.toml:222` leaves 8 passed.**

So the docstring's claim — *"self-test for the shipped configuration: if
pytest-timeout goes missing or `timeout_method` regresses to signal, this
goes red"* — is half false. It is a self-test for the *plugin*, not for the
configuration. This is **v1's B5a, unclosed**: the same test, the same
gap, a different mechanism.

Failure scenario: `timeout_method` is dropped in a merge. `pixi run
test-launcher` is unaffected (it passes the flag on the CLI), but the
documented inner loop `pytest launcher/tests/test_x.py -k …` arms markers
using `signal` — which cannot interrupt a Qt modal. The exact hang class
this slug exists to kill returns silently, with a green suite.

**Fix (verified):** replace `:134` with
`shutil.copy(<repo>/pyproject.toml, tmp_path / "pyproject.toml")` and drop
the now-unavailable `b"timeout method: thread"` assertion at `:154`; the
generated file's `@pytest.mark.timeout(10)` supplies the budget. Measured
green on base (10.4 s), red when `pyproject.toml:222` is deleted.

## B3 — BLOCKING: the self-test deterministically creates an unkillable orphan

Reproduced 2/2 on the **shipped** version: start the outer pytest,
`SIGKILL` it at t=3 s, and the inner `test_hangs.py` survives
indefinitely — two were observed alive at 14 minutes. The inner run's own
`timeout = 10` and `@pytest.mark.timeout(10)` do **not** save it.

Mechanism, isolated by bisection:

| inner output target | outer SIGKILLed | inner process |
|---|---|---|
| `capture_output=True` (pipes) — **as shipped, `:147`** | yes | **alive at +25 s, orphan** |
| redirected to a file | yes | self-exits at ~+9 s |

When the parent dies its pipe read-ends close, and `pytest_timeout`'s
`finally: terminal.flush()` raises `BrokenPipeError` **before** its
`os._exit(1)` — which is therefore never reached, leaving the inner
process blocked in Qt's C++ event loop forever. The same path applies when
the parent dies by OOM-kill or by its own thread-method timeout (`os._exit`
in the parent skips `subprocess.run`'s child cleanup).

**There is live proof on this machine.** PID 840666 has been spinning at
99.9% CPU for ~47 minutes running
`.../pytest-81/test_timeout_backstop_fires0/test_hangs.py`. Its tree shows
`pytest.ini` = `timeout_method = thread` with **no `timeout =`**, and
`test_hangs.py` = `while True: pass` — an *earlier* iteration that armed no
timeout at all. The shipped version is better but not fixed.

**A test whose failure mode is an unbounded orphan is not acceptable in the
slug chartered to eliminate unbounded orphans.**

**Fix (verified, 2 lines):** write inner output to files in `tmp_path`
rather than pipes, and wrap the call in `try/finally: proc.kill()`.
Optionally `start_new_session=True` + `os.killpg`. Measured: file mode
self-exits ~9 s after a parent SIGKILL; pipe mode never does.

## Advisory — carry into v3

- **The entire teardown rewrite is untested** (`conftest.py:48-66`):
  reverting to v1's shape, moving the identity restores after the drain,
  and deleting `sendPostedEvents(DeferredDelete)` each yield 8 passed.
- **And the drain's stated rationale is measurably wrong on this env.** It
  is justified by "one QApplication is shared by the whole session", but a
  **new QApplication is created for every test** (3 distinct objects across
  3 tests), because the fixture drops the only Python reference at
  teardown — so Qt's own destructor already frees everything, confirmed
  under the old conftest and with the drain removed entirely. The drain is
  currently unobservable, which is *why* the mutations above cannot be
  caught. This is the same "stated justification contradicted by
  measurement" shape v1 was rejected for, and it is the third cycle running
  in which it appears. **Either add a consequence test or correct the
  rationale — do not leave an unmeasured justification in the record.**
- `_ROOT_POISON` (`test_harness.py:15`) is inert (deleting it changes
  nothing) and actively misleading: after the first fixture teardown it is
  destroyed, and its comment contradicts the slug's own premise that
  `setPath` is consulted per construction. Delete it or make it
  load-bearing.
- The fixture's `setPath` is never undone, so after the first fixture-using
  test the process root is stuck at *that test's* `tmp_path`. Isolation
  from `~/.config` still holds, so this is not a leak — but the docstring's
  "keeps every test inside a scratch tree" is true only until first
  fixture use.
- `_SCRATCH_ROOT` is never cleaned: **86** `launcher-tests-settings-*`
  directories are currently on this host, one per process importing the
  conftest. Use an `atexit` `rmtree`.
- Cost: the timeout self-test is 10 s of a 10.5 s suite, and its *red*
  costs 100 s. Drop the inner budget to ~3 s and the outer to ~30 s.
- `no_qfiledialog` is autouse while `no_qmessagebox` is opt-in, so a
  sibling test that forgets it still blocks at `roi_selector.py:230`/`:2038`
  — T1's own sites — protected only by the 120 s marker. Worth stating in
  the PR body, since the commit describes the patch as covering "every
  blocking site in `launcher/`".
- From the advisory ui-aspects review, all real-mechanism/no-current-trigger:
  the drain calls `close()` outside modal cover (a `closeEvent` that execs
  would hang teardown — reproduced; no launcher widget has one today); the
  drain destroys module/session-scoped widgets a longer-scoped fixture owns
  (reproduced; nothing trips it today); `QDialog.Accepted` is the wrong
  return for `QMessageBox` *instances*, whose `exec_` contractually returns
  a `StandardButton` (forward-looking); `open()` cover is partial because
  `QMessageBox`/`QFileDialog` define their own; and the patched `exec_`
  runs no event loop, so `dlg.exec_(); dlg.result()` inverts.

## Suggested v3 order

1. **B3** first — it is actively degrading the machine and is 2 lines.
2. **B1** — the two verified tests; this is the finding that keeps F1 from
   silently regressing.
3. **B2** — copy the real `pyproject.toml` instead of hand-writing an ini.
4. The teardown rationale (advisory 1–2): add a consequence test **or**
   correct the claim to match the measured per-test-QApplication behavior.
5. `_ROOT_POISON`, `_SCRATCH_ROOT` cleanup, and the cost reduction.

**Note for the Analyst on scope:** v2 closed everything v1 named. The v3
findings are all of the form "the fix is right but nothing pins it" —
which is the natural end state of a hardening slug and is worth one more
pass, but if v3 does not converge, the escalation path is the honest
outcome rather than a v4 under a different name.
