# Integrator review gate — `harness-hardening` v1 REJECTED

**Tests are green** (`test-launcher` 5 passed; `test-reduction` 107 passed,
539.02 s, exit 0; pre-commit clean; diff exactly the three planned files;
`pixi.lock` untouched at `version: 6`). This rejection is entirely from the
review gate.

**Disposition: blocking — a flagged departure from the plan's advisory
declaration.** The plan declares both domains advisory, and its stated
rationale for the test-reviewer is procedural: *"advisory — it authored the
findings; this closes them."* That premise is what failed. Measured
verification says **F3 is closed for the stated inner loop, F1 is only
partially closed, and F2 is not closed for a single production call
site** — and all three new self-tests are defective. A slug whose entire
deliverable is *"the harness's guarantees are real and self-tested"* has
not delivered when the guarantees are documented but absent.

Both reviewers reached this independently, by probe against this exact
binding (qtpy 2.4.3 → PyQt5 5.15.11, Qt 5.15.15), and the Integrator
verified the load-bearing claims separately. **Analyst override remains
cheap and I will not re-litigate it:** re-tag `qa/harness-hardening` at an
unchanged SHA and I will open the draft PR with all of this recorded as
advisory instead.

## B1 — F1: a write to the real config was reproduced against the NEW fixture

The mechanism is right and deserves credit: `QSettings.setPath()` is
genuinely **not** subject to the cached-root defect — measured, a second
`setPath` is honored per construction where a second `XDG_CONFIG_HOME`
redirect is ignored. That was the correct choice.

But the redirect covers **`IniFormat`/`UserScope` only**, and
`conftest.py:46` restores `defaultFormat` to `NativeFormat` at teardown.
On Linux `NativeFormat` keeps its own path entry that
`setPath(IniFormat, …)` never touches. End-to-end reproduction against the
new fixture (real config stood in for by a redirected fake HOME):

```
t0 (a test WITHOUT isolated_qapp) -> <realconfig>/ORNL/lr_reduction_new_launcher.conf
t1 (inside isolated_qapp)         -> <tmp_path>/test-org-t1/test-app-t1.ini   OK
t2 (after teardown)               -> <realconfig>/ORNL/lr_reduction_new_launcher.conf
REAL-config file created? True    contents: [probe] t0=1  t2=1
```

So the protection is **window-scoped, not process-wide**, and the plan's
failure-mode matrix claim — *"setPath applies process-wide once set by the
first fixture use"* — is false as implemented. This is the same shape as
the S3 defect: a stated justification contradicted by the code's own
ordering. It is not hypothetical here: **all 15 launcher app classes
construct `QtCore.QSettings()` in `__init__`**, so any sibling-slug test
that instantiates a widget without `isolated_qapp` writes to the
developer's real config. Sibling branches already contain fixture-less
tests; they merely happen not to construct `QSettings` yet.

**Fix:** do `setDefaultFormat(IniFormat)` + `setPath` for **both**
`IniFormat` and `NativeFormat` at conftest **import** scope (before any
test or collection-time construction), and do **not** restore
`defaultFormat`. Restore org/app/domain only. Note `setPath` has no Qt
undo, so the permanent half is the harmless half — the current code
restores exactly the wrong one.

## B2 — F2: not closed for any production call site, and the patched surface is dead

Measured under the new fixture, offscreen, with an 8 s hard kill:

| call form | result | production sites |
|---|---|---|
| `QMessageBox().exec_()` | returns 1024 | **0** |
| `QDialog().exec_()` | **BLOCKS** | `roi_selector.py:230`, `:2038`, `direct_beam.py:113`, `:185` |
| `QMessageBox` subclass with a Python `exec_` override | **BLOCKS** | `direct_beam.py:355`, `:363` |
| `QFileDialog.getExistingDirectory(...)` (no opt-in fixture) | **BLOCKS** | 26 sites |

Root cause: `'exec_' in QMessageBox.__dict__` is `False` — it is inherited
from `QDialog`. Patching `QMessageBox.exec_` installs a shadowing attribute
on `QMessageBox` alone; `QDialog.exec_` is untouched, and
`CdSettingsDialog`/`ModeratorDialog` define their own `exec_` in the
subclass `__dict__`, so the MRO never reaches the patch.

**The patch covers a surface with zero production sites.** Independently
confirmed by the Integrator: `launcher/` constructs 8 `QMessageBox`
instances (`dynamic_30Hz`, `refracted`, `xrr`, `off_spec`, `quick_reduce`,
`dynamic_60Hz`, `sld_calculator`, `reduction`) and **never calls `exec_()`
or `show()` on any of them** — they are configured and dropped, which is a
separate pre-existing bug (the error dialog never appears). And the
direction of travel makes it worse: S1 converts the last `QMessageBox`
dialogs to `QDialog`, after which this cover protects nothing at all.
`test_instance_exec_neutralized` is a live test of a dead surface.

`roi_selector.py` is the campaign's named subject; **T1 is the slug that
hits `:230`/`:2038` and gets a hang, not a failure** — precisely the
10-hour-orphan class this harness exists to kill.

**Fix:** patch `QtWidgets.QDialog.exec_`/`exec`/`open` (which `QMessageBox`
and `QFileDialog` inherit, covering all forms in one stroke). Keep the
`hasattr` guard over both names — in PyQt5 they are distinct slots.
`monkeypatch` restores this cleanly: because the attribute is inherited it
records `notset` and undoes with `delattr`, so no copy is pinned onto the
subclass.

**Return-value hazard that must be fixed with it:** the two halves inject
opposite semantics. Statics return `None`; instance `exec_` returns
`QMessageBox.Ok` (1024). Every real comparison is `== 1`
(`direct_beam.py:355,363`), `== QDialog.Accepted` (`roi_selector.py:230`)
or `!= QDialog.Accepted` (`:2038`). **1024 ≠ 1.** Latent today only because
the patch misses those sites — it goes live the moment the obvious fix
lands, at which point every accept path silently becomes a cancel path and
the tests stay green. The replacement must return `QDialog.Accepted`, or
better, be parametrizable so a test states which branch it wants.

Also: `no_qfiledialog` has the four getter names **complete and correctly
shaped** (verified: `(str, str)`, `(list, str)`, bare `str`; `**_kwargs`
absorbs the `options=` keyword used at several sites) — good work — but it
is opt-in and **nothing opts in**, so the 26 sites remain unprotected by
default. Make it `autouse=True` within `launcher/tests`; an author wanting
a chosen path overrides locally, as its own docstring already says. One
nit: `getOpenFileNames` binds a single list object reused across calls; a
caller that mutates it poisons the next call.

## B3 — the timeout hook is session-global, contradicting an acceptance criterion

`conftest.py:101` states *"This conftest governs only launcher/tests, so
the reduction suite's timings are untouched."* False.
`pytest_collection_modifyitems` dispatches through the **global** hook
relay, so once this conftest is registered it receives the entire session's
item list, and the loop has no path filter. Reproduced with a two-tree
scaffold: both `tests/` items came back carrying `timeout(120)`.

With `timeout_method = "thread"`, a fired timeout calls `os._exit()` — the
run dies mid-suite with no FAILURES section, no summary, no `coverage.xml`.
It does not fire under today's two separate pixi invocations, which is why
the commit's verification could not see it: that check confirms the ini is
inert, it never exercises the hook. It fires on `pytest .`,
`pytest tests launcher/tests`, or the moment `launcher/tests` is added to
`testpaths` — the natural next step as S1–S3 land.

This makes the plan's acceptance criterion *"reduction suite green (its
runtime profile unchanged — **no new per-test timeouts there**)"* true only
by accident of invocation shape, not structurally.

**Fix (one line):** filter on
`item.path.is_relative_to(Path(__file__).parent)` before adding the marker.
The docstring then becomes true. Also make the marker conditional on
`item.config.getoption("timeout", None) is None`, or it silently overrides
`--timeout` on the command line — which breaks the very RED-phase technique
the RED commit documents (*"run it with `--timeout=10` to keep RED cheap"*).

## B4 — the teardown can skip its own restores, and the drain is a measured no-op

`conftest.py:39-46` runs the widget drain **before** the four identity
restores, unguarded. `topLevelWidgets()` includes parented dialogs, so if
closing one entry destroys another in the same snapshot the loop raises
`RuntimeError: wrapped C/C++ object ... has been deleted` and the restores
never run — leaving `test-org-…` and `IniFormat` installed process-globally.
The leak then self-propagates: the next test snapshots the polluted org and
faithfully restores *that*. Reachable via `direct_beam.py:57,133`
(`self.dlg = QDialog(parent)`). The setup half is likewise unguarded — if
`QApplication([])` raises, teardown never runs at all.

And the drain does not drain. `processEvents()` does not dispatch
`DeferredDelete`; measured:

```
after drain (close + deleteLater + processEvents)   topLevelWidgets=2, isdeleted=False
after sendPostedEvents(None, QEvent.DeferredDelete) topLevelWidgets=0, isdeleted=True
```

So every widget any test creates stays alive for the whole session and each
teardown re-iterates the accumulated list (O(N²)). The tail item "no
top-level-widget drain between tests" is recorded as closed but is not.

**Fix:** `try/finally` with the restores in the `finally` and ahead of the
drain; guard each widget with `sip.isdeleted()`/`except RuntimeError`; add
`sendPostedEvents(None, QEvent.DeferredDelete)` after `processEvents()`.

## B5 — all three self-tests are defective

This is the finding that decides the disposition. The slug's deliverable is
*self-tested* guarantees.

- **`test_timeout_backstop_fires` cannot fail for what it claims.** It
  spawns pytest with `--timeout=5 --timeout-method=thread` **on the
  subprocess CLI**, and the subprocess's rootdir is `/tmp`, so the repo's
  `pyproject.toml` is never loaded. **Delete both F3 changes and this test
  still passes.** It detects only "pytest-timeout is missing/broken." The
  commit's framing ("or the method silently reverts to signal, this goes
  red") is not true. *Fix:* copy `conftest.py` into `tmp_path` and pass no
  `--timeout`, asserting the conftest-armed `timeout: 120.0s / method:
  thread` header.
- **`test_qsettings_file_lands_in_tmp` is order-dependent.** Measured both
  ways against a simulated pre-fix fixture: run **after** its neighbour →
  FAIL (bug caught); run **alone or first** (`-k`, the exact inner loop this
  slug exists to protect) → **PASS, bug invisible**. Its redness is a side
  effect of a neighbouring test constructing `QSettings` first. *Fix:*
  construct a `QSettings` at module import so the cached root is always
  poisoned.
- **`test_instance_exec_neutralized` is falsifiable but guards nothing** —
  zero production sites, going to zero-forever after S1 (B2).

## What is genuinely good (verified, not assumed)

- **F3 is really closed for the inner loop.** For `pytest
  launcher/tests/test_x.py -k …` from anywhere in the tree, rootdir
  resolves to the repo, the conftest loads, the marker is added, and
  `pytest_timeout.get_params()` takes `timeout=120` from the marker and the
  method from the new ini key. Both halves load-bearing and correct. And
  `timeout_method` alone genuinely arms nothing, so CI and the pixi tasks
  are safe — the reasoning in the commit is right.
- **`setPath` over the XDG env var was the correct mechanism choice**, and
  the docstring's explanation of the caching defect is accurate.
- **`QFileDialog` getter names are complete** for today's code (16
  `getExistingDirectory` + 10 `getOpenFileName` = the 26 sites F2 named),
  with correct return shapes.
- **`ARG001` per-file-ignore is correctly scoped and load-bearing**
  (verified with ruff directly, with and without the config); it cannot
  reach `src/` or `tests/`.
- **Single-QApplication handling is correct** — created once, reused via
  `instance()`, never destroyed. Destroy/recreate is a known crash source
  and was rightly avoided.
- **`monkeypatch` teardown is clean** for every name patched, both the
  inherited and own-`__dict__` cases. No cross-test leak from the patching
  mechanism.
- **No regression risk to S1/S2/S3**: no autouse fixture, names and
  signatures unchanged, additive contract holds. S3's
  `test_isolation_from_legacy_launcher` now lands in `tmp_path` — the
  specific leak cited in that rejection **is** closed by this slug.

## Suggested v2 order

1. **B2** — patch `QDialog.exec_`/`exec`/`open`, return `QDialog.Accepted`,
   make `no_qfiledialog` autouse. Highest value: it is the only finding that
   already has a named victim (T1 / `roi_selector.py:230`, `:2038`).
2. **B1** — module-scope `setPath` for both formats; drop the
   `defaultFormat` restore.
3. **B3** — one-line path filter; make the marker respect an explicit
   `--timeout`.
4. **B4** — `try/finally`, restores first, `sendPostedEvents`.
5. **B5** — make each self-test able to fail: mirrored ini for the timeout
   test, import-time poisoning for the settings test, and retarget the
   exec test at `QDialog` once B2 lands.

Tail item worth carrying: the `ARG001` ignore now blesses the fixture-argument
form, but 4 of 5 tests in `test_harness.py` still use `usefixtures` marks.
Convert them, and say which form is blessed in the PR body so S2/S3 authors
follow one convention.
