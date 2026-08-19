# Integrator review gate — `pixi-lock-format-guard` v2 REJECTED

Supersedes the v1 rejection recorded in this file at `3c1808b`
(`git log -- todo.md`).

**Gate: green** (`test-launcher` 2 passed, `test-reduction` 107 passed,
exit 0, 526.34 s). **Disposition: blocking, per the plan's own
declaration** — the Analyst upgraded `design-reviewer` to blocking for
this attempt ("when the slug's entire deliverable is a guard, the domain
that judges whether the guard guards must gate it"), and that domain
returned two blocking findings. No Integrator departure was needed this
time. Retry budget: **attempt 2; charter N = 3, so v3 is the last.**

## All four v1 findings are fixed, and the B4 fix is exemplary

Verified independently, not taken on report:

- **B1**: a ≥0.68 shim now yields exit 1 with a v7-specific message naming
  `pixi self-update --version 0.67.2`; the generic drift path says
  "regenerate UNDER A pixi <= 0.67.x". Nothing tells anyone to run
  `pixi lock` any more.
- **B2**: from a subdirectory → exit 0, **no stray file**, repo lock
  untouched; outside any repo → exit 1, "failing closed".
- **B4**: the header now separates MEASURED from UNVERIFIED and adds a
  durable guard — *"Do not restate this as measured until someone measures
  it."* The reviewer diffed the real committed lock against the real
  check-mutated tree and confirmed the documented mutation triple
  (`version`, `sha256`, `editable: true`) matches line for line. **This is
  the right way to answer a stated-vs-measured finding** — not by deleting
  the claim but by scoping it and blocking its recurrence.
- **The self-package classifier is narrowed** (v1 stripped the whole
  stanza): `requires_dist` added → blocked; `requires_python` changed →
  blocked; `name` changed → blocked; benign restamp → passes.
- **The deviation from the plan's literal wording is correct and flagged**:
  stripping only `version`/`sha256`, as the plan said, is non-empty on a
  benign push (`editable: true` alone) and would block every ordinary push.
  Measured, deviated, flagged — correct-and-flag done right.
- **CI placement verified** by parsing the YAML: the format check precedes
  `setup-pixi` in all three pixi jobs, and it checks **format**, not
  consistency. The escape hatch is the right one: `SKIP=pixi-lock-check`
  disables only the network-dependent wrapper while the `always_run`
  format tripwire stays armed.

## B5 — BLOCKING: the trap is armed over an empty snapshot and can truncate the tracked lock

`scripts/pixi_lock_check.sh:44-49`:

```
44  snap=$(mktemp) || exit 1
48  trap 'cp "$snap" pixi.lock 2>/dev/null; rm -f "$snap"' EXIT INT TERM
49  cp pixi.lock "$snap" || exit 1
```

`mktemp` has already created a **0-byte** file. The trap is armed at `:48`
**before** `:49` populates it, so any exit in that two-line window copies
emptiness onto the real repo-root `pixi.lock`. Reproduced twice by the
reviewer, two independent ways:

| Trigger | Result |
|---|---|
| snapshot `cp` fails with a healthy repo FS (ENOSPC/quota on `/tmp`) | `pixi.lock` **327276 → 0 bytes**, `M pixi.lock` |
| `ulimit -f` (shim-free; a full or quota'd TMPDIR) | `pixi.lock` **327276 → 1024 bytes** |

**This is strictly worse than the v1 defect it replaced.** v1's fail-open
wrote a 0-byte lock into a *subdirectory*, leaving the real lock intact.
v2 truncates the **tracked repo-root lock** — and in the partial-truncation
variant `head -1` still reads `version: 6`, so **all three tripwires this
slug ships wave the corpse through**: the script's own format check
(`:56`), the `pixi-lock-format` pre-commit hook
(`.pre-commit-config.yaml:29`), and the new CI step. The guard can damage
the very artifact it exists to protect, in a way its own guards cannot see.

Failure scenario: `/tmp` fills on a shared facility workstation; the dev
pushes; the hook silently truncates `pixi.lock`; `git status` shows a
modification they did not make; a `git commit -am` carries it; pre-commit
passes; CI's format check passes; the job dies deep inside `setup-pixi` on
an inscrutable parse error — the exact outcome the CI step was added to
prevent. The natural recovery instinct is `pixi lock`, which on a ≥0.68
machine is the v7 rewrite.

**Attribution — this is not a Developer slip.** The plan's own v2 fix
sketch prescribes this exact ordering verbatim. **The plan needs the same
correction, or v3 reproduces it.**

## B6 — BLOCKING: the INT/TERM handler does not exit, so the script runs on with a deleted snapshot

A bash `trap … INT` handler returns to the interrupted script unless it
exits. This one restores the lock **and `rm -f "$snap"`**, then execution
continues at `:52`; by `:95`, `cmp -s pixi.lock "$snap"` compares against a
file that no longer exists. Measured, SIGINT during the solve with the lock
already converted to v7:

```
pixi-lock-check: restored pixi.lock (the check mutated it)
awk: fatal: cannot open file `/tmp/tmp.i0ueVZZ4k9' for reading
pixi-lock-check: pixi.lock does not match pyproject.toml — regenerate it UNDER A pixi <= 0.67.x
0a1,8283
> version: 6
… 8,287 lines — the entire lock dumped as a "drift hunk"
```

The lock itself **is** correctly restored, so the v2 transcript is accurate
as far as it measured; what it did not measure is the script's *verdict*
after the interrupt, which is a fabricated drift accusation. Two
consequences in this domain: the comment at `:45` ("One restore path, and
it survives an interrupt") is a documented guarantee only half provided —
the same stated-vs-measured shape B4 was rejected for, now one line below
the B4 fix; and a Ctrl-C answered with a false accusation plus 8k lines of
spew is exactly what drives a maintainer to `--no-verify`, which
`developer.rst` correctly identifies as the thing that disarms the format
tripwire.

## Fix for both (reviewer-verified, no regressions)

```bash
snap=$(mktemp) || exit 1
# Snapshot BEFORE arming the trap: an EXIT trap installed over an empty
# mktemp file would copy that emptiness onto pixi.lock if this cp failed.
cp pixi.lock "$snap" || { rm -f "$snap"; exit 1; }
restore() { cp "$snap" pixi.lock 2>/dev/null; rm -f "$snap"; }
trap restore EXIT
trap 'restore; trap - EXIT; exit 130' INT
trap 'restore; trap - EXIT; exit 143' TERM
```

Measured with this patch: snapshot-`cp` failure → lock intact at 327276
bytes, exit 1; SIGINT → exit 130, **zero output**, lock byte-identical;
benign / v7 / drift stubs unchanged (0 / 1 / 1). `bash -n` and `shellcheck`
clean.

## Advisory — carry into v3

- **`version: 6` is written 5×, `0.67.2` 5×, across 4 files.** No single
  source; the CI snippet is copy-pasted three times. A
  `scripts/check_lock_format.sh` invoked by CI (`run:`) and pre-commit
  (`entry:`) makes retirement one deletion.
- **Retirement is currently archaeology.** The trigger *condition* is in
  the script header, but the checklist lives only in the campaign plan on
  the analysis branch — which the maintainer who eventually lifts the pin
  will not be reading. Put a five-line "How to retire this" block in the
  script header naming all four files and the trigger.
- **`:41` — `cd ""` is a bash no-op returning 0**, so `|| exit 1` never
  fires if `git rev-parse --show-toplevel` fails; the `[ -f pixi.lock ]`
  backstop is what actually fails closed. Behavior is right, but not via
  the guard the code appears to rely on. Use
  `root=$(git rev-parse --show-toplevel) && [ -n "$root" ] && cd "$root" || exit 1`.
- **The escape hatch is invisible at the moment it is needed**: on the
  wedged-network path the script prints only pixi's error. Print the
  `SKIP=pixi-lock-check` hint when `rc != 0` and the lock is unchanged and
  v6 — the one path where the developer's next move is otherwise
  `--no-verify`.
- **The docs fix is right but not self-enforcing.** Already-onboarded
  developers will not re-run `pre-commit install`.
  `default_install_hook_types: [pre-commit, pre-push]` in the config makes
  plain `pre-commit install` do the right thing and single-sources the fact.
- **At pre-push the guard checks the working tree, not the pushed ref**
  (`always_run`/`pass_filenames: false`), so `git push origin featureX:main`
  from another branch validates the wrong artifact. The new CI format check
  is the real backstop here — which is an argument that it earns its place
  beyond the pin.
- `.pre-commit-config.yaml:33-38` — `pixi-lock-check` lacks
  `always_run: true` while its sibling has it, so a push whose changed
  files all fall inside the repo's global `exclude` glob skips the wrapper.
  Minor (the tripwire still runs), but the asymmetry looks accidental.
- `:63` prints "has been restored" *before* the restore happens (the trap
  runs at exit); net effect correct.

## Suggested v3 order

1. **B5 and B6** — the six-line trap rewrite above. Ask the Analyst to
   amend the plan's fix sketch at the same time, since it prescribes B5's
   ordering.
2. Single-source `version: 6` / `0.67.2` and add the retirement block —
   together these are what make this guard removable in one commit when the
   facility upgrades.
3. The `cd ""` hardening, the SKIP hint on the network path, and
   `default_install_hook_types`.

**Note on scope for the Analyst:** v2 closed every v1 finding and closed
them well. Both new findings come from the *same* fix (B3's trap) and are
six lines apart. This is a converging slug, not a floundering one — but v3
is the last retry, so the trap rewrite should land with the plan
correction, not ahead of it.
