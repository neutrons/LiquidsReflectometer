# Integrator review gate — `pixi-lock-format-guard` v1 REJECTED

**Gate: green** (`test-launcher` 2 passed, `test-reduction` 107 passed,
exit 0, 510.60 s). This rejection is entirely from the review gate.

**Disposition: blocking — a flagged departure from the plan's advisory
declaration.** The plan declares `design-reviewer (advisory)`, and I have
honored advisory declarations throughout this campaign where the findings
were quality opinions. I am departing here because this slug's *entire
deliverable is a guard*, and three of the four findings describe the guard
**permitting, causing, or instructing the exact harm it exists to
prevent**. All four were verified in source by the Integrator before
rejecting. Retry budget: attempt 1 of N = 3.

**Analyst override remains cheap and I will not re-litigate it:** re-tag
`qa/pixi-lock-format-guard` at an unchanged SHA and I will open the draft
PR with all of this recorded as advisory instead.

## What is right, and better than the plan — state this first

Two judgment calls here improve on the plan's own sketch, and should
survive into v2 unchanged:

- **`cp` back rather than `git checkout -- pixi.lock`** (`:48-50`). This is
  the case most likely to be got wrong, and it is right: a deliberate
  uncommitted lock edit survives the check byte-for-byte. Verified at
  review by appending a WIP stanza and running the wrapper against a
  restamping pixi — the WIP survived.
- **Deciding from file content rather than from `rc`** (`:41-51`), which is
  strictly stronger than the plan's rc-trusting design and the correct
  response to an exit code whose behaviour is unverified.
- The `awk` self-package strip is precisely scoped (removes exactly the 5
  stanza lines on the real lock; the two indented `- pypi: ./` environment
  references survive), and putting the pure-shell tripwire on **both**
  `pre-commit` and `pre-push` stages is what keeps the guard alive at all
  (see A2).

## B1 — the failure message prescribes the destructive action

`scripts/pixi_lock_check.sh:45`:

```
pixi-lock-check: pixi.lock is out of date with pyproject.toml — run 'pixi lock' and commit the result
```

On a machine with pixi ≥ 0.68 — the only machines where the v7 hazard
exists at all — `pixi lock` produces the v7 rewrite this slug exists to
prevent. The push is correctly blocked, and then the guard tells the
maintainer to go and do the thing. Reproduced at review with a ≥0.68 shim:
the wrapper blocks with `EXIT=1`, prints the `version: 6` → `version: 7`
hunk, and issues that instruction.

**Fix:** after `pixi lock --check`, test `head -1 pixi.lock` and branch to
a v7-specific message *before* the generic drift path — "your pixi wrote
lock-format v7; do NOT commit it; pin pixi < 0.68" plus the pin command.

## B2 — the missing-lock path fails OPEN and writes a 0-byte lock

`:24-25` does `cp pixi.lock "$snap"` unchecked, and the script has **no
repo-root assertion** (verified: zero occurrences of `show-toplevel` or a
`-f pixi.lock` test). Run anywhere `pixi.lock` is not in the cwd and the
chain is: `cp` fails → `awk` fails → both stripped streams are empty →
`diff -q` reports them identical → `drift=0` → `cp "$snap" pixi.lock`
writes the **empty** snapshot → `exit $rc`.

Reproduced at review from a subdirectory: `EXIT=0`, and a **0-byte
`pixi.lock` created in the cwd**. This is the one path where the guard
silently passes. With real pixi it is worse: pixi searches parent
directories for the manifest, so it rewrites the *real* `../pixi.lock`
while the wrapper restores nothing and the push proceeds.

**Fix:** `cd "$(git rev-parse --show-toplevel)" || exit 1` and
`[ -f pixi.lock ] || { echo ...; exit 1; }` at the top.

## B3 — no `trap`, so the file's own contract is false under interrupt

`:20` states: *"The lock is restored either way, so a push never dirties
the tree."* There is no `trap` in the file (verified: zero). Reproduced by
killing the hook mid-check: the tree is left **converted to v7**, `git
status` shows `M pixi.lock`, and the only backup is an unnamed
`/tmp/tmp.XXXXXXXX` the developer has no way to identify.

This is not a corner case. `pixi lock --check` performs a network solve,
and this commit's own body reports that solve failing on network in this
environment today — so a Ctrl-C or a dropped terminal during it is an
everyday event, and the outcome is the exact state the slug exists to
prevent, caused by the slug.

**Fix:** `trap 'cp "$snap" pixi.lock 2>/dev/null; rm -f "$snap"' EXIT INT TERM`.

## B4 — the script asserts a measurement the same commit says was never made

`scripts/pixi_lock_check.sh:12-15` states as fact:

> *"Measured on pixi 0.67.2 (2026-08-19): with a dependency added to
> pyproject.toml and absent from the lock, `pixi lock --check` UPDATED the
> lock and exited 0."*

The body of commit `2b04921` states:

> *"I have NOT demonstrated that. Three attempts each hit a different wall
> … the rc's behaviour under real drift remains unknown here."*

Both landed in the same commit. That measurement is the **sole
justification for the entire content-classifier design**, so this is not a
comment nit: either the claim is unverified and must not be stated as
measured, or the measurement happened and the commit body is wrong — and
per the plan's own pathological row it would then need recording through
the review loop, which has not happened.

This is the campaign's signature defect shape once more — *a documented
guarantee the code does not provide* — here in the form of a documented
**measurement that was not taken**. Reconcile it before merge.

## Advisory — carry into v2

- **The CI change is a pin, not a check.** The workflow has no `pre-commit`
  step at all, so this guard has zero presence in GitHub Actions; pinning
  `setup-pixi` to `v0.67.2` means a v7 lock breaks CI as an inscrutable
  parse failure rather than as "your lock is v7". A two-line step placed
  **before** the first `setup-pixi` closes it properly:
  `head -1 pixi.lock | grep -qx 'version: 6' || { echo "..."; exit 1; }`.
- **The documented onboarding never installs the pre-push hook.**
  `docs/developer/developer.rst:106` says `pre-commit install`, not
  `--hook-type pre-push`, and this clone has no `.git/hooks/pre-push` — so
  the wrapper is dead code for anyone following the docs. The dual-stage
  tripwire rescues the design; the doc line should become
  `pre-commit install --hook-type pre-commit --hook-type pre-push`.
- **The classifier's blind spot is the whole stanza, not just the restamp.**
  Any real change confined to the `- pypi: ./` record is invisible —
  reproduced by adding a `requires_dist` entry inside it (`EXIT=0`,
  silently passed). Exposure is low today (deps are conda-side), but
  `editable: true` lives in that stanza and the plan itself flags its loss
  as a semantic concern. Strip only the `version:`/`sha256:` lines *within*
  the stanza rather than the stanza itself.
- **A network failure wedges the push with no message and no escape.**
  Fail-closed is right for *format*, but the pure-shell tripwire already
  covers that offline; a failed network solve blocking a push with only
  pixi's stderr invites `--no-verify`, which disables the tripwire too.
  Document `SKIP=pixi-lock-check git push` in `developer.rst`.
- **The restore is silent.** The plan's
  `echo "restored pixi.lock (the check mutated it)"` was dropped in
  implementation, so a developer on an unpinned pixi never learns their
  toolchain is mutating the lock — and the commit body cites that line as
  verification output it cannot have produced.
- **Hard-coded values are not single-sourced.** `version: 6` lives only in
  `.pre-commit-config.yaml:29` — the wrapper never checks the format at all,
  which is the mechanical root of B1. `v0.67.2` appears three times in the
  workflow with no comment. Both failure messages cite "charter amendment
  14", a document not present in this repository, so the citation is a dead
  reference for any reader who is not the author. Name the condition inline
  at each pin site.
- Nits: `diff` is computed twice (`:43`, `:46`); `head -20` will truncate a
  real v7 diff misleadingly; the script is unlinted
  (`.pre-commit-config.yaml:5` excludes `scripts/.*`, and there is no
  shellcheck hook).

## A design-level observation worth acting on

During this review, `pixi.lock` in the working tree was restamped **with no
push involved** — most likely an ordinary `pixi install` or direnv's
`watch_file pixi.lock` activation. A pre-push-only wrapper structurally
cannot keep the tree clean, because the mutation arrives through everyday
environment activation. The dual-stage tripwire already acknowledges this;
the wrapper's own framing ("a push never dirties the tree") does not.
Worth deciding explicitly in v2 whether the goal is *"pushes never carry a
v7 lock"* (achievable) or *"the tree is never dirty"* (not achievable at
this layer).

## Suggested v2 order

1. **B2** and **B3** — the fail-open path and the missing `trap`. Together
   they are ~4 lines and they are the two ways this guard currently permits
   or causes the harm.
2. **B1** — branch on `head -1 pixi.lock` before the generic drift message.
3. **B4** — reconcile the comment with what was actually measured; if the
   measurement stands, record it, and if it does not, say so in the comment
   and keep the content-classifier design on its *other* justification
   (an unverified rc is reason enough not to trust it).
4. The CI format check and the `developer.rst` hook-type line — without
   them the guard does not reach the machines that most need it.
