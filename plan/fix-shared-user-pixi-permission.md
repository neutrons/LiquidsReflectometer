# Plan: Fix Shared-User Pixi Permission Issue in lr_reduction

## Context

When multiple developers share a single project clone on a shared analysis
system (e.g., one group checkout at `/SNS/REF_L/shared/dev/lr_reduction/`),
running `pixi install` as one user creates `.pixi/envs/` and
`.pixi/solve-group-envs/` with that user's ownership. Other users then get
permission errors on `pixi run {command}` because pixi writes to those
directories even when the environment is already installed.

The nsd-app-wrap project shows two complementary mechanisms:
1. **`--frozen --manifest-path`** on `pixi run` (used in `pixi_launch()`) —
   prevents lock-file updates; bypasses local `.pixi/` discovery entirely.
2. **`detached-environments` in `.pixi/config.toml`** — moves the conda
   environment out of the project's `.pixi/envs/` into a managed location so
   pixi never needs to write to the project directory during normal use.

The `.pixi/.gitignore` that pixi auto-generates already contains:
```
*
!config.toml
```
meaning `.pixi/config.toml` IS tracked by git if it exists; everything else in
`.pixi/` is ignored. So committing a `config.toml` requires no root-`.gitignore`
changes (though adding a root-level pattern mirrors nsd-app-wrap best practice).

---

## Scenario Breakdown

### A. Individual developer on their own machine
**Status: no problem.** `pixi install` → `.pixi/envs/` owned by that one user.
`pixi run {task}` works fine. No changes needed.

### B. Shared team clone (collaborative development) ← _the question_
Multiple users share one checkout directory. One user's `pixi install`
breaks others. **Three options; see below.**

### C. Shared system deployment (production/QA)
Already solved by nsd-app-wrap's `pixi_launch` pattern in the existing
`nr_launcher.sh` / `nr_launcher_qa.sh` scripts. Environments live at
`/usr/local/pixi/{env}/`, managed by admins.

---

## Options for Scenario B (Shared Team Clone)

### Option 1 — `detached-environments = true` (per-user home)  ✅ RECOMMENDED

Commit `.pixi/config.toml` to the project:

```toml
# .pixi/config.toml
[workspace]
detached-environments = true
```

**Effect:** When any user runs `pixi install`, their conda environment is stored
in their own `~/.pixi/envs/{hash}/` (where the hash is derived from the
project manifest path). Each user has an **independent private copy** of the
environment. The shared project clone's `.pixi/` directory contains only this
`config.toml` and the auto-generated housekeeping files — nothing user-owned.

**Workflow:**
```bash
# Each developer, once after cloning or after pixi.lock changes:
pixi install
# Then normal usage — no flags needed:
pixi run test-reduction
pixi run docs-build
```

**Is it feasible?** Yes — completely. Pixi handles it transparently.
**Is it recommended?** Yes — it's the simplest fix with the widest benefit.
No admin involvement, no discipline required, works identically on any machine
(shared cluster or personal laptop).

**Trade-offs:**
- Each user needs ~2–3 GB of disk in their home directory for the extracted
  environment (large due to Mantid). However, the **conda package cache**
  (tarballs) lives at `~/.pixi/cache/` by default and can be pointed at a
  shared location via `PIXI_CACHE_DIR=/shared/pixi-cache` to avoid
  re-downloading packages.
- When `pixi.lock` is updated by one developer and pushed, all others must
  re-run `pixi install`. This is normal git workflow behavior.
- `pixi.lock` itself is still at the project root (version-controlled); only
  one person should run `pixi update`/`pixi add` at a time.

---

### Option 2 — Group-writable `.pixi/` + `--frozen` discipline

Admin sets up `.pixi/` once with group permissions:
```bash
pixi install                         # as admin/owner
chgrp -R <devgroup> .pixi/
chmod -R 2775 .pixi/                 # setgid bit preserves group on new files
```

All developers then use `--frozen` on every run:
```bash
pixi run --frozen test-reduction
```

**Is it feasible?** Yes, technically.
**Is it recommended?** No — it is fragile. A single user running `pixi run`
*without* `--frozen` (habit, CI script, IDE integration) will attempt to
rewrite solver state and may corrupt `.pixi/` for everyone. It requires
ongoing discipline and re-setup after every `pixi.lock` update.

---

### Option 3 — Shared group path (`detached-environments = "/shared/path/"`)

Create a LOCAL (not committed) `.pixi/config.toml` at the shared checkout:
```toml
[workspace]
detached-environments = "/SNS/REF_L/shared/pixi-envs/lr_reduction_dev"
```
Admin installs once; all users share one environment via that group-writable
path.

**Is it feasible?** Yes, but this is essentially the same model as Scenario C
(deployment) applied to development.
**Is it recommended?** Only if the team explicitly wants a single shared
environment under admin control. For development it creates coordination
overhead (one bad `pixi add` affects everyone).

---

## Recommended Implementation (Options 1 + C)

### Step 1: Commit `.pixi/config.toml` with `detached-environments = true`

File: `lr_reduction/.pixi/config.toml` (new — pixi's auto `.pixi/.gitignore`
already allows this file to be tracked)

```toml
[workspace]
detached-environments = true
```

This fixes Scenario B immediately and is invisible to developers on personal
machines (their environment just lands in `~/.pixi/envs/` instead of the
project's `.pixi/`).

### Step 2: Add `launcher/nr_launcher_next.sh`

Completes the environment-naming coverage (prod=`refred`, qa=`lr_reduction_qa`,
next=missing):

```bash
#!/usr/bin/bash

# import library to do the real work
# shellcheck disable=SC1091
. /bin/nsd-app-wrap.sh

# put together arguments - env, application, argv
args=("lr_reduction_next" "python /SNS/REF_L/shared/launcher/launcher.py" "$@")
# launch the tool
( cd /SNS/REF_L/shared/launcher && pixi_launch "${args[@]}" )
```

### Step 3: Update `CLAUDE.md` — add Shared Development and Deployment sections

```markdown
## Pixi Environment Management

### Shared development clone (multi-user)
`.pixi/config.toml` is committed with `detached-environments = true`. This
means each developer's conda environment is installed to their own
`~/.pixi/envs/{hash}/` — never into the shared project directory. Normal
`pixi run {task}` works for all users without permission conflicts.

Each developer runs `pixi install` once after cloning, and again after any
`pixi.lock` update is pulled.

To share the conda *package cache* (tarballs, ~several GB) and avoid
re-downloading for each user, set in your shell profile:
```bash
export PIXI_CACHE_DIR=/SNS/REF_L/shared/pixi-cache   # or any shared path
```

### Deployment (analysis.sns.gov)
Production and QA environments are centrally managed at `/usr/local/pixi/`:

| Branch | Environment path                      | Launcher script          |
|--------|---------------------------------------|--------------------------|
| next   | `/usr/local/pixi/lr_reduction_next/`  | `nr_launcher_next.sh`    |
| qa     | `/usr/local/pixi/lr_reduction_qa/`    | `nr_launcher_qa.sh`      |
| prod   | `/usr/local/pixi/refred/`             | `nr_launcher.sh`         |

Admin deploys a deployment `.pixi/config.toml` (not in git) pointing to the
system path, installs once with proper group permissions, then users access
the environment read-only via the launcher scripts (which use nsd-app-wrap's
`pixi_launch --frozen --manifest-path`).
```

---

## Files to Modify

| File | Action |
|------|--------|
| `lr_reduction/.pixi/config.toml` | Create (new, will be tracked by git) |
| `lr_reduction/CLAUDE.md` | Add Pixi Environment Management section |
| `lr_reduction/launcher/nr_launcher_next.sh` | Create (new) |

---

## Verification

1. **Scenario B (shared clone)**: Two users on the same machine share a clone.
   - User A: `pixi install` — confirm env appears in `~/.pixi/envs/`, NOT in
     project `.pixi/envs/`
   - User B: `pixi install` — confirm env appears in User B's `~/.pixi/envs/`
   - User B: `pixi run test-reduction` — should succeed with no permission error
2. **Scenario A (local dev)**: `pixi install` on a personal machine — confirm
   env goes to `~/.pixi/envs/` (same behavior, no regression)
3. **nr_launcher_next.sh**: source it, verify it calls
   `pixi_launch lr_reduction_next python ...`
