# Plan: Fix Shared-User Pixi Permission Issue in lr_reduction

## Context

When the lr_reduction project is deployed to a shared analysis system (e.g.,
`/SNS/REF_L/shared/launcher/`), running `pixi install` as one user creates
`.pixi/envs/` with that user's ownership/permissions. Subsequent `pixi run
{command}` calls by other users fail because pixi attempts to write to the
project's `.pixi/` directory (lock state, cache timestamps, etc.) and is
denied.

The nsd-app-wrap project (used system-wide at analysis.sns.gov) has solved this
via two complementary mechanisms drawn from the pattern in
`nsd-app-wrap/src/nsd-app-wrap.sh`:

1. **`detached-environments` in `.pixi/config.toml`**: moves the conda
   environment out of the project's `.pixi/envs/` into an admin-managed system
   path (e.g., `/usr/local/pixi/lr_reduction_next/`). Pixi no longer needs to
   write into the project directory at all.

2. **`--frozen --manifest-path /central/path/`**: the key flags on every
   `pixi run` call via `pixi_launch()` in `nsd-app-wrap.sh`. `--frozen`
   prevents lock-file updates; `--manifest-path` bypasses local `.pixi/`
   discovery and points straight to the pre-installed central environment.

Together, these mean a user running `pixi run {task}` (via a wrapper or
directly) never writes to the project checkout's `.pixi/`.

nsd-app-wrap's `.gitignore` reflects the intent:
```
.pixi/*
!.pixi/config.toml
```
— gitignore everything in `.pixi/` *except* the config, so a deployment-site
`config.toml` can be tracked or preserved without polluting git history.

---

## Proposed Changes

### 1. Update `.gitignore` (lr_reduction root)

Add the nsd-app-wrap pattern so `.pixi/config.toml` is not swallowed by the
auto-generated `.pixi/.gitignore`:

```
# pixi environments
.pixi/*
!.pixi/config.toml
```

File: `lr_reduction/.gitignore`

### 2. Create `.pixi/config.toml` template (NOT committed to repo)

Document (in CLAUDE.md, step 3) that the deployment admin creates this file
**at the deployment site only**:

```toml
# .pixi/config.toml — created by admin at deployment time; NOT in git
[workspace]
detached-environments = "/usr/local/pixi/lr_reduction_<tier>"
```

Where `<tier>` is `next`, `qa`, or `prod` (matching the environment naming
already used by the launch scripts).

**Why not commit this?** Committing a hard-coded system path would break
developer machines where `/usr/local/pixi/` may not exist or be writable.
The `.gitignore` exception allows the deployment site to maintain this file
locally without accidentally committing it.

### 3. Add `launcher/nr_launcher_next.sh`

The production and QA launchers already use the correct nsd-app-wrap pattern.
A `next`/development launcher is missing. Create it:

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

This mirrors `nr_launcher_qa.sh` (which uses `lr_reduction_qa`) for the
`next` branch deployment at `/usr/local/pixi/lr_reduction_next/`.

### 4. Update `CLAUDE.md` — add Deployment section

Add a section documenting the complete multi-user deployment workflow:

```markdown
## Shared-System Deployment (analysis.sns.gov)

Pixi environments on shared analysis systems are managed centrally at
`/usr/local/pixi/` to avoid permission conflicts between users.

### Environment naming convention
| Branch | Pixi environment path                      |
|--------|--------------------------------------------|
| next   | `/usr/local/pixi/lr_reduction_next/`       |
| qa     | `/usr/local/pixi/lr_reduction_qa/`         |
| prod   | `/usr/local/pixi/refred/` (shared w/ REFL)|

### Admin-only: initial deployment
1. Deploy the project checkout to `/usr/local/pixi/lr_reduction_<tier>/`
2. Create `.pixi/config.toml` there:
   ```toml
   [workspace]
   detached-environments = "/usr/local/pixi/lr_reduction_<tier>"
   ```
3. Run `pixi install` as admin → environment installs into the detached path
   with admin-controlled permissions (e.g., group-readable).
4. Set permissions: `chmod -R 755 /usr/local/pixi/lr_reduction_<tier>/.pixi/`

### User access — run pixi tasks without write access
```bash
pixi run --frozen --manifest-path /usr/local/pixi/lr_reduction_<tier>/ test-reduction
```
Or via the launcher scripts (which use `nsd-app-wrap`'s `pixi_launch`
automatically):
```bash
./launcher/nr_launcher_next.sh       # next
./launcher/nr_launcher_qa.sh         # qa
./launcher/nr_launcher.sh            # prod (refred env)
```

### Rule: never run bare `pixi install` in the shared checkout
Always use the `--frozen --manifest-path` pattern or go through the launcher
wrappers. Plain `pixi install` creates/modifies `.pixi/` with the current
user's ownership and breaks other users.
```

---

## Files to Modify

| File | Action |
|------|--------|
| `lr_reduction/.gitignore` | Add `.pixi/*` + `!.pixi/config.toml` |
| `lr_reduction/CLAUDE.md` | Add Deployment section |
| `lr_reduction/launcher/nr_launcher_next.sh` | Create (new file) |

## Verification

1. On the shared deployment system:
   - Admin creates `.pixi/config.toml` with `detached-environments`
   - Admin runs `pixi install` — confirm environment appears at central path,
     NOT in project `.pixi/envs/`
   - Second user runs `pixi run --frozen --manifest-path /central/path/ test-reduction`
     — should succeed without write errors
2. Locally (developer machine): confirm `pixi run test-reduction` still works
   normally (no `.pixi/config.toml` → default local behavior unchanged)
3. `nr_launcher_next.sh` smoke test: source it and verify it calls
   `pixi_launch lr_reduction_next python ...`
