#!/usr/bin/env bash
# Push-time guard for the committed pixi.lock.
#
# Condition being defended (self-contained, so this script needs no external
# document to be understood): analysis.sns.gov runs a pixi older than 0.68 and
# cannot read lock-format v7. pixi >= 0.68 writes only v7, and no newer pixi can
# write v6 — so a lock regenerated on an unpinned machine is unusable at the
# facility. The version test itself lives in scripts/check_lock_format.sh.
#
# CONTRACT: a push never carries a non-v6 lock, and never carries dependency
# drift. "The working tree is never dirty" is NOT the contract — ordinary pixi
# activity restamps the self-package record, and the standing discipline is to
# restore rather than stage that.
#
# KNOWN LIMIT: at the pre-push stage these hooks inspect the working tree, not
# the ref being pushed, so `git push <remote> other-branch:main` validates the
# wrong artifact. The CI format check is the structural backstop for that case,
# which is most of why it exists.
#
# RETIREMENT — when analysis.sns.gov's pixi is upgraded to >= 0.68, delete this
# guard in ONE commit touching:
#     scripts/pixi_lock_check.sh, scripts/check_lock_format.sh,
#     .pre-commit-config.yaml, .github/workflows/test_and_deploy.yml
# (the campaign control plane's amendment 14 tracks the upgrade trigger).
#
# What is known about `pixi lock --check`, stated exactly:
#
#   * It is NOT read-only. MEASURED on pixi 0.67.2 (2026-08-19, twice, plus
#     repeated live observation during the campaign): on an up-to-date lock it
#     still re-stamps the `- pypi: ./` self-package record — versioningit
#     derives that version from `git describe`, so it changes after every
#     commit — and exits 0.
#   * Whether its exit code reports genuine dependency drift is UNVERIFIED
#     here. Three probe attempts hit three different walls: an unsolvable
#     package (non-zero for the wrong reason), a package already present in the
#     lock as a transitive dependency (no drift to detect), and — for both a
#     real package and a hand-tampered lock — `failed to fetch conda-pypi
#     mapping from remote source`, which fails every re-solve in this
#     environment. Do not restate this as measured until someone measures it.
#
# That uncertainty is precisely why the verdict below comes from the resulting
# file rather than from the exit code.
set -u

# git hooks run with a minimal PATH; pixi's default install location is not on
# it. (Dropping this line is what made an earlier revision unrunnable as a hook.)
PATH="$HOME/.pixi/bin:$PATH"

# Always operate on the repository root: pre-commit may invoke this from a
# subdirectory, and pixi searches parent directories for a manifest — so a
# relative path would have it rewrite the real lock while restoring nothing.
# The emptiness test is load-bearing: `cd ""` succeeds as a no-op, so a failed
# rev-parse would otherwise slip through the `|| exit 1`.
root=$(git rev-parse --show-toplevel) && [ -n "$root" ] && cd "$root" || exit 1
[ -f pixi.lock ] || { echo "pixi-lock-check: no pixi.lock at the repository root (failing closed)"; exit 1; }

snap=$(mktemp) || exit 1
# Snapshot BEFORE arming the trap: an EXIT trap installed over an empty mktemp
# file would copy that emptiness onto pixi.lock if this cp failed — and a
# truncated lock whose first line still reads `version: 6` passes every
# tripwire this guard ships.
cp pixi.lock "$snap" || { rm -f "$snap"; exit 1; }
restore() { cp "$snap" pixi.lock 2>/dev/null; rm -f "$snap"; }
# A bash INT/TERM handler returns to the script unless it exits; without the
# explicit exits below, an interrupt would restore, delete the snapshot, and
# then run the classifier against a file that no longer exists.
trap restore EXIT
trap 'restore; trap - EXIT; exit 130' INT
trap 'restore; trap - EXIT; exit 143' TERM

pixi lock --check
rc=$?

# The format check comes first and is decisive: a v7 lock is unusable at the
# facility no matter what else the check reported.
if ! head -1 pixi.lock | grep -qx 'version: 6'; then
    cat >&2 <<'MSG'
pixi-lock-check: your pixi wrote lock-format v7 — analysis.sns.gov's pixi
(< 0.68) cannot read it. Do NOT commit this lock. Pin your toolchain instead:

    pixi self-update --version 0.67.2

Your working-tree pixi.lock will be restored to the committed v6 version.
MSG
    exit 1
fi

# Drift classification from content, not from rc (see the header).
#
# Ignored inside the `- pypi: ./` stanza: exactly the three things the check
# itself rewrites on an up-to-date lock. MEASURED on pixi 0.67.2 (2026-08-19),
# committed lock vs post-check, whole-stanza diff:
#
#     - version: 2.10.0.dev20260610180700+gf5b8e45   (versioningit restamp)
#     + version: 2.10.0.dev20260819073821+g3c1808b
#     - sha256:  c96e2ca…                            (follows the version)
#     + sha256:  6781a23…
#     - editable: true                               (dropped; the installed
#                                                     env stays editable)
#
# Anything else in the stanza — requires_dist, requires_python, name — and
# anything anywhere else in the file is real drift. Narrowing this set further
# (for instance treating the `editable` removal as drift) makes the guard fail
# on every ordinary push, since the check performs that removal every time.
strip_self_package_stamp() {
    awk '
        /^- pypi: \.\/$/ { inself = 1; print; next }
        inself && /^- / { inself = 0 }
        inself && /^  (version|sha256): / { next }
        inself && /^  editable: true$/ { next }
        { print }
    ' "$1"
}

if ! cmp -s pixi.lock "$snap"; then
    echo "pixi-lock-check: pixi.lock will be restored (the check mutated it)"
    hunk=$(diff <(strip_self_package_stamp "$snap") <(strip_self_package_stamp pixi.lock) || true)
    if [ -n "$hunk" ]; then
        echo "pixi-lock-check: pixi.lock does not match pyproject.toml — regenerate it UNDER A pixi <= 0.67.x and commit the result:"
        echo "$hunk"
        exit 1
    fi
elif [ "$rc" -ne 0 ]; then
    # The check failed without touching the lock: an unreachable index or a
    # wedged network, not drift. Name the targeted skip, because the
    # alternative a blocked developer reaches for is --no-verify, which would
    # also disable the format tripwire.
    echo "pixi-lock-check: the check failed but pixi.lock is unchanged and v6 (network or index problem?)." >&2
    echo "If you must push now, skip this one hook rather than all of them:  SKIP=pixi-lock-check git push" >&2
fi

exit $rc
