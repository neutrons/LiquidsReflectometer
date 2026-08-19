#!/usr/bin/env bash
# Guard the committed pixi.lock at push time (charter amendment 14, campaign
# exp-settings-roi).
#
# Two behaviours of `pixi lock --check` make it unusable on its own:
#
#   1. It is not read-only. Even when the lock is current it re-stamps the
#      editable self-package record (versioningit derives that version from
#      `git describe`, so it changes on every commit), and on pixi >= 0.68 it
#      rewrites the whole file to lock-format v7, which analysis.sns.gov's pixi
#      cannot read.
#   2. Its exit code does not report drift. Measured on pixi 0.67.2
#      (2026-08-19): with a dependency added to pyproject.toml and absent from
#      the lock, `pixi lock --check` UPDATED the lock and exited 0 — so a stale
#      lock would sail straight through a gate that trusts its rc.
#
# So: snapshot, run the check, and decide from the resulting file rather than
# from the exit code. Changes confined to the self-package stanza are the
# benign restamp and are ignored; anything else is real drift and blocks the
# push. The lock is restored either way, so a push never dirties the tree.
set -u
PATH="$HOME/.pixi/bin:$PATH"

snap=$(mktemp) || exit 1
cp pixi.lock "$snap"

# Drop the `- pypi: ./` self-package stanza, whose version/sha256 restamp is
# expected on every commit and says nothing about dependency drift.
strip_self_package() {
    awk '
        /^- pypi: \.\/$/ { skip = 1; next }
        skip && /^- / { skip = 0 }
        skip { next }
        { print }
    ' "$1"
}

pixi lock --check
rc=$?

drift=0
if ! cmp -s pixi.lock "$snap"; then
    if ! diff -q <(strip_self_package "$snap") <(strip_self_package pixi.lock) >/dev/null; then
        drift=1
        echo "pixi-lock-check: pixi.lock is out of date with pyproject.toml — run 'pixi lock' and commit the result"
        diff <(strip_self_package "$snap") <(strip_self_package pixi.lock) | head -20
    fi
    cp "$snap" pixi.lock
    # cp back rather than `git checkout --`: restores only what the check
    # changed, preserving any deliberate uncommitted lock edit.
fi
rm -f "$snap"

if [ "$drift" -ne 0 ]; then
    exit 1
fi
exit $rc
