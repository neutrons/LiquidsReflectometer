#!/usr/bin/env bash
# The single source of the lock-format condition.
#
# analysis.sns.gov runs a pixi older than 0.68 and cannot read lock-format v7.
# pixi >= 0.68 writes only v7, and no newer pixi can write v6 — so a lock
# regenerated on an unpinned machine is unusable at the facility.
#
# Invoked from the pixi-lock-format pre-commit/pre-push hook and from every CI
# job that touches the lock, so the version test and its rationale live here
# once. See scripts/pixi_lock_check.sh for the retirement checklist.
set -u

lock="${1:-pixi.lock}"

if [ ! -f "$lock" ]; then
    echo "check_lock_format: $lock not found" >&2
    exit 1
fi

if ! head -1 "$lock" | grep -qx 'version: 6'; then
    echo "$lock is not lock-format v6 — analysis.sns.gov's pixi (< 0.68) cannot read v7." >&2
    echo "Regenerate it under a pixi <= 0.67.x and commit the result." >&2
    exit 1
fi
