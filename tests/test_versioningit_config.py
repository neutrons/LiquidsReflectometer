"""Guard the versioningit configuration against the 2026-06-05 autodeploy regression.

Commit 5f1c70d changed the dev-build version format from `{next_version}.dev{distance}`
(-> 2.10.0.dev80) to `{version}+{distance}.{vcs}{rev}` (-> 2.9.0rc1+86.g57f255c) and the
next-version method from `minor` to `minor-release`. Both lowered the version: a
`2.9.0rc1+N` string sorts BELOW the already-published `2.10.0.dev80`, so the deployment
environment (`lr_reduction = {version="*"}` from the `neutrons/label/exp` channel) kept
selecting the OLD build and the new code never deployed.

These tests pin the invariant that fixed it: dev builds derive from {next_version} (the
bumped version) via the `minor` method, while retaining the git hash for provenance.
See tasking/plan/lr_reduction_exp/autodeploy-investigation.md for the full analysis.
"""
import os

try:
    import tomllib
except ImportError:  # Python < 3.11
    import tomli as tomllib

_HERE = os.path.dirname(__file__)
_PYPROJECT = os.path.join(_HERE, os.pardir, "pyproject.toml")


def _versioningit_cfg():
    with open(_PYPROJECT, "rb") as fh:
        return tomllib.load(fh)["tool"]["versioningit"]


def test_distance_format_uses_next_version_base():
    """Dev builds must be derived from {next_version}, not the raw {version} tag.

    Using the raw {version} (e.g. 2.9.0rc1) regresses the version below a prior
    publish; {next_version} (e.g. 2.10.0) keeps it monotonically increasing.
    """
    fmt = _versioningit_cfg()["format"]
    for key in ("distance", "distance-dirty"):
        assert "{next_version}" in fmt[key], (
            f"[tool.versioningit.format].{key} = {fmt[key]!r} must use {{next_version}} "
            "so dev builds sort ABOVE prior releases. Using the raw {version} tag base "
            "caused the 2026-06-05 autodeploy regression."
        )
        assert "{version}+" not in fmt[key], (
            f"[tool.versioningit.format].{key} = {fmt[key]!r} uses the raw {{version}} "
            "tag base directly; that is the regression. Use {next_version}."
        )


def test_next_version_method_bumps_minor():
    """`minor` bumps 2.9.0rc1 -> 2.10.0; `minor-release` would only strip the rc -> 2.9.0
    (lower). The monotonic-version fix requires the bumping method."""
    assert _versioningit_cfg().get("next-version", {}).get("method") == "minor", (
        "next-version method must be 'minor' (bumps the minor and drops the prerelease); "
        "'minor-release' only strips the prerelease and sorts lower."
    )


def test_format_retains_git_hash():
    """Keep the provenance commit 5f1c70d intended: the git revision in the version."""
    assert "{vcs}{rev}" in _versioningit_cfg()["format"]["distance"], (
        "the dev-build version should embed the git hash ({vcs}{rev}) for traceability."
    )


def test_computed_version_does_not_regress():
    """Belt-and-suspenders: if versioningit is importable, the computed version's release
    must be >= (2, 10, 0) so it sorts at/above the published 2.10.0.dev80 baseline."""
    import pytest

    versioningit = pytest.importorskip("versioningit")
    from packaging.version import Version

    raw = versioningit.get_version(project_dir=os.path.join(_HERE, os.pardir))
    v = Version(raw.replace(".dirty", ""))
    assert v.release >= (2, 10, 0), (
        f"computed version {raw!r} (release {v.release}) regressed below the published "
        "2.10.0.dev80 baseline; dev builds would be un-selectable by `version=\"*\"`."
    )


def test_match_glob_present():
    """The agentic workflow pushes lightweight coordination tags (qa/*, review/*,
    triage/*, analysis/*). Without a release-tag allow-list, `git describe` could feed
    one of those to versioningit -> InvalidVersionError at build/env-prep time. A
    positive match glob keeps only release tags as version candidates. See
    tasking/plan/lr_reduction_exp/findings.md F3."""
    vcs = _versioningit_cfg()["vcs"]
    assert "match" in vcs, (
        "[tool.versioningit.vcs].match is required so versioningit ignores the workflow's "
        "qa/*, review/*, triage/*, analysis/* coordination tags."
    )
    assert any("v" in g and "[0-9]" in g for g in vcs["match"]), (
        f"unexpected match glob {vcs['match']!r}; expected a release-tag glob like ['v[0-9]*']."
    )


def test_match_glob_finds_a_release_tag():
    """Sanity: the configured glob matches at least one tag reachable from HEAD, so
    versioning still resolves (catches a future glob typo)."""
    import re
    import subprocess

    globs = _versioningit_cfg()["vcs"]["match"]
    cmd = ["git", "describe", "--tags", *(f"--match={g}" for g in globs), "HEAD"]
    out = subprocess.run(cmd, capture_output=True, text=True, cwd=os.path.join(_HERE, os.pardir))
    assert out.returncode == 0, f"git describe with match={globs} failed: {out.stderr!r}"
    assert re.match(r"^v\d", out.stdout.strip()), (
        f"git describe returned {out.stdout.strip()!r}; expected a v-prefixed release tag."
    )
