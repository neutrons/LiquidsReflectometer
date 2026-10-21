"""noqa hygiene (noqa-sweep) — every `noqa` tag names the rule it silences, and none is dead.

A tag that silences nothing hides nothing today, but it will hide the next real finding on its line. A tag spelled
without a colon (`# noqa E722`) is not code-scoped at all: ruff reads it as a blanket `# noqa` and silences every rule
on its line, whatever code it names. These tests are the guard:
- ruff's RUF100 reports no tag in the files it lints (none names a rule that is not enabled there, none silences
  nothing);
- no tag sits in a file ruff does not lint (`scripts/` is excluded in pyproject.toml): such a tag is dead by
  construction;
- the pre-commit `ruff-check` hook, where this repository runs ruff, adds RUF100 to the selection, so the hook fails when
  a tag goes dead. pyproject.toml's selection is left alone: that file is also pixi's self-package manifest, and any
  edit to it leaves pixi.lock out of date;
- no tag is a blanket: each names its codes after a colon (`# noqa: BLE001 -- the reason`).

Frozen files are not ours to lint: the legacy-reduction files kept at upstream's content (PR #219 review, 2026-10-20), and
the template system's two (2026-10-21), are excluded by the `ruff-check` hook's `exclude`, so these tests skip what that `exclude` names, from the same config (one
source), and a guard checks that it names exactly those files and nothing else.
"""

import io
import re
import shutil
import subprocess
import tokenize
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
TREES = ("src", "launcher", "tests", "scripts")
#: A noqa directive, as ruff reads one: `noqa` right after a `#` (or after `ruff:` / `flake8:` for a file's
#: exemption), then a colon, whitespace or the end. Prose that mentions noqa is not one.
DIRECTIVE = re.compile(r"#\s*(?:(?:ruff|flake8)\s*:\s*)?noqa(?=[\s:]|$)(?P<rest>.*)", re.IGNORECASE)
SCOPED = re.compile(r"\s*:\s*[A-Z]+[0-9]+(?:\s*,\s*[A-Z]+[0-9]+)*")
#: The legacy-reduction files kept at upstream's content (PR #219 review, 2026-10-20), and the template system's two
#: (`template.py`, `new_reduction_from_template.py`; 2026-10-21). The hook's `exclude` is what the tests read; this list is
#: only the guard's authority on what that `exclude` may name.
FROZEN = (
    "launcher/apps/roi_selector.py",
    "launcher/scripts/off_spec.py",
    "launcher/scripts/quick_reduce.py",
    "launcher/scripts/refracted_beam.py",
    "launcher/scripts/xrr_reduction.py",
    "src/lr_reduction/background.py",
    "src/lr_reduction/dead_time_correction.py",
    "src/lr_reduction/new_reduction_from_template.py",
    "src/lr_reduction/new_reduction_template_reader.py",
    "src/lr_reduction/output.py",
    "src/lr_reduction/peak_finding.py",
    "src/lr_reduction/reduction_template_reader.py",
    "src/lr_reduction/template.py",
    "src/lr_reduction/time_resolved.py",
)


def ruff(*args, stdin=None):
    """`ruff check` from the repository's root, so that its pyproject.toml configuration applies. pixi's default
    environment, which runs the tests, carries ruff."""
    executable = shutil.which("ruff")
    assert executable, "ruff is not on PATH: run the tests in pixi's default environment, which carries it"
    return subprocess.run([executable, "check", "--no-fix", "--output-format", "concise", *args], cwd=ROOT,
                          input=stdin, capture_output=True, text=True, check=False, timeout=300)


def ruff_hook():
    """The pre-commit `ruff-check` hook, the one place this repository runs ruff."""
    config = yaml.safe_load((ROOT / ".pre-commit-config.yaml").read_text(encoding="utf-8"))
    hooks = [hook for repo in config["repos"] for hook in repo.get("hooks", []) if hook.get("id") == "ruff-check"]
    assert len(hooks) == 1, hooks
    return hooks[0]


def hook_args():
    """The args of the pre-commit `ruff-check` hook."""
    return [str(arg) for arg in ruff_hook().get("args", [])]


def hook_excluded_paths():
    """The file paths the hook's `exclude` names literally, one per line of its verbose regex (`tests/data/.*`, a
    pattern, is not a path)."""
    paths = []
    for line in str(ruff_hook().get("exclude", "")).splitlines():
        entry = line.strip().rstrip("|")
        if not entry or entry.startswith("(?x)") or entry.startswith(")") or ".*" in entry:
            continue
        paths.append(entry.replace("\\.", "."))
    return paths


def python_files():
    """The Python files ruff is ours to run on: not test data, not what the hook excludes."""
    excluded = re.compile(str(ruff_hook().get("exclude", "")) or r"(?!)")
    for tree in TREES:
        for path in sorted((ROOT / tree).rglob("*.py")):
            if (ROOT / "tests" / "data") not in path.parents and not excluded.search(path.relative_to(ROOT).as_posix()):
                yield path


def ruff_on_ours(*args):
    """`ruff check` on the trees with the hook's excluded files left out, as the hook leaves them out."""
    frozen = hook_excluded_paths()
    return ruff(*args, *(["--extend-exclude", ",".join(frozen)] if frozen else []), *TREES)


def noqa_directives(path):
    """(line, comment, the text after `noqa`) for each directive in the file's comments; strings are not read."""
    source = path.read_text(encoding="utf-8")
    for token in tokenize.generate_tokens(io.StringIO(source).readline):
        if token.type == tokenize.COMMENT:
            for match in DIRECTIVE.finditer(token.string):
                yield token.start[0], token.string, match.group("rest")


def test_ruff_reports_no_unused_noqa_tag():
    """U2: with the hook's selection (the repository's, plus RUF100), RUF100 reports no tag in src, launcher, tests or
    scripts, and the tree is otherwise clean."""
    result = ruff_on_ours("--extend-select", "RUF100")
    unused = [line for line in result.stdout.splitlines() if "RUF100" in line]
    assert unused == [], "\n".join(unused)
    assert result.returncode == 0, result.stdout + result.stderr


def test_no_noqa_tag_where_ruff_does_not_look():
    """U2: a tag in a file ruff does not lint silences nothing. The files ruff lints are ruff's own list
    (--show-files), so this follows pyproject.toml's exclusions and the hook's (passed as --extend-exclude); the
    files the hook excludes are skipped, as they are not ours to lint."""
    listed = ruff_on_ours("--show-files")
    assert listed.returncode == 0, listed.stdout + listed.stderr
    linted = {Path(line.strip()).resolve() for line in listed.stdout.splitlines() if line.strip()}
    unlinted = [f"{path.relative_to(ROOT)}:{line}: {comment}" for path in python_files() if path.resolve() not in linted
                for line, comment, _rest in noqa_directives(path)]
    assert unlinted == [], "\n".join(unlinted)


def test_the_pre_commit_hook_selects_ruf100():
    """U2b: the pre-commit hook adds RUF100 to the selection (`--extend-select RUF100` in its args), and ruff given the
    hook's args reports a dead tag: here one naming F401, which is enabled but does not fire on its line, on stdin
    under a path in src/. Without it, the hook would not see a tag go dead."""
    args = hook_args()
    selected = [args[i + 1] for i, arg in enumerate(args[:-1]) if arg == "--extend-select"]
    selected += [arg.split("=", 1)[1] for arg in args if arg.startswith("--extend-select=")]
    assert any("RUF100" in value.split(",") for value in selected), args
    as_the_hook = [arg for arg in args if arg not in ("--fix", "--exit-non-zero-on-fix")]
    result = ruff(*as_the_hook, "--stdin-filename", "src/lr_reduction/noqa_probe.py", "-",
                  stdin="VALUE = 1  # noqa: F401\n")
    assert "RUF100" in result.stdout, result.stdout + result.stderr


def test_every_noqa_tag_names_its_codes():
    """U3: no tag is a blanket. `# noqa`, `# noqa E722` (no colon) and `# noqa:` without a code each silence every
    rule on their line; a tag names its codes after a colon."""
    blankets = [f"{path.relative_to(ROOT)}:{line}: {comment}" for path in python_files()
                for line, comment, rest in noqa_directives(path) if not SCOPED.match(rest)]
    assert blankets == [], "\n".join(blankets)


def test_the_hook_excludes_exactly_the_frozen_files():
    """The skip above follows the hook's `exclude`, so the `exclude` must not grow. Read two ways, as pre-commit reads it
    and line by line: compiled, it matches every frozen file and no other tracked Python file outside tests/data/ (a
    pattern-shaped widening — a `.*` line, an entry on the `(?x)` line, an alternative after the close, a wider
    tests/data/ — fails here); literally, every path it names exists and is a frozen file, and every frozen file is
    named (a typo'd extra path, which hides nothing, fails here)."""
    if not (ROOT / ".git").exists():
        pytest.skip("not a git checkout: the tracked files cannot be listed")
    excluded = re.compile(str(ruff_hook().get("exclude", "")))
    assert [path for path in FROZEN if not excluded.search(path)] == []
    tracked = subprocess.run(["git", "ls-files", "*.py"], cwd=ROOT, capture_output=True, text=True, check=True, timeout=60)
    hidden = [path for path in tracked.stdout.split()
              if excluded.search(path) and path not in FROZEN and not path.startswith("tests/data/")]
    assert hidden == [], hidden
    named = hook_excluded_paths()
    missing = [path for path in named if not (ROOT / path).is_file()]
    assert missing == [], missing
    assert sorted(named) == sorted(FROZEN), sorted(set(named) ^ set(FROZEN))
