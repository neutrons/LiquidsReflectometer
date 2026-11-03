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
"""

import io
import re
import shutil
import subprocess
import tokenize
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
TREES = ("src", "launcher", "tests", "scripts")
#: A noqa directive, as ruff reads one: `noqa` right after a `#` (or after `ruff:` / `flake8:` for a file's
#: exemption), then a colon, whitespace or the end. Prose that mentions noqa is not one.
DIRECTIVE = re.compile(r"#\s*(?:(?:ruff|flake8)\s*:\s*)?noqa(?=[\s:]|$)(?P<rest>.*)", re.IGNORECASE)
SCOPED = re.compile(r"\s*:\s*[A-Z]+[0-9]+(?:\s*,\s*[A-Z]+[0-9]+)*")


def ruff(*args, stdin=None):
    """`ruff check` from the repository's root, so that its pyproject.toml configuration applies. pixi's default
    environment, which runs the tests, carries ruff."""
    executable = shutil.which("ruff")
    assert executable, "ruff is not on PATH: run the tests in pixi's default environment, which carries it"
    return subprocess.run([executable, "check", "--no-fix", "--output-format", "concise", *args], cwd=ROOT,
                          input=stdin, capture_output=True, text=True, check=False, timeout=300)


def hook_args():
    """The args of the pre-commit `ruff-check` hook, the one place this repository runs ruff."""
    config = yaml.safe_load((ROOT / ".pre-commit-config.yaml").read_text(encoding="utf-8"))
    hooks = [hook for repo in config["repos"] for hook in repo.get("hooks", []) if hook.get("id") == "ruff-check"]
    assert len(hooks) == 1, hooks
    return [str(arg) for arg in hooks[0].get("args", [])]


def python_files():
    for tree in TREES:
        for path in sorted((ROOT / tree).rglob("*.py")):
            if (ROOT / "tests" / "data") not in path.parents:
                yield path


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
    result = ruff("--extend-select", "RUF100", *TREES)
    unused = [line for line in result.stdout.splitlines() if "RUF100" in line]
    assert unused == [], "\n".join(unused)
    assert result.returncode == 0, result.stdout + result.stderr


def test_no_noqa_tag_where_ruff_does_not_look():
    """U2: a tag in a file ruff does not lint silences nothing. The files ruff lints are ruff's own list
    (--show-files), so this follows pyproject.toml's exclusions."""
    listed = ruff("--show-files", *TREES)
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
