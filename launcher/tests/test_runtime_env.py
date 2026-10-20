"""matplotlib's cache kept node-local, outside $XDG_CACHE_HOME, however the launcher is started.

The analysis nodes' login hook sets XDG_CACHE_HOME=/var/tmp/xdgcache-$USER and deletes it on every login shell,
so anything the launcher keeps there is gone when a scientist opens a second terminal. The launcher therefore
chooses matplotlib's cache directory itself (launcher/runtime_env.py), before anything imports matplotlib.

U1-U3 and V3 call the deciding function with a mapping and a base under tmp_path. V1 starts the launcher by
the two commands that run (python -m launcher.new_launcher, and the installed new_launcher gui-script), because
the property that matters is the order of imports, which only a real start shows. `python launcher/new_launcher.py`
does not run at the base: launcher/launcher.py shadows the package (plan F3).
"""

import getpass
import logging
import os
import pathlib
import shutil
import stat
import subprocess
import sys
import time

import pytest

_REPO = pathlib.Path(__file__).resolve().parents[2]
_LOGGER = "launcher.runtime_env"
# Permission checks do not bind root, so the legs that rely on one are skipped there.
_NOT_ROOT = pytest.mark.skipif(os.geteuid() == 0, reason="root is not refused by file modes")


def _prepare(environ, base):
    from launcher.runtime_env import prepare_runtime_env

    return prepare_runtime_env(environ, base=str(base))


def _own_directory(base):
    return base / f"mpl-{getpass.getuser()}"


def _one_refusal(caplog, named):
    """The one line a refusal logs: a warning that names the path (or the reason), and no traceback (an expected
    refusal, not an error caught by the broad except)."""
    (record,) = [r for r in caplog.records if r.name == _LOGGER]
    assert record.levelno == logging.WARNING
    assert str(named) in record.getMessage()
    assert record.exc_info is None


# --- U1: the types table ------------------------------------------------------------------------------------


@pytest.mark.parametrize("xdg", [None, "xdgcache-u", "elsewhere"], ids=["no-xdg", "hook-xdg", "other-xdg"])
def test_an_unset_cache_directory_is_created_private_under_the_base(tmp_path, xdg):
    """U1, L2: wherever XDG_CACHE_HOME points, or with it unset, MPLCONFIGDIR becomes <base>/mpl-<user>,
    created 0700. Nothing is made under XDG_CACHE_HOME."""
    environ = {} if xdg is None else {"XDG_CACHE_HOME": str(tmp_path / xdg)}

    # A common umask, so that the mode is the call's doing: under 077 a mkdir with no mode makes 0700 as well.
    previous = os.umask(0o022)
    try:
        chosen = _prepare(environ, tmp_path)
    finally:
        os.umask(previous)

    expected = _own_directory(tmp_path)
    assert chosen == str(expected)
    assert environ.get("MPLCONFIGDIR") == str(expected)
    assert expected.is_dir() and not expected.is_symlink()
    assert stat.S_IMODE(expected.lstat().st_mode) == 0o700
    if xdg is not None:
        assert not (tmp_path / xdg).exists()


def test_my_existing_directory_is_reused_and_its_mode_left_alone(tmp_path):
    """U1: present, mine and writable: reused as it is, with no chmod."""
    existing = _own_directory(tmp_path)
    existing.mkdir()
    existing.chmod(0o750)
    environ = {}

    assert _prepare(environ, tmp_path) == str(existing)
    assert environ == {"MPLCONFIGDIR": str(existing)}
    assert stat.S_IMODE(existing.lstat().st_mode) == 0o750


def _make_file(path, _tmp_path):
    path.write_text("not a directory")


def _make_symlink_to_my_directory(path, tmp_path):
    target = tmp_path / "mine"
    target.mkdir(mode=0o700)
    path.symlink_to(target)


def _make_unwritable_directory(path, _tmp_path):
    path.mkdir(mode=0o700)
    path.chmod(0o500)


@pytest.mark.parametrize(
    "make",
    [_make_file, _make_symlink_to_my_directory, pytest.param(_make_unwritable_directory, marks=_NOT_ROOT)],
    ids=["a-file", "a-symlink", "not-writable"],
)
def test_a_path_that_is_not_my_writable_directory_is_refused(tmp_path, caplog, make):
    """U1, L4: a file, a symlink (never followed, even to a directory of mine) or a directory I cannot write:
    nothing is set, one line is logged, and the path is left exactly as it was."""
    path = _own_directory(tmp_path)
    make(path, tmp_path)
    before = path.lstat()
    environ = {"XDG_CACHE_HOME": str(tmp_path / "xdgcache-u")}

    with caplog.at_level(logging.DEBUG, logger=_LOGGER):
        assert _prepare(environ, tmp_path) is None

    assert environ == {"XDG_CACHE_HOME": str(tmp_path / "xdgcache-u")}
    _one_refusal(caplog, path)
    after = path.lstat()
    assert (after.st_mode, after.st_ino, after.st_mtime_ns) == (before.st_mode, before.st_ino, before.st_mtime_ns)


def test_another_users_directory_is_refused(tmp_path, caplog, monkeypatch):
    """U1, L4: a directory at the path owned by someone else is refused (the owner is compared with my uid)."""
    path = _own_directory(tmp_path)
    path.mkdir(mode=0o700)
    real_uid = os.getuid()
    monkeypatch.setattr(os, "getuid", lambda: real_uid + 1)
    environ = {}

    with caplog.at_level(logging.DEBUG, logger=_LOGGER):
        assert _prepare(environ, tmp_path) is None

    assert environ == {}
    _one_refusal(caplog, path)


@_NOT_ROOT
def test_a_base_that_cannot_be_written_is_refused_without_raising(tmp_path, caplog):
    """U1, L4: mkdir fails under the base: nothing is set, one line is logged, no exception leaves."""
    base = tmp_path / "read-only"
    base.mkdir()
    base.chmod(0o500)
    environ = {}
    try:
        with caplog.at_level(logging.DEBUG, logger=_LOGGER):
            assert _prepare(environ, base) is None
    finally:
        base.chmod(0o700)

    assert environ == {}
    assert not _own_directory(base).exists()
    _one_refusal(caplog, _own_directory(base))


def test_without_user_variables_the_name_comes_from_the_uid(tmp_path, monkeypatch):
    """U1: $USER (and getpass's other variables) unset: the name is the uid's, with the same outcome."""
    import pwd

    for name in ("LOGNAME", "USER", "LNAME", "USERNAME"):
        monkeypatch.delenv(name, raising=False)
    environ = {}

    expected = tmp_path / f"mpl-{pwd.getpwuid(os.getuid()).pw_name}"
    assert _prepare(environ, tmp_path) == str(expected)
    assert expected.is_dir()


def test_a_user_name_that_is_not_one_path_component_is_refused(tmp_path, caplog, monkeypatch):
    """U1, L4: a name that is not one path component ("a/b", "..") is refused rather than joined. getpass
    honours $LOGNAME first. <base>/mpl-a exists and is mine, so without the guard "a/b" would make a directory
    inside it."""
    (tmp_path / "mpl-a").mkdir(mode=0o700)
    environ = {}
    for name in ("a/b", ".."):
        monkeypatch.setenv("LOGNAME", name)
        with caplog.at_level(logging.DEBUG, logger=_LOGGER):
            assert _prepare(environ, tmp_path) is None
    assert environ == {}
    assert sorted(p.name for p in tmp_path.iterdir()) == ["mpl-a"]
    assert list((tmp_path / "mpl-a").iterdir()) == []


def test_an_unexpected_error_is_logged_with_its_traceback_and_goes_no_further(tmp_path, caplog, monkeypatch):
    """L4: the call runs while the launcher is imported, so no error may leave it. One that nothing expected is
    logged with its traceback, and nothing is set."""
    from launcher import runtime_env

    def broken(_base):
        raise RuntimeError("probe")

    monkeypatch.setattr(runtime_env, "_node_local_directory", broken)
    environ = {}
    with caplog.at_level(logging.DEBUG, logger=_LOGGER):
        assert _prepare(environ, tmp_path) is None

    assert environ == {}
    (record,) = [r for r in caplog.records if r.name == _LOGGER]
    assert record.exc_info is not None and record.exc_info[0] is RuntimeError


def test_the_directory_name_is_spelled_in_one_module():
    """U1, L5: the cache directory's name ("mpl-" and the user) appears in launcher/runtime_env.py and in no other
    launcher module, in code or in comments."""
    spelled = sorted(
        str(path.relative_to(_REPO))
        for path in (_REPO / "launcher").rglob("*.py")
        if "tests" not in path.relative_to(_REPO).parts and "mpl-" in path.read_text()
    )
    assert spelled == ["launcher/runtime_env.py"]


@pytest.mark.parametrize("error", [KeyError, OSError], ids=["no-passwd-entry", "no-user-name"])
def test_no_user_name_is_refused_without_raising(tmp_path, caplog, monkeypatch, error):
    """V3 (v2, design advisory A5): getpass.getuser() raising (KeyError from pwd on Python 3.11, OSError on 3.13):
    nothing set, nothing created, one line, no exception."""

    def no_user():
        raise error("no user name for this uid")

    monkeypatch.setattr(getpass, "getuser", no_user)
    environ = {}
    with caplog.at_level(logging.DEBUG, logger=_LOGGER):
        assert _prepare(environ, tmp_path) is None

    assert environ == {}
    assert list(tmp_path.iterdir()) == []
    _one_refusal(caplog, "no user name")


def test_a_path_that_cannot_be_inspected_is_refused_without_raising(tmp_path, caplog, monkeypatch):
    """V3 (v2, design advisory A5): os.lstat raising OSError on the path: nothing set, one line naming it, no
    exception."""
    target = str(_own_directory(tmp_path))
    real_lstat = os.lstat

    def lstat(path, *args, **kwargs):
        if os.fspath(path) == target:
            raise PermissionError(13, "Permission denied", target)
        return real_lstat(path, *args, **kwargs)

    environ = {}
    with monkeypatch.context() as patch, caplog.at_level(logging.DEBUG, logger=_LOGGER):
        patch.setattr(os, "lstat", lstat)
        assert _prepare(environ, tmp_path) is None

    assert environ == {}
    _one_refusal(caplog, target)


# --- U2: a preset MPLCONFIGDIR ---------------------------------------------------------------------------------


def test_an_empty_mplconfigdir_counts_as_unset(tmp_path):
    """U2: matplotlib treats an empty MPLCONFIGDIR as unset, and so does the launcher."""
    environ = {"MPLCONFIGDIR": ""}
    assert _prepare(environ, tmp_path) == str(_own_directory(tmp_path))
    assert environ["MPLCONFIGDIR"] == str(_own_directory(tmp_path))


def test_a_preset_mplconfigdir_is_left_alone(tmp_path):
    """U2, L1: the wrapper's or the user's choice wins. It is returned as the directory in force, the mapping is
    not touched, and neither it nor <base>/mpl-<user> is created."""
    preset = tmp_path / "preset"
    environ = {"MPLCONFIGDIR": str(preset), "XDG_CACHE_HOME": str(tmp_path / "xdgcache-u")}

    assert _prepare(environ, tmp_path) == str(preset)
    assert environ == {"MPLCONFIGDIR": str(preset), "XDG_CACHE_HOME": str(tmp_path / "xdgcache-u")}
    assert not preset.exists() and not _own_directory(tmp_path).exists()


# --- U3: the deciding module stays light ----------------------------------------------------------------------


def test_the_deciding_module_imports_neither_qt_nor_matplotlib():
    """U3, L3: importing launcher.runtime_env loads no Qt binding and no matplotlib, so calling it first can come
    before anything that fixes matplotlib's cache directory."""
    code = (
        "import sys, launcher.runtime_env\n"
        "heavy = {'matplotlib', 'qtpy', 'PyQt5', 'PyQt6', 'PySide2', 'PySide6'}\n"
        "print(sorted(heavy & {name.split('.')[0] for name in sys.modules}))\n"
    )
    result = subprocess.run([sys.executable, "-c", code], cwd=_REPO, capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().splitlines()[-1] == "[]"


# --- V1: the two start paths that run, by their actual commands -----------------------------------------------

# Test-only, on the child's PYTHONPATH: at interpreter start, before the launcher is imported, it points the
# deciding module's base at the test's scratch directory. runtime_env.py reads no switch for this.
_SITECUSTOMIZE = """import launcher.runtime_env
launcher.runtime_env.NODE_LOCAL_BASE = {base!r}
"""
_GUI_SCRIPT = shutil.which("new_launcher")
_START_COMMANDS = [
    pytest.param([sys.executable, "-m", "launcher.new_launcher"], id="python-m"),
    pytest.param(
        [_GUI_SCRIPT], id="gui-script",
        marks=pytest.mark.skipif(_GUI_SCRIPT is None, reason="no new_launcher gui-script installed in this environment"),
    ),
]


def _font_cache(directory):
    """matplotlib's font cache in `directory`: its sign of having chosen that directory."""
    return sorted(directory.glob("fontlist-*.json")) if directory.is_dir() else []


@pytest.mark.parametrize("command", _START_COMMANDS)
@pytest.mark.parametrize("state", ["hook-xdg", "preset", "no-xdg"])
def test_every_start_path_that_runs_keeps_the_cache_out_of_the_hooks_directory(tmp_path, command, state):
    """V1 (v2), L1-L3: the launcher started by its actual command, offscreen, and stopped once matplotlib has
    written its font cache somewhere. Then, on the filesystem:
    - under the hook's XDG_CACHE_HOME: <base>/mpl-<user> holds the cache, 0700, and the hook's directory has none.
      A call made in main(), after the imports, would still create <base>/mpl-<user>, but only after matplotlib
      had written under XDG;
    - with MPLCONFIGDIR preset: the preset holds it, and <base>/mpl-<user> is not created;
    - with no XDG_CACHE_HOME: <base>/mpl-<user> holds it, and <home>/.cache/matplotlib has none."""
    assert command[1:] == ["-m", "launcher.new_launcher"] or os.path.basename(command[0]) == "new_launcher"
    base, home, site = tmp_path / "var-tmp", tmp_path / "home", tmp_path / "site"
    for directory in (base, home, site):
        directory.mkdir()
    (site / "sitecustomize.py").write_text(_SITECUSTOMIZE.format(base=str(base)))
    hook, preset, mine = tmp_path / "xdgcache-u", tmp_path / "preset", _own_directory(base)
    env = {k: v for k, v in os.environ.items() if k not in ("MPLCONFIGDIR", "XDG_CACHE_HOME", "PYTHONPATH")}
    env.update(HOME=str(home), XDG_CONFIG_HOME=str(tmp_path / "config"), QT_QPA_PLATFORM="offscreen",
               PYTHONPATH=str(site))
    if state != "no-xdg":
        env["XDG_CACHE_HOME"] = str(hook)
    if state == "preset":
        env["MPLCONFIGDIR"] = str(preset)
    places = [mine, hook / "matplotlib", home / ".cache" / "matplotlib", preset]

    with open(tmp_path / "stderr.txt", "wb") as stderr:
        child = subprocess.Popen(command, cwd=_REPO, env=env, stdout=subprocess.DEVNULL, stderr=stderr)
        try:
            deadline = time.monotonic() + 90
            while time.monotonic() < deadline and child.poll() is None and not any(map(_font_cache, places)):
                time.sleep(0.2)
            exited = child.poll()
        finally:
            child.terminate()
            try:
                child.wait(timeout=30)
            except subprocess.TimeoutExpired:
                child.kill()
                child.wait()

    log = (tmp_path / "stderr.txt").read_text(errors="replace")[-2000:]
    assert exited is None, f"the launcher exited ({exited}) before matplotlib wrote its cache:\n{log}"
    if state == "preset":
        assert _font_cache(preset) and not mine.exists()
    else:
        assert _font_cache(mine), [str(path) for place in places for path in _font_cache(place)]
        assert stat.S_IMODE(mine.lstat().st_mode) == 0o700
        assert not (hook / "matplotlib").exists()
        assert not (home / ".cache" / "matplotlib").exists()
