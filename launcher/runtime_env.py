"""Where the launcher keeps matplotlib's cache: node-local, outside $XDG_CACHE_HOME.

The analysis nodes' login hook (/etc/profile.d/fontcache-workaround.sh) sets
XDG_CACHE_HOME=/var/tmp/xdgcache-$USER and deletes that directory on every login
shell, so anything kept there is gone when a scientist opens a second terminal.
matplotlib keeps its cache (the font list) under $XDG_CACHE_HOME unless
MPLCONFIGDIR is set, and fixes the directory the first time anything loads its
font manager. The launcher's tabs import matplotlib when they are imported, so
the choice is made by this module, which imports neither Qt nor matplotlib, and
is called from the top of launcher/new_launcher.py before those imports.

The directory is <base>/mpl-<user>, the deploy wrapper's choice too. It is
created 0700 when absent, and used only when it is a directory owned by this
user that this user can write. It is never followed through a symlink, chmod'ed
or deleted. When it cannot be used, nothing is set, one line is logged, and
matplotlib applies its own fallback.
"""

import getpass
import logging
import os
import stat

#: Node-local and kept across logins; the login hook's own directory is under it too.
NODE_LOCAL_BASE = "/var/tmp"
CACHE_DIR_PREFIX = "mpl-"

_log = logging.getLogger(__name__)


def prepare_runtime_env(environ=os.environ, base=None):
    """
    Set MPLCONFIGDIR to a node-local per-user directory outside $XDG_CACHE_HOME, unless it is set already.

    Parameters
    ----------
    environ : MutableMapping[str, str]
        The environment to read and update; the process's own by default.
    base : str, optional
        The directory to make it in; NODE_LOCAL_BASE when None.

    Returns
    -------
    str or None
        The cache directory in force: MPLCONFIGDIR as preset, or as set here. None when nothing could be set;
        the reason is logged.
    """
    # The wrapper's or the user's choice wins. matplotlib reads an empty value as unset, and so does this.
    if environ.get("MPLCONFIGDIR"):
        return environ["MPLCONFIGDIR"]
    try:
        path, reason = _node_local_directory(NODE_LOCAL_BASE if base is None else base)
    except Exception:
        # This runs while the launcher is imported: nothing here may stop it starting.
        _log.warning("matplotlib cache: MPLCONFIGDIR not set after an unexpected error", exc_info=True)
        return None
    if path is None:
        _log.warning("matplotlib cache: %s; MPLCONFIGDIR not set, so matplotlib chooses", reason)
        return None
    environ["MPLCONFIGDIR"] = path
    return path


def _node_local_directory(base):
    """Return (<base>/mpl-<user>, None), created 0700 if absent, or (None, why it cannot be used)."""
    try:
        user = getpass.getuser()
    except (KeyError, OSError) as error:
        return None, f"no user name ({error})"
    if user in ("", os.curdir, os.pardir) or os.sep in user or (os.altsep and os.altsep in user):
        return None, f"the user name {user!r} is not one path component"
    path = os.path.join(base, CACHE_DIR_PREFIX + user)
    try:
        os.mkdir(path, 0o700)
    except FileExistsError:
        pass  # judged below as found: never followed, chmod'ed or replaced
    except OSError as error:
        return None, f"cannot create {path} ({error.strerror})"
    try:
        info = os.lstat(path)
    except OSError as error:
        return None, f"cannot inspect {path} ({error.strerror})"
    if not stat.S_ISDIR(info.st_mode):
        return None, f"{path} is not a directory"
    if info.st_uid != os.getuid():
        return None, f"{path} belongs to uid {info.st_uid}"
    if not os.access(path, os.W_OK | os.X_OK):
        return None, f"{path} is not writable"
    return path, None
