from contextlib import contextmanager
from copy import deepcopy
from pathlib import Path
from typing import Union

from mantid.api import Workspace
from mantid.kernel import ConfigService
from mantid.simpleapi import mtd

from lr_reduction.typing import MantidWorkspace


def mantid_algorithm_exec(algorithm_class, **kwargs):
    algorithm_instance = algorithm_class()
    assert algorithm_instance.PyInit, "str(algorithm_class) is not a Mantid Python algorithm"
    algorithm_instance.PyInit()
    for name, value in kwargs.items():
        algorithm_instance.setProperty(name, value)
    algorithm_instance.PyExec()
    if 'OutputWorkspace' in kwargs:
        return algorithm_instance.getProperty('OutputWorkspace').value


def workspace_handle(workspace: MantidWorkspace) -> Workspace:
    r"""
    Utility function to get a workspace handle from either a workspace name or a workspace object.

    Parameters
    ----------
    workspace
        Name of the workspace or the workspace object

    Returns
    -------
        The workspace object
    """
    if isinstance(workspace, str):
        return mtd[workspace]
    else:
        return workspace


@contextmanager
def amend_config(
    new_config: dict = None, data_dir: Union[str, list] = None, data_dir_insert_mode: str = "prepend"
) -> None:
    r"""
    Context manager to safely modify Mantid Configuration Service while
    the function is executed.

    Parameters
    ----------
    new_config
        (key, value) pairs to substitute in the configuration service
    data_dir
        prepend one (when passing a string) or more (when passing a list)
        directories to the list of data search directories. Alternatively, replace instead of prepend.
    data_dir_insert_mode
        How to insert the data directories. Options are: "prepend" (default) and "replace".
    """
    modified_keys = list()
    backup = dict()
    config = ConfigService.Instance()
    if new_config is not None:
        SEARCH_ARCHIVE = "datasearch.searcharchive"
        if SEARCH_ARCHIVE not in new_config:
            new_config[SEARCH_ARCHIVE] = "hfir, sns"
        DEFAULT_FACILITY = "default.facility"
        if DEFAULT_FACILITY not in new_config:
            new_config[DEFAULT_FACILITY] = "SNS"
        for key, val in new_config.items():
            backup[key] = config[key]
            config[key] = val  # config does not have an 'update' method
            modified_keys.append(key)
    if data_dir is not None:
        data_dirs = (
            [
                data_dir,
            ]
            if isinstance(data_dir, str)
            else data_dir
        )
        # make sure the data_dirs exists and are directories
        for path in data_dirs:
            if Path(path).is_dir() is False:
                raise ValueError(f"Data directory: {path} does not exist or is not a directory")
        key = "datasearch.directories"
        backup[key] = deepcopy(config[key])
        # prepend or replace our custom data directories to the list of data search directories
        if data_dir_insert_mode == "prepend":
            config.setDataSearchDirs(data_dirs + list(config.getDataSearchDirs()))
        elif data_dir_insert_mode == "replace":
            config.setDataSearchDirs(data_dirs)
        else:
            raise ValueError(f"Invalid data_dir_insert_mode: {data_dir_insert_mode}")
    try:
        yield
    finally:
        for key in modified_keys:
            config[key] = backup[key]
