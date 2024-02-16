# standard imports
from contextlib import contextmanager
from copy import deepcopy
from typing import Union

# third-party libraries
from mantid.kernel import ConfigService


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
