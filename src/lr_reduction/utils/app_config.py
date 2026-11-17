"""
Application configuration utilities.

Configuration is managed by the `neutrons_standard` package, which provides a centralized way to load and access
application settings from a YAML configuration file.
"""

import logging
from pathlib import Path
from typing import Any, Literal, get_args, get_origin

from neutrons_standard import Config

logger = logging.getLogger(__name__)


def _matches(value: Any, expected_type: Any) -> bool:
    """`isinstance` that also accepts `Literal[...]` and rejects `bool` where `int` is expected."""
    if get_origin(expected_type) is Literal:
        return value in get_args(expected_type)
    if expected_type is int and isinstance(value, bool):  # bool is a subclass of int
        return False
    return isinstance(value, expected_type)


def validate_config():
    """
    Validate the application configuration.

    This function checks that all required configuration parameters are present and valid. If any required parameters
    are missing or invalid, it raises a ValueError with a descriptive message.

    Raises:
        ValueError: If any required configuration parameters are missing or invalid.
    """

    logger.info("Validating application configuration...")

    params_type_map = {
        "diagnostic_plots": bool,
        "rq4": bool,
        "logging": dict,
        "logging.level": Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        "logging.rotating_log_files": bool,
        "logging.num_backups": int,
        "logging.log_file": str,
        "logging.mantid_log_level": Literal["DEBUG", "INFORMATION", "NOTICE", "WARNING", "ERROR", "CRITICAL"],
    }

    invalid_params = []
    missing_params = []
    for param, expected_type in params_type_map.items():
        parent = param.rpartition(".")[0]
        if parent in missing_params or parent in invalid_params:
            continue  # parent already reported; Config can't look up keys inside a non-dict value
        if not Config.exists(param):
            missing_params.append(param)
        elif not _matches(Config[param], expected_type):
            invalid_params.append(param)
    bad_params = missing_params + invalid_params

    # check log file path exists
    if (
        all(p not in bad_params for p in ["logging", "logging.log_file", "logging.rotating_log_files"])
        and Config["logging.rotating_log_files"] is True
    ):
        log_file_path = Path(Config["logging.log_file"])
        if not log_file_path.parent.exists():
            logger.warning(
                f"Log file directory does not exist: {log_file_path.parent}. Logging will fall back to stdout."
            )

    error_message = ""
    if missing_params:
        error_message += f"Missing required configuration parameters: {missing_params}. "

    if invalid_params:
        error_message += f"Invalid configuration parameters: {invalid_params}. "

    if error_message:
        logger.error(error_message)
        raise ValueError(error_message)
