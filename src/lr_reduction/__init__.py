import neutrons_standard

try:
    from lr_reduction._version import __version__
except ImportError:
    __version__ = "unknown"


### Initialize the Neutrons Standard package for lr_reduction
neutrons_standard.init("lr_reduction")


### Configure logging

# Load the application configuration from the YAML file and re-export it as `app_config` for use throughout the package.
# This allows for centralized configuration management and easy access to application settings.
from mantid.utils.logging import log_to_python  # noqa: E402
from neutrons_standard import Config  # noqa: E402

from lr_reduction.utils.logging import setup_logging  # noqa: E402

app_config = Config  # Re-export the loaded configuration for use in other modules of the package

setup_logging(
    level=Config["logging.level"],
    reset_handlers=True,
    rotate=Config["logging.rotating_log_files"],
    num_backups=Config["logging.num_backups"],
)

# Configure Mantid logging to integrate with the Python logging system,
# ensuring that logs from Mantid are captured and formatted consistently with the rest of the application.
log_to_python(
    level="notice",
    pattern="%s | %t",
)

#######################################################################################
__all__ = ["__version__", "app_config"]
