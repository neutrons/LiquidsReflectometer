import neutrons_standard

try:
    from lr_reduction._version import __version__
except ImportError:
    __version__ = "unknown"


### Initialize the Neutrons Standard package for lr_reduction
neutrons_standard.init("lr_reduction")

# Load the application configuration from the YAML file and re-export it as `app_config` for use throughout the package.
# This allows for centralized configuration management and easy access to application settings.
from lr_reduction.utils.logging import setup_logging # noqa: E402, I001
from neutrons_standard import Config  # noqa: E402, I001

app_config = Config

setup_logging(
    level=Config["logging.level"],
    reset_handlers=True,
    rotate=Config["logging.rotating_log_files"],
    num_backups=Config["logging.num_backups"],
)


#######################################################################################
__all__ = ["__version__", "app_config"]
