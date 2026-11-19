import neutrons_standard

try:
    from lr_reduction._version import __version__
except ImportError:
    __version__ = "unknown"


### Initialize the Neutrons Standard package for lr_reduction
neutrons_standard.init("lr_reduction")

# Load the application configuration from the YAML file and re-export it as `app_config` for use throughout the package.
# This allows for centralized configuration management and easy access to application settings.
from neutrons_standard import Config  # noqa: E402

from lr_reduction.utils.app_config import validate_config  # noqa: E402

app_config = Config  # Re-export the loaded configuration for use in other modules of the package
validate_config()

#######################################################################################
__all__ = ["__version__", "app_config"]
