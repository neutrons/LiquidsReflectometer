# Logging in LR Reduction

LR Reduction uses the python `logging` module to provide a flexible logging system.
The logging configuration is specified in the application configuration file, which can be found at `resources/application.yml`. Users can override the default logging settings by creating a user-specific configuration file at `~/.config/lr_reduction/lr_reduction-user.yml`.

It also uses Mantid's `log_to_python` function to redirect Mantid's logging output to the Python logging system. The Mantid log level can also be configured in the application configuration file, separately from the standard Python logging configuration.

See {ref}`app_config` for more information on the application configuration system.

**NOTES**:

- Logging is only configured when the application is run via the command line interface (CLI).
  If you are using the API directly (e.g. from a Python script or within Mantid workbench),
  you will need to configure logging yourself. For example:

  ```python
  import logging

  from lr_reduction.utils.logging import setup_logging

  setup_logging()
  logger = logging.getLogger(__name__)

  logger.info("This is an info message.")
  ```
