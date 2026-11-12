# Application Configuration

The `lr_reduction` application uses a YAML-based configuration system with environment layering, deep merging, and user configuration support. This page describes the application-level configuration, which is separate from the reduction workflow configuration described in {ref}`reduction_config`.

## Default Configuration

The default application configuration is provided by `resources/application.yml`:

```yaml
application:
  name: lr_reduction

diagnostic_plots: true
rq4: false

logging:
  level: INFO
  rotating_log_files: true
  num_backups: 15
```

This file includes:

- **`application.name`** — the application identifier
- **`diagnostic_plots`** — enable or disable diagnostic plot generation (default: `true`)
- **`rq4`** — enable or disable RQ4-specific behavior (default: `false`)
- **`logging.level`** — the logging level for the application (default: `INFO`)
- **`logging.rotating_log_files`** — whether to use rotating log files (default: `true`)
- **`logging.num_backups`** — the number of backup log files to keep (default: `15`)

## Environment-Specific Configurations

In addition to the base `application.yml`, you can provide environment-specific configuration files (e.g., `production.yml`, `qa.yml`, `development.yml`). When an environment is specified, the base configuration and environment overrides are merged recursively, with environment values taking precedence over base values.

Specify the environment via:
- The `env` environment variable (e.g., `env=production`)
- The `environment` property in YAML
- User configuration (see below)

## User Application Configuration

You can override the default application configuration by creating a user-specific configuration file at:

```
~/.config/lr_reduction/lr_reduction-user.yml
```

The user configuration file:
- **Auto-loads** if present and no explicit environment is specified
- **Deep merges** with the base `application.yml`, recursively combining dictionaries while replacing scalar values
- **Is archived and regenerated** if the `application.version` in your user config doesn't match the installed version
- **Preserves backups** in `~/.config/lr_reduction/` with versioned names like `lr_reduction-user-X.Y.Z-2026-03-04T16:10:34.817Z.yml.bak`

The number of archived user configs kept is controlled by `application.user_yml_archive_limit` (default: 5, or `-1` for unlimited).

## Configuration Access in Code

Within the application, configuration values are accessed using the `Config` singleton with dot-delimited notation:

```python
from lr_reduction import app_config

# Access values
app_config["logging.level"] # Dot-delimited lookup

# Check if a key exists
app_config.exists("optional.key")
```

Nested values can be referenced using dot notation, and token substitution (e.g., `${key.subkey}`) is supported to avoid duplication.

## Configuration Backups

The system automatically maintains backups:
- **Latest backup:** `~/.config/lr_reduction/application.yml.bak` — updated on every config reload
- **User config archives:** Versioned backups of outdated user configs, with the newest `application.user_yml_archive_limit` kept
