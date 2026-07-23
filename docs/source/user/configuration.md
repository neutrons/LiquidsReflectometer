# Configuration of Reduction Workflow

This describes the new, configuration-driven reduction workflow, which is separate from the
legacy RefRed/XML-template workflow described in {ref}`workflow`. The two workflows have separate
entry points and coexist during a transition period; this page covers only the new one.

## Global vs. per-run parameters

A configuration file describes a set of runs to reduce — not necessarily all from the same
sequence — and distinguishes:

- **Global parameters** — the same for every run in the configuration (e.g. `instrument`,
  `assembly` stitching options, `output` settings, `corrections` such as dead-time and
  emission-time).
- **Per-run parameters** — the scientific parameters that can vary per run. ROI-selection
  parameters (peak/background pixel ranges) are shared by direct beam and reflected runs alike and
  modeled by `RunParameters` ({class}`lr_reduction.models.config.RunParameters`); parameters that
  only make sense for a reflected run (Q binning, acceptance window, etc.) are modeled by
  `ReflectedRunParameters` ({class}`lr_reduction.models.config.ReflectedRunParameters`), which
  extends `RunParameters`.

Per-run parameters use an inheritance model: a global `defaults` block sets the values
shared by every run, and each run may override any subset of those fields (`None` means "inherit
from defaults"). A typical configuration is therefore short: most runs only override
`direct_beam`, `peak`, and maybe `background`. `ReductionConfig.effective()` resolves the final,
merged parameters for a given run.

## Composite direct beams

Composite direct beams are defined once, in their own named key space (`direct_beams`), separate
from the `sequence_number` key space used for reflected runs. A reflected run references its
direct beam by name, and the same direct beam may be referenced by more than one reflected run:

```yaml
instrument: BL4B

defaults:
  q_binning:
    q_min: 0.005
    q_max: 0.5
    q_step: 0.01
    method: constantQ

direct_beams:
  db_8mm:
    db_runs: [12345, 12346, 12347]
  db_20mm:
    db_runs: [12350, 12351]

runs:
  - sequence_number: 1
    direct_beam: db_8mm
    run_number: 12360
    peak: { min: 140, max: 150 }
  - sequence_number: 2
    direct_beam: db_20mm
    run_number: 12361
    peak: { min: 138, max: 152 }
  - sequence_number: 3
    direct_beam: db_8mm       # sharing db_8mm with run 1
    run_number: 12362
```

Referential integrity (every `direct_beam` reference resolves to a defined composite direct beam)
and range checks (e.g. pixel ranges, Q ranges) are enforced by Pydantic at load time, with clear
error messages.

A reflected run's data source is either a single `run_number`, or a `source_runs` list naming
multiple runs to sum (exactly one of the two must be set); an optional `filter` restricts the
source run(s) by time and/or a log value. Both are recorded as part of the configuration and
therefore carried into provenance and re-applied on re-reduction.

## Loading a configuration

`lr_reduction.io.config_loader.ConfigLoader` loads a native YAML configuration file and resolves
it to a `ReductionConfig`:

```python
from lr_reduction.io.config_loader import ConfigLoader

config = ConfigLoader().load("config.yaml")
```

### Single-run extraction

Selecting one run from a configuration yields a valid, self-contained one-run
configuration — there is no separate single-run format:

```python
single_run_config = config.for_run(sequence_number=2)
```

### Overrides

Caller-supplied overrides can be applied on top of a loaded configuration for a single invocation,
without ever modifying the configuration file on disk. There are two ways to do this:

Directly setting an attribute is the simplest option for a one-off override, and is still
validated — every model in the schema sets `validate_assignment=True`, so an out-of-range or
wrong-type value is rejected immediately, the same as at load time:

```python
config.assembly.q_norm = 0.02          # validated against AssemblyConfig's own constraints
config.runs[0].peak = {"min": 100, "max": 110}   # a plain dict is coerced into a PixelRange
```

The caveat: this only re-validates the model whose attribute you set, not a *different* model
that depends on it — for example, setting `config.runs[0].direct_beam` to an unknown name does
not re-run `ReductionConfig`'s referential-integrity check, since that check lives on the parent
model. For overrides that could affect cross-model checks like referential integrity, or to
express a whole batch of overrides as one dict, use `apply_overrides`, which re-validates the
entire configuration tree:

```python
from lr_reduction.models.config import apply_overrides

config = apply_overrides(config, {"assembly": {"q_norm": 0.02}})
```
