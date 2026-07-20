# Configuration of Reduction Workflow

This describes the new, configuration-driven reduction workflow, which is separate from the
legacy RefRed/XML-template workflow described in {ref}`workflow`. The two workflows have separate
entry points and coexist during a transition period; this page covers only the new one.

## Sequence-wide vs. per-run parameters

A configuration file describes an entire sequence of runs (identified by `sequence_id`), and
distinguishes:

- **Sequence-wide parameters** — the same for every run in the sequence (e.g. `instrument`,
  `assembly` stitching options, `output` settings).
- **Per-run parameters** — the scientific parameters that can vary per run (peak/background pixel
  ranges, Q binning, corrections, etc.), modeled by `RunParameters`
  ({class}`lr_reduction.models.config.RunParameters`).

Per-run parameters use an inheritance model: a sequence-wide `defaults` block sets the values
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
sequence_id: 12340
experiment_id: IPTS-12345

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

`lr_reduction.io.config_loader.ConfigLoader` is designed to accept either form below, both
resolving to a `ReductionConfig`, so a caller need not distinguish between them:

```python
from lr_reduction.io.config_loader import ConfigLoader

# (a) a native YAML configuration file
config = ConfigLoader().load("config.yaml")

# (b) a prior ORSO output file, whose header carries the recorded configuration (see Provenance)
config = ConfigLoader().load("reduced_output.ort")
```

A legacy XML template can also be loaded this way (`ConfigLoader().load("template.xml")`) for
backward compatibility during the coexistence period, but that path is deprecated: the new
workflow's own input forms are just YAML and a prior ORSO output, per (a)/(b) above. Converting a
legacy template to YAML first (see Migrating a legacy XML template, below) is the sanctioned way
to bring it into the new workflow.

### Single-run extraction

Selecting one run from a sequence configuration yields a valid, self-contained one-run
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

## Provenance

The complete configuration used to reduce a run is embedded in the reduced output's ORSO header
(mirroring `ReductionResult.reduction_config`, {class}`lr_reduction.models.results.ReductionResult`),
so a reduced file is self-describing and can be handed back to `ConfigLoader` to reproduce or
re-run the reduction (§4.2.3). Reading it back (`lr_reduction.io.orso.read_config`) is implemented;
writing it out in the first place depends on the reduction/assembly operations that produce a
`ReductionResult`, which are still stubs, so `lr_reduction.io.orso.write` is too.

## Migrating a legacy XML template

The new workflow does not read legacy RefRed XML templates directly (§3.3.10). A standalone
utility, kept separate from `ConfigLoader` so the new workflow's own loader never reads XML,
converts an XML reduction template into this YAML format — injecting schema defaults for
parameters the XML doesn't describe, and wrapping each legacy single direct-beam run as a one-run
composite direct beam (§3.3.10.1–.2):

```bash
$ lr-xml-to-yaml template.xml config.yaml
```

or, from a repository checkout without installing the package first,
`python -m lr_reduction.io.xml_to_yaml template.xml config.yaml`. Both the command-line entry
point and the `lr_reduction.io.xml_to_yaml.convert()` library call are implemented.
