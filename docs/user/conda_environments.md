# Conda Environments

Three conda environments are available in the analysis nodes, beamline machines, as well as the
jupyter notebook severs. On a terminal:

```bash
$ conda activate <environment>
```

where `<environment>` is one of `lr_reduction`, `lr_reduction-qa`, and `lr_reduction-dev`

## lr_reduction Environment

Activates the latest stable release of `lr_reduction`. Typically users will reduce their data in this environment.

## lr_reduction-qa Environment

Activates a release-candidate environment.
Instrument scientists and computational instrument scientists will carry out testing on this environment
to prevent bugs being introduce in the next stable release.

## lr_reduction-dev Environment

Activates the environment corresponding to the latest changes in the source code.
Instrument scientists and computational instrument scientists will test the latest changes to `lr_reduction` in this
environment.
