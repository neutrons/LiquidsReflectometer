- should the config loader:
  - search a default directory for config or orso files?
  - accept a directory path as an argument and look for config or orso files there?
  - accept a file path directly to a config or orso file?

- what does the cd attenuation table look like?

- what are partials?
  - load one or all partials from a directory?

ORSO:

- according to the [orso spec](https://www.reflectometry.org/advanced_and_expert_level/file_format), it is recommended to use the `.ort` extension - i would suggest we follow that recommendation rather than use `.txt` for partials. we can use `.txt` for the legacy partials, but we should use `.ort` for new partials.

- facility = sns or ornl?

- what fields do we want to populate?
  - mainly in Experiment, Sample, InstrumentSettings, Measurement

- which values should be single values vs. ranges?

- are liquid experiments polarized or unpolarized? if polarized, what are the polarization states?

```
Meta = {
  "wl_min": 13.05504715524623,
  "wl_max": 16.39347348797596,
  "q_min": 0.008027891007325118,
  "q_max": 0.010080776946106923,
  "theta": 0.010472986057602888,
  "start_time": "2023-01-16T10:47:39.627455667",
  "experiment": "IPTS-29196",
  "run_number": "201282",
  "run_title": "Expt 8 Cu-B BF4 noEtOH Full OCV 1-201282-1.",
  "norm_run": 201043,
  "time": "Fri Feb 20 11:31:18 2026",
  "dq0": 0,
  "dq_over_q": 0.022671471200417816,
  "sequence_number": 1,
  "sequence_id": 201282,
  "q_summing": false,
  "specular_pixel": 141.0,
  "use_functional_bck": false,
  "scaling_factors": {"a": 1,
  "err_a": 0,
  "b": 0,
  "err_b": 0},
  "tof_weighted": false,
  "bck_in_q": false,
  "theta_offset": 0
}
```
