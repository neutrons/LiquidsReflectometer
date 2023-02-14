## Event-based reduction for the Liquids Reflectometer

### template.py
`template.py` is used to load a template and launch the reduction for a single run.
It will also apply the scaling factors.

`process_from_template()` will load the data, but `process_from_template_ws` can
also be called on a Mantid workspace. This makes it possible to reduce event-filtered data.
`reduce_30Hz_slices` will reduce time-resolved data by splitting the run in time intervals.

### event_reduction.py
`event_reduction.py` implements the data reduction for a single Mantid workspace.
It takes the scattering and direct beam workspaces, along with all the necessary reduction
parameters. It can do specular and off-specular. It can reduced specular reflectivity in
constant-Q binning or not, and can estimate the background in both pixel-TOF space or pixel-Q space.

### reduction_template_reader.py
Implements a template reader for the xml templates used by RefRed and the automated reduction.
It can read and write templates.

### workflow.py
Abstracts out the automated reduction workflow. From a run number and an output directory,
it will call the reduction, assemble results from all the runs in the corresponding set,
write a template that can be loaded in RefRed with the proper run numbers, and writes the
final output.

### output.py
Implements the logic to handle a collection of runs that maps a single R(q) curve.
It can read and write the LR standard output files.

### To-do list
 - Abstract out distances
 - Create time-resolved config file as part of the output so that we can easily run it again
 - Add extra options to template, like q-summing, pre/post-cut, etc...
 - Launcher: Add alternate output directory
 - Launcher: Add time-resolved option to use the new reduction
 - Test bck in Q
 - Add option to skip background subtraction
 - Remove direct beam background
 - Constant-q binning: sign of delta pixel may depend on whether we reflect up or down
 - Capture tthd when calibrating for fixed tthd so that we can use the calibration if tthd is moved
 - In template.process_from_template_ws(), we should fit the peak center when doing const-q binning
   instead of taking the center of the peak range.