## Event-based reduction for the Liquids Reflectometer

### template.py
`template.py` is used to load a template and launch the reductiuon for a single run.
`It will also apply the scaling factors.

`process_from_template()` will load the data, but `process_from_template_ws` can
also be called on a Mantid workspace. This makes it possible to reduce event-filtered data.
`reduce_30Hz_slices` will reduce time-resolved data by splitting the run in time intervals.

### event_reduction.py
`event_reduction.py` implements the data reduction for a single Mantid workspace.
It takes the scattering and direct beam workspaces, along with all the necessary reduction
parameters. It can do specular and off-specular. It can reduced specular reflectivity in
constant-Q binning or not, and can estimate the background in both pixel-TOF space or pixel-Q space.





### To-do list

 - Write something to create proper output with options to combine overlapping points
 - Write ORSO writer
 - Output xml/json with parameters used (like the autoreduction)
 - Write auto-reduction wrapper to leave partial results on disk and combine them