.. _release_notes:

Release Notes
=============

Notes for major or minor releases. Notes for patch releases are deferred.

..
    2.3.0
    -----
    (date of release, format YYYY-MM-DD)

    **Of interest to the User**:

    - PR #XYZ one-liner description

    **Of interest to the Developer:**

    - PR #73 merge test and deploy workflows
    - PR #70 Introduce pixi to the project
..


2.2.0
------
2025-07-22

**Of interest to the Developer:**

- PR #76 Adds new option, dead-time threshold - an upper limit for dead-time correction ratios
- PR #74 Update actions, documentation, and versioningit
- PR #67 Make the documentation work with readthedocs
- PR #66 Modernizes conda packaging


2.1.13
------
2025-05-01

**Of interest to the User**:

- PR #61 Add support for custom instrument geometry settings

**Of interest to the Developer:**

- PR #65 Update mantid to 6.12.0


2.1.0
-----

**Of interest to the User**:

- PR #33 enable dead time correction for runs with skipped pulses
- PR #26 add dead time correction to the computation of scaling factors
- PR #23 add dead time correction
- PR #19 Functionality to use two backgrounds
- PR #15 Ability to fit a background with a polynomial function

**Of interest to the Developer:**

- PR #40 documentation to create a patch release
- PR #37 documentation conforming to that of the python project template
- PR #36 versioning with versioningit
- PR #25 Read in error events when computing correction
- PR #21 switch dependency from mantidworkbench to mantid
- PR #20 allow runtime initialization of new attributes for ReductionParameters
- PR #14 add first GitHub actions
- PR #12 switch from mantid to mantidworkbench conda package

2.1.6
------
09/2024

- Implement emission time correction

2.0.26
------
05/2024

- Implement a better dead time correction using weighted events

2.0.25
------
04/2024

- Use dead time parameters from template

2.0.24
------
04/2024

- Fix issue with errors when using dead time correction

2.0.22
------
04/2024

- Add dead time correction to scaling factor calculation

2.0.20
------
03/2024

- Add error events to dead time correction

2.0.19
------
03/2024

- Add dead time correction and clean up functional background option

2.0.14
------
02/2024

- Add functional background

2.0.13
------
08/2023

- Get correct angle with free liquids

2.0.9
-----
04/2023

- Subtract normalization background & add x-direction option

2.0.8
-----
04/2023

- Move pre/post cuts to template.

2.0.7
-----
03/2023

- Update parameters that will be read from file by Mantid.
