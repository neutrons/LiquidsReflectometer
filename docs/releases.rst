.. _release_notes:

Release Notes
=============

Notes for major or minor releases. Notes for patch releases are deferred.

..
    Release notes are written in reverse chronological order, with the most recent release at the top,
    using the following format:

    ## <Next Major or Minor Release>
    (date of release, format YYYY-MM-DD)

    **Of interest to the User**:

    - PR #XYZ one-liner description

    **Of interest to the Developer:**

    - PR #XYZ one-liner description

2.1.13
------
2025-05-01

**Of interest to the User**:

- PR #61 Add support for custom instrument geometry settings

**Of interest to the Developer:**

- PR #66 Update mantid to 6.12.0

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
