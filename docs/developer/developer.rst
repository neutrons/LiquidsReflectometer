.. _developer_documentation:

Developer Documentation
=======================

.. contents::
   :local:
   :depth: 1

Local Environment
-----------------
For purposes of development, create conda environment `lr_reduction` with file `environment.yml`, and then
install the package in development mode with `pip`:

.. code-block:: bash

   $> cd /path/to/lr_reduction/
   $> conda create env --solver libmamba --file ./environment.yml
   $> conda activate lr_reduction
   (lr_reduction)$> pip install -e ./

By installing the package in development mode, one doesn't need to re-install package `lr_reduction` in conda
environment `lr_reduction` after every change to the source code.

pre-commit Hooks
----------------

Activate the hooks by typing in the terminal:

.. code-block:: bash

   $> cd /path/to/lr_reduction/
   $> conda activate lr_reduction
   (lr_reduction)$> pre-commit install

Development procedure
---------------------

1. A developer is assigned with a task during neutron status meeting and changes the task's status to **In Progress**.
2. The developer creates a branch off *next* and completes the task in this branch.
3. The developer creates a pull request (PR) off *next*.
4. Any new features or bugfixes must be covered by new and/or refactored automated tests.
5. The developer asks for another developer as a reviewer to review the PR.
   A PR can only be approved and merged by the reviewer.
6. The developer changes the task’s status to **Complete** and closes the associated issue.

Updating mantid dependency
--------------------------
The mantid version and the mantid conda channel (`mantid/label/main` or `mantid/label/nightly`) **must** be
synchronized across these files:

- environment.yml
- conda.recipe/meta.yml
- .github/workflows/package.yml

Using the Data Repository liquidsreflectometer-data
---------------------------------------------------
To run the integration tests in your local environment, it is necessary first to download the data files.
Because of their size, the files are stored in the Git LFS repository
`lr_reduction-data <https://code.ornl.gov/sns-hfir-scse/infrastructure/test-data/liquidsreflectometer-data>`_.

It is necessary to have package `git-lfs` installed in your machine.

.. code-block:: bash

   $> sudo apt install git-lfs

After this step, initialize or update the data repository:

.. code-block:: bash

   $> cd /path/to/lr_reduction
   $> git submodule update --init

This will either clone `liquidsreflectometer-data` into `/path/to/lr_reduction/tests/liquidsreflectometer-data` or
bring the `liquidsreflectometer-data`'s refspec in sync with the refspec listed within file
`/path/to/liquidsreflectometer/.gitmodules`.

An intro to Git LFS in the context of the Neutron Data Project is found in the
`Confluence pages <https://ornl-neutrons.atlassian.net/wiki/spaces/NDPD/pages/19103745/Using+git-lfs+for+test+data>`_
(login required).


Coverage reports
----------------

GitHuh actions create reports for unit and integration tests, then combine into one report and upload it to
`Codecov <https://app.codecov.io/gh/neutrons/lr_reduction>`_.


Building the documentation
--------------------------
A repository webhook is setup to automatically trigger the latest documentation build by GitHub actions.
To manually build the documentation:

.. code-block:: bash

   $> conda activate lr_reduction
   (lr_reduction)$> make docs

After this, point your browser to
`file:///path/to/lr_reduction/docs/build/html/index.html`


Creating a stable release
-------------------------

- *patch* release, it may be allowed to bypass the creation of a candidate release.
  Still, we must update branch `qa` first, then create the release tag in branch `main`.
  For instance, to create patch version "v2.1.1":

.. code-block:: bash

   VERSION="v2.1.2"
   # update the local repository
   git fetch --all --prune
   git fetch --prune --prune-tags origin
   # update branch qa from next, possibly bringing work done in qa missing in next
   git switch next
   git rebase -v origin/next
   git merge --no-edit origin/qa  # commit message is automatically generated
   git push origin next  # required to "link" qa to next, for future fast-forward
   git switch qa
   git rebase -v origin/qa
   git merge --ff-only origin/next
   # update branch main from qa
   git merge --no-edit origin/main  # commit message is automatically generated
   git push origin qa  # required to "link" main to qa, for future fast-forward
   git switch main
   git rebase -v origin/main
   git merge --ff-only origin/qa
   git tag $VERSION
   git push origin --tags main

- *minor* or *major* release, we create a stable release *after* we have created a Candidate release.
  For this customary procedure, follow:

  + the `Software Maturity Model <https://ornl-neutrons.atlassian.net/wiki/spaces/NDPD/pages/23363585/Software+Maturity+Model>`_ for continous versioning as well as creating release candidates and stable releases.
  + Update the :ref:`Release Notes <release_notes>` with major fixes, updates and additions since last stable release.
