.. _developer_documentation:

Developer Documentation
=======================

.. contents::
   :local:
   :depth: 1


Local Environment
-----------------
We use `pixi <https://pixi.sh/latest/>`_ to create the conda environment for local development of the project,
as well as for creating conda packages for this project.

If you don't have `pixi` installed in your Linux machine, you can install it with:

.. code-block:: bash

   $> curl -fsSL https://pixi.sh/install.sh | sh

See the `pixi installation page <https://pixi.sh/latest/installation/>`_ for more options.

Then, you can create the conda environment for local development with:

.. code-block:: bash

   $> cd /path/to/lr_reduction/
   $> pixi install

This command will also install project `lr_reduction` in editable mode.
By installing the project in editable mode,
one doesn't need to re-install the project after a change in the source code.

To activate the conda environment, type in the terminal:

.. code-block:: bash

   $> pixi shell

This command will start a bash shell with the conda environment activated. To exit the shell, type:

.. code-block:: bash

   $> exit

Pixi also offers a set of commands to help with the development of the project,
like building the documentation, and creating conda packages.
To see the list of available commands, type in the terminal:

.. code-block:: bash

   $> pixi run
   Available tasks:
      build-conda
      build-docs
      clean-all
      clean-conda
      clean-docs
      conda-builder
      reset-version
      sync-version

Each task has a brief description in file pyproject.toml, under the section `[tool.pixi.tasks]`.


Activating the Environment Automatically
----------------------------------------
Wouldn't be nice if every time you enter the project directory, the conda environment is activated automatically?
To achieve this, install `direnv <https://direnv.net/docs/installation.html>`_
and create file `.envrc` in the project root directory with the following content:

.. code-block:: bash

   watch_file pixi.lock
   eval "$(pixi shell-hook)"
   unset PS1

Then, in the terminal, type:

.. code-block:: bash

   $> direnv allow

Now direnv activates the environment when you enter the project directory,
and deactivates it when you leave the directory.

Line `watch_file pixi.lock` directs direnv to re-evaluate the environment whenever file `pixi.lock` changes.
Line `unset PS1` prevents direnv from
`reporting on a nagging, albeit harmless, error message <https://github.com/direnv/direnv/wiki/PS1>`_.


pre-commit Hooks
----------------

Activate the hooks by typing in the terminal:

.. code-block:: bash

   $> cd /path/to/lr_reduction/
   $> pixi shell
   $> pre-commit install


Development procedure
---------------------

1. A developer is assigned with a task during neutron status meeting and changes the task's status to **In Progress**.
2. The developer creates a branch off *next* and completes the task in this branch.
3. The developer creates a pull request (PR) off *next*.
4. Any new features or bugfixes must be covered by new and/or refactored automated tests.
5. The developer asks for another developer as a reviewer to review the PR.
   A PR can only be approved and merged by the reviewer.
6. The developer changes the task’s status to **Complete** and closes the associated issue.


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

GitHub actions create reports for unit and integration tests, then combine into one report and upload it to
`Codecov <https://app.codecov.io/gh/neutrons/lr_reduction>`_.


Building the documentation
--------------------------
A repository webhook is setup to automatically trigger the latest documentation build by GitHub actions.
To manually build the documentation:

.. code-block:: bash

   $> pixi run build-docs

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
