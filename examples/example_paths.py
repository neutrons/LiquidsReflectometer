"""
Where the examples read their inputs and write their results.

No example names a folder of its own. Every input is a path below one data folder, laid out as the instrument's
folder is (IPTS-<n>/nexus, IPTS-<n>/shared/..., shared/...), and every result goes to one output folder:

    LR_REDUCTION_EXAMPLE_DATA      the data folder. Unset: the example data of this repository's test-data
                                   submodule, tests/data/liquidsreflectometer-data/examples.
    LR_REDUCTION_EXAMPLE_OUTPUT    the output folder. Unset: lr_reduction_example_output in the folder the
                                   example is started from. Created when missing.

At the facility, LR_REDUCTION_EXAMPLE_DATA=/SNS/REF_L reads the experiments' own files.
"""

import os
import sys
from pathlib import Path

DATA_VARIABLE = "LR_REDUCTION_EXAMPLE_DATA"
OUTPUT_VARIABLE = "LR_REDUCTION_EXAMPLE_OUTPUT"
REPOSITORY_DATA = Path(__file__).resolve().parent.parent / "tests" / "data" / "liquidsreflectometer-data" / "examples"
DEFAULT_OUTPUT_NAME = "lr_reduction_example_output"


def data_folder() -> Path:
    """The folder the examples read from."""
    chosen = os.environ.get(DATA_VARIABLE)
    folder = Path(chosen) if chosen else REPOSITORY_DATA
    if not folder.is_dir():
        raise FileNotFoundError(
            f"The example data folder {folder} does not exist. Set {DATA_VARIABLE} to a folder laid out as the "
            f"instrument's is (/SNS/REF_L at the facility), or fetch this repository's test data: "
            f"git submodule update --init && git submodule foreach git lfs pull"
        )
    return folder


def output_folder() -> Path:
    """The folder the examples write to, created when missing."""
    chosen = os.environ.get(OUTPUT_VARIABLE)
    folder = Path(chosen) if chosen else Path.cwd() / DEFAULT_OUTPUT_NAME
    folder.mkdir(parents=True, exist_ok=True)
    return folder


def run_examples(examples, arguments=None):
    """
    Run the examples named on the command line, or all of them, in the order given:

        python example_nr_reduction.py                              every example of the file
        python example_nr_reduction.py example_constant_tof_reduction

    examples: the example functions of the file
    arguments: the names to run (the command line's, when None)
    """
    by_name = {example.__name__: example for example in examples}
    names = sys.argv[1:] if arguments is None else list(arguments)
    unknown = [name for name in names if name not in by_name]
    if unknown:
        raise SystemExit(f"No example named {', '.join(unknown)}. This file has: {', '.join(by_name)}")
    return {name: by_name[name]() for name in (names or by_name)}
