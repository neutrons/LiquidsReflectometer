"""
Every example runs, start to end, on the example data.

The examples are the published description of how the reduction is called. A change that breaks one of them breaks
what the instrument's users do, so each example function is a test: it is called as written and must return and
write its results. An example whose input files are not on this machine is skipped; the summary at the end of the run names the example and the first
missing file, and the header of the run names the folder that was searched.
"""

import importlib
import importlib.util
import inspect
import os
from pathlib import Path

import pytest

import example_paths
from example_inputs import INPUTS

INSTRUMENT_FOLDER = "/SNS/REF_L/"  # the instrument's folder at the facility, where the reduction's defaults point

CASES = [(module, function) for module, functions in INPUTS.items() for function in functions]

# An example that warns needs looking at, so a warning fails its test. One kind is not the examples' and not the
# reduction's: drawing text, matplotlib calls pyparsing by names that pyparsing has since renamed.
WARNINGS = ["error"]
if importlib.util.find_spec("pyparsing.warnings"):
    WARNINGS.append("ignore::pyparsing.warnings.PyparsingDeprecationWarning")


def data_folder():
    return Path(os.environ.get(example_paths.DATA_VARIABLE) or example_paths.REPOSITORY_DATA)


def skip_reason(module, function, missing, total):
    return f"{module}.{function}: {len(missing)} of {total} input files not in the example data, first {missing[0]}"


@pytest.mark.parametrize("module", sorted(INPUTS))
def test_every_example_of_a_file_has_its_inputs_listed(module):
    """A new example gets a test only by being listed: an example that is not listed fails here."""
    written = {
        name
        for name, value in vars(importlib.import_module(module)).items()
        if name.startswith("example_") and inspect.isfunction(value) and value.__module__ == module
    }
    assert written == set(INPUTS[module])


def test_the_skip_reason_names_the_example_and_fits_a_report_line():
    """
    The reason is read in the summary at the end of the run (-rs), where pytest prints it whole but groups equal
    reasons under one line without the test's name. So it carries the example's name, and no absolute path: the data
    folder is in the header of the run.
    """
    reason = skip_reason("example_direct_beam", "example_direct_beam", ["IPTS-1/nexus/REF_L_1.nxs.h5", "b"], 4)

    assert reason.startswith("example_direct_beam.example_direct_beam: ")
    assert "2 of 4" in reason
    assert "IPTS-1/nexus/REF_L_1.nxs.h5" in reason
    assert str(data_folder()) not in reason


def read_outside_the_data_folder(opened, folder):
    """The files of the instrument's folder among those opened that are not below the example data folder."""
    return sorted({path for path in opened if path.startswith(INSTRUMENT_FOLDER) and not Path(path).is_relative_to(folder)})


def test_only_the_instruments_folder_counts_as_read_at_the_facility():
    """
    Home folders at the facility are below /SNS too: a checkout there, its settings.json and the Python environment
    are not the instrument's files. With the data folder at the instrument's folder nothing is outside it.
    """
    opened = [
        "/SNS/users/someone/checkout/src/lr_reduction/settings.json",
        "/SNS/users/someone/.cache/matplotlib/fontlist.json",
        "/SNS/REF_L/shared/Cd_transmission/Cd_mu_2025.dat",
        "/tmp/data/shared/Cd_transmission/Cd_mu_2025.dat",
    ]

    assert read_outside_the_data_folder(opened, Path("/tmp/data")) == ["/SNS/REF_L/shared/Cd_transmission/Cd_mu_2025.dat"]
    assert read_outside_the_data_folder(opened, Path("/SNS/REF_L")) == []


def test_every_example_file_is_listed():
    files = {path.stem for path in Path(example_paths.__file__).parent.glob("example_*.py")} - {"example_paths"}
    assert files == set(INPUTS)


@pytest.mark.datarepo
@pytest.mark.filterwarnings(*WARNINGS)
@pytest.mark.parametrize(("module", "function"), CASES, ids=[f"{module}.{function}" for module, function in CASES])
def test_example_runs(module, function, example_output, files_opened):
    missing = [path for path in INPUTS[module][function] if not (data_folder() / path).is_file()]
    if missing:
        pytest.skip(skip_reason(module, function, missing, len(INPUTS[module][function])))

    result = getattr(importlib.import_module(module), function)()

    assert result is not None
    written = [path for path in example_output.rglob("*") if path.is_file()]
    assert written, "the example wrote nothing to the output folder"
    elsewhere = read_outside_the_data_folder(files_opened, data_folder())
    assert not elsewhere, f"read in the instrument's folder and not below the example data folder, so CI does not have it: {elsewhere}"
