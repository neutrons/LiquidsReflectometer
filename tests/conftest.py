# standard imports
import os
from pathlib import Path

# third-party imports
import pytest
from mantid import simpleapi as mtd_api


@pytest.fixture(scope="session", autouse=True)
def mantid_facility():
    """Set the facility and instrument for the whole session.

    Mantid needs these to turn a run number into a file. `amend_config` only
    sets them when it is passed a `new_config`, and several modules pass just
    `data_dir` — so those modules were green only when collected after a module
    that happened to set them. Hoisting it here means no module can be silently
    coupled to another's collection order again.
    """
    mtd_api.config["default.facility"] = "SNS"
    mtd_api.config["default.instrument"] = "REF_L"


@pytest.fixture(autouse=True)
def _no_cwd_leak():
    """Fail the test that leaks a working directory, not one months later.

    A bare `os.chdir` in a test silently rebases every relative path in every
    test after it, which is how this suite's cwd-dependent failures stayed
    hidden and how a whole-suite green became unfalsifiable.
    """
    before = os.getcwd()
    yield
    assert os.getcwd() == before, (
        f"test changed the working directory and did not restore it: {before} -> {os.getcwd()}; "
        f"use monkeypatch.chdir, which pytest restores"
    )


@pytest.fixture(scope="session")
def datarepo_dir() -> str:
    r"""Absolute path to the event nexus files"""
    return str(Path(__file__).parent.parent / "tests/data/liquidsreflectometer-data")


@pytest.fixture(scope="session")
def nexus_dir() -> str:
    r"""Absolute path to the event nexus files"""
    return str(Path(__file__).parent.parent / "tests/data/liquidsreflectometer-data/nexus")


@pytest.fixture(scope="session")
def template_dir() -> str:
    r"""Absolute path to reduction/data/ directory"""
    return str(Path(__file__).parent / "data")
