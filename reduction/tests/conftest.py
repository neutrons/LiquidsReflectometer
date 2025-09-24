# standard imports
from pathlib import Path

# third-party imports
import pytest


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
    return str(Path(__file__).parent.parent / "data")
