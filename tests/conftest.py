# standard imports
from pathlib import Path

# third-party imports
import pytest


def pytest_collection_modifyitems(items):
    r"""Mark every test that uses the data repository as ``datarepo``"""
    for item in items:
        if {"datarepo_dir", "nexus_dir"} & set(item.fixturenames):
            item.add_marker(pytest.mark.datarepo)


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
