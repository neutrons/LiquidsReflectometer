import logging
import os
import re
import sys
from pathlib import Path

repo_dir: str = str(Path(__file__).resolve().parents[2])


def extract_specifier(content: str, package: str) -> str:
    """
    Extract the full version specifier for a given package (e.g. 'mantid >=6.12.0').
    Handles spacing around operators.
    """
    pattern = rf"{package}\s*([<>=!~]+)\s*(\d+(?:\.\d+)*)"
    match = re.search(pattern, content)
    if match:
        return match.group(1) + match.group(2)
    else:
        raise ValueError(f"{package} specifier not found")


def check_dependencies_synced():
    """
    Check that the dependencies of environment.yml, meta.yaml and pyproject.toml are in sync
    """
    conda_env = Path(repo_dir, "environment.yml").read_text()
    pyproject_toml = Path(repo_dir, "pyproject.toml").read_text()
    conda_recipe = Path(repo_dir, "conda.recipe", "meta.yaml").read_text()

    # Extract full version specifiers (e.g., '>=6.12.0', '==6.12.0')
    mantid_env_spec = extract_specifier(conda_env, "mantid")
    mantid_pyproj_spec = extract_specifier(pyproject_toml, "mantid")
    mantid_recipe_spec = extract_specifier(conda_recipe, "mantid")

    if mantid_env_spec != mantid_pyproj_spec:
        raise RuntimeError("environment.yml and pyproject.toml ask different versions of mantid")
    if mantid_env_spec != mantid_recipe_spec:
        raise RuntimeError("environment.yml and meta.yaml ask different versions of mantid")


if __name__ == "__main__":
    try:
        check_dependencies_synced()
    except RuntimeError as e:
        logging.error(f"{e}")
        sys.exit(1)
    sys.exit(0)
