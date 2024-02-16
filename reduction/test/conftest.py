# standard imports
import logging
from pathlib import Path
import os


def pytest_sessionstart(session):
    r"""invoked by pytest at the very beginning"""
    logger = logging.getLogger()
    # Insert data directory in Mantid config file
    mantid_config_dir = Path(os.environ["HOME"]) / ".mantid"
    os.makedirs(str(mantid_config_dir), exist_ok=True)  # create directory if it doesn't exists
    mantid_config_file = mantid_config_dir / "Mantid.user.properties"
    write_mode = "a" if mantid_config_file.exists() else "w"  # append or create-then-write
    with open(mantid_config_file, write_mode) as file_handle:
        data_path = str(Path(__file__).parent.parent / "tests/data/liquidsreflectometer-data/nexus")
        file_handle.write(f"datasearch.directories={data_path}")
        logger.info(f"Appending data directory {data_path} to mantid config file {str(mantid_config_file)}")
