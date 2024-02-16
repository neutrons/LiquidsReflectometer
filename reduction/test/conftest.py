# standard imports
import logging
from pathlib import Path
import os

data_dir = Path(__file__).parent.parent / "data"


def pytest_sessionstart(session):
    r"""invoked by pytest at the very beginning"""
    logger = logging.getLogger()
    # Insert data directory in Mantid config file
    mantid_config_dir = Path(os.environ["HOME"]) / ".mantid"
    os.makedirs(str(mantid_config_dir), exist_ok=True)  # create directory if it doesn't exists
    mantid_config_file = mantid_config_dir / "Mantid.user.properties"
    write_mode = "a" if mantid_config_file.exists() else "w"  # append or create-then-write
    with open(mantid_config_file, write_mode) as file_handle:
        data_path = data_dir / "liquidsreflectometer-data" / "nexus"
        file_handle.write(f"datasearch.directories={str(data_path)}")
        logger.info("Appending data directory to mantid config file")
