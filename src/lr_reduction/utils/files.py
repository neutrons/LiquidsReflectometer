from pathlib import Path


def get_sequence_id_from_path(path: str | Path) -> str:
    """Extracts the sequence ID from the configuration file path.

    Looks for the first underscore-delimited, all-digit token in the file stem, so both
    `template_201282.xml` and the standard SNS `REF_L_201282_auto_template.xml` naming
    convention (a multi-token instrument prefix before the run number) resolve correctly.
    """
    path = Path(path)
    for token in path.stem.split("_"):
        if token.isdigit():
            return token
    raise ValueError(f"Invalid configuration file name: {path.name}")
