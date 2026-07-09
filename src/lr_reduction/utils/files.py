from pathlib import Path


def get_sequence_id_from_path(path: str | Path) -> str:
    """Extracts the sequence ID from the configuration file path."""
    path = Path(path)
    try:
        return path.stem.split("_")[1]
    except IndexError:
        raise ValueError(f"Invalid configuration file name: {path.name}")
