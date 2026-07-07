import csv
from pathlib import Path

import numpy as np

cd_attenuation_file = Path(__file__).parent / "Cd-attenuation.csv"


def cd_attenuation() -> tuple[np.ndarray, np.ndarray]:
    """Reads the Cadmium attenuation data from a CSV file."""
    data = np.loadtxt(cd_attenuation_file, delimiter="\t", skiprows=1)
    wl, mu = data[:, 0], data[:, 1]
    wl_final = (2 * wl[-1]) - wl[-2]
    wl = np.append(wl, wl_final)
    return wl, mu
