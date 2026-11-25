"""
reduce_from_template finds the file it has just written.

It saves one run as a "partial" file and then collects every partial file of the sequence from the output folder.
The name it writes and the name it looks for are made in two places, and must agree.
"""

import numpy as np

from lr_reduction import new_reduction_from_template as template
from lr_reduction.nr_reduction_config import NRReductionConfig


def config_of_run(folder, subname=""):
    config = NRReductionConfig()
    config.Spath = folder
    config.Sname = "REFL_211029_1_211029"
    config.subname = subname  # NR_Reduction leaves it empty, or starting with an underscore
    return config


def test_the_partial_file_of_a_run_is_found_when_the_sequence_is_assembled(tmp_path):
    config = config_of_run(tmp_path)
    q = np.linspace(0.01, 0.05, 5)
    np.savetxt(tmp_path / f"{template.partial_name(config)}.dat", np.array([q, q * 2, q * 0.1, q * 0.01]).T)

    sequence_numbers, runs, combined, _, _ = template.assemble_results(211029, tmp_path, autoscale=False, plot=False)

    assert (sequence_numbers, runs) == ([1], [211029])
    np.testing.assert_allclose(combined["Q"], q)


def test_the_names_carry_the_subname_when_there_is_one(tmp_path):
    config = config_of_run(tmp_path, subname="_trial")

    assert template.partial_name(config) == "REFL_211029_1_211029_trial_partial"
    assert template.combined_name(211029, config) == "REFL_211029_trial_combined_data"


def test_the_combined_file_is_not_taken_for_a_partial_one(tmp_path):
    config = config_of_run(tmp_path)

    assert template.combined_name(211029, config) == "REFL_211029_combined_data"
