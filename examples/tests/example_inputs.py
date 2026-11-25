"""
The files each example reads, as paths below the example data folder (examples/example_paths.py).

The tests use this list to skip an example whose inputs are not there, and it is the list of what the example
data of the test-data repository must hold. Printed one path per line, it feeds a copy:

    python examples/tests/example_inputs.py | rsync -a --files-from=- <instrument folder>/ <example data folder>/
"""


def nexus(experiment, *runs):
    return [f"{experiment}/nexus/REF_L_{run}.nxs.h5" for run in runs]


def files(folder, *names):
    return [f"{folder}/{name}" for name in names]


_SHARED = "shared/lr_reduction/new_workflow_test_outputs"
_SETTINGS_JSON = files(_SHARED, "REFL_227158_settings.json")
_DB_36776 = files("IPTS-36776/shared/transmission", "DB_A1_Cd.txt", "DB_A2_Cd.txt", "DB_A3_Cd.txt")
_DB_36681 = files("IPTS-36681/shared/transmission", "DB_A1_Cd.txt", "DB_A2_Cd.txt", "DB_A3_Cd.txt")
_FROM_JSON = _SETTINGS_JSON + _DB_36776 + nexus("IPTS-36776", 227164, 227165, 227166)
_DIRECT_CONFIG_36681 = _DB_36681 + nexus("IPTS-36681", 228220, 228221, 228222)
_DIRECT_CONFIG_36776 = _DB_36776 + nexus("IPTS-36776", 227147, 227148, 227149)
_TEMPLATE_DB = files("shared/EBW/Cd_DB_processing/DBs", "A1_air_div1_Cd_DB.dat", "A2_air_div1_Cd_DB.dat", "A3_air_div1_Cd_DB.dat")
_CD_ATTENUATION = files("shared/Cd_transmission", "Cd_mu_2025.dat")
_TEMPLATE_RUNS = nexus("IPTS-30101", 211029, 211030, 211031)

INPUTS = {
    "example_direct_beam": {
        "example_direct_beam": _CD_ATTENUATION + nexus("IPTS-36119", 225318, 225319, 225320, 225321),
        "example_direct_beam_pixel_offset": _CD_ATTENUATION + nexus("IPTS-36119", 235850),
        "example_direct_beam_with_cd_list": nexus("IPTS-36119", 231187, 231189, 231191),
    },
    "example_nr_reduction": {
        "example_mean_theta_reduction": _DIRECT_CONFIG_36681,
        "example_constant_tof_reduction": _DIRECT_CONFIG_36681,
        "example_mean_theta_reduction_8col": _DIRECT_CONFIG_36776,
        "example_mean_theta_reduction_savepdf": _DIRECT_CONFIG_36776,
        "example_mean_theta_reduction_savepdf_noshow": _DIRECT_CONFIG_36776,
    },
    "example_nr_reduction_from_file": {
        "example_from_dat": files("IPTS-36776/shared/EBW_reduced", "NR_meanTheta_8test_combined_8col.dat")
        + _DB_36776
        + nexus("IPTS-36776", 227158, 227159, 227160, 227161, 227162),
        "example_from_json": _FROM_JSON,
        "example_wrong_order": _FROM_JSON,
        "example_one_missing": _FROM_JSON,
        "example_with_search": _FROM_JSON,
        "example_with_savepdf": _FROM_JSON,
        "example_with_savepdf_noshow": _FROM_JSON,
        "example_with_search_savepdf": _FROM_JSON,
        "example_with_8col": _FROM_JSON,
        "example_with_8col_search": _FROM_JSON,
        "example_for_autoreduction": files("IPTS-35331/shared/autoreduce", "reduce_settings_down.json")
        + files("IPTS-35331/shared/transmission", "A1_Si_Cd.txt", "A2_Si_Cd.txt", "A3_Si_Cd.txt")
        + nexus("IPTS-35331", 230395),
        "example_incidenttheta": files("IPTS-36119/shared/reduced/Aug2026", "REFL_231035_settings.json")
        + files("IPTS-36119/shared/transmission/Aug2026", *(f"A{angle}_div5_Cd.txt" for angle in range(1, 7)))
        + nexus("IPTS-36119", 231035),
    },
    "example_nr_from_template": {
        "example_template_reduction": files(_SHARED, "test_template.xml") + _TEMPLATE_DB[:1] + _TEMPLATE_RUNS[:1],
        "example_template_reduction_3ang": files(_SHARED, "test_template_3ang.xml") + _TEMPLATE_DB + _TEMPLATE_RUNS,
        "example_template_reduction_new": files(_SHARED, "REFL_211029_template_new.xml") + _TEMPLATE_DB + _TEMPLATE_RUNS,
    },
    "example_time_slices": {
        "example_slices_from_json_2": _FROM_JSON,
        "example_slices_from_json_5": _FROM_JSON,
        "example_times_from_json": _FROM_JSON,
    },
}


def all_inputs():
    return sorted({path for examples in INPUTS.values() for paths in examples.values() for path in paths})


if __name__ == "__main__":
    print("\n".join(all_inputs()))
