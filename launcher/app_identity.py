#!/usr/bin/python3
"""The launcher's QSettings identity, in one place.

Every settings layer — the per-tab UI preferences this slug persists, and the
reduction- and global-settings layers that follow — must resolve to the same
store. QSettings() derives its path from the QCoreApplication identity, so the
identity has to be established before the first QSettings is constructed. Both
shipped entry points (launcher, new_launcher) install it; tabs call
ensure_identity() themselves so a non-GUI entry point resolves the same store
rather than a nameless one.

Installing the identity **repoints the store for the whole process**. Before it,
every launcher tab read and wrote the identity-less
`Unknown Organization.conf`; after it, they resolve to
ORNL/lr_reduction_new_launcher. LEGACY_KEYS below is what that costs if it is
not carried across — measured from the source, not estimated.
"""

from qtpy import QtCore

ORG_NAME = "ORNL"
ORG_DOMAIN = "ornl.gov"
APP_NAME = "lr_reduction_new_launcher"

# Written into the new store once the one-shot copy has run, so the migration
# is gated on a recorded fact rather than inferred from whether data looks
# present.
MIGRATION_SENTINEL = "settings_migrated_from_legacy"

# Every key the launcher modules read or write, enumerated per module from
# `settings.value(...)` / `settings.setValue(...)` call sites under launcher/.
# Enumerated deliberately, never a wholesale allKeys() copy: the identity-less
# store is the shared dumping ground for every Qt application on the machine
# that never set an organization, and several of these names are generic enough
# (output_dir, db_output_dir) that a blind copy would import a stranger's value.
LEGACY_KEYS = (
    # dynamic_30Hz.py
    "30Hz_data_run_number",
    "30Hz_output_dir",
    "30Hz_reference",
    "30Hz_ref_run_number",
    "30Hz_template",
    "30Hz_time_slice",
    # dynamic_60Hz.py
    "60Hz_data_run_number",
    "60Hz_fix_offset",
    "60Hz_output_dir",
    "60Hz_scan_index",
    "60Hz_template",
    "60Hz_time_slice",
    "60Hz_use_fix_offset",
    # file_batch.py
    "settings_datapath",
    "settings_DBpath",
    "settings_dir",
    "settings_enable_browse",
    "settings_experiment_id",
    "settings_file",
    "settings_max_plots",
    "settings_plot",
    "settings_runs",
    "settings_save_summary",
    "settings_Spath",
    "settings_subname",
    # off_spec.py
    "offspec_output_dir",
    "offspec_run_number",
    "offspec_wl_step",
    # overplot.py
    "overplot_folder",
    "overplot_xscale",
    "overplot_ytransform",
    # quick_reduce.py
    "db_output_dir",
    "quick_db_pixel",
    "quick_db_run_number",
    "quick_output_dir",
    "quick_pixel",
    "quick_run_number",
    # reduction.py
    "fit_first_peak",
    "reduction_avg_overlap",
    "reduction_const_q",
    "reduction_first_run_number",
    "reduction_fix_offset",
    "reduction_last_run_number",
    "reduction_template",
    "reduction_use_fix_offset",
    "reduction_use_old",
    # refracted.py
    "refracted_material",
    "refracted_output_dir",
    "refracted_run_number",
    # sld_calculator.py
    "sld_composition",
    "sld_wavelength",
    # template_batch.py
    "template_datapath",
    "template_DBpath",
    "template_dir",
    "template_enable_browse",
    "template_experiment_id",
    "template_file",
    "template_max_plots",
    "template_plot",
    "template_runs",
    "template_save_summary",
    "template_Spath",
    "template_subname",
    # xrr.py
    "output_dir",
    "xrr_data_file",
)


def ensure_identity():
    """Install the launcher identity if it is not already in place.

    Guarded on the organization name alone. applicationName() is auto-derived
    from argv[0] by QApplication, so a script constructed as
    QApplication(sys.argv) already has one — an `or applicationName()` guard
    would early-return there and leave the process resolving to
    `Unknown Organization/<script>.conf`, which is exactly the split this
    module exists to prevent (and T2/T3 are contractually pointed here).

    Non-clobbering by design: a test fixture's throwaway organization survives,
    which is what keeps tab construction safe under an isolated root.
    """
    app = QtCore.QCoreApplication
    if app.organizationName():
        return
    app.setOrganizationName(ORG_NAME)
    app.setOrganizationDomain(ORG_DOMAIN)
    app.setApplicationName(APP_NAME)


def migrate_legacy_settings():
    """One-shot copy of the pre-identity keys into the new store.

    Installing the identity repoints the store for the entire process, so
    without this every launcher tab would appear to forget its settings once —
    66 keys across 11 modules, in a slug whose whole symptom statement is "the
    GUI forgets typed values".

    Gated on MIGRATION_SENTINEL, written after a successful copy, and it never
    overwrites a value already present in the new store. Not re-entrant: it
    blanks the identity statics to read the nameless store, so it must run once
    on the main thread before any other thread constructs QSettings.
    """
    current = QtCore.QSettings()
    if current.value(MIGRATION_SENTINEL, False) in (True, "true", "True"):
        return

    app = QtCore.QCoreApplication
    org, domain, name = app.organizationName(), app.organizationDomain(), app.applicationName()
    # The try opens BEFORE the clears: a raise between them would otherwise
    # leave the process nameless for the rest of its life.
    try:
        app.setOrganizationName("")
        app.setOrganizationDomain("")
        app.setApplicationName("")
        legacy = QtCore.QSettings()
        carried = {key: legacy.value(key) for key in LEGACY_KEYS if legacy.contains(key)}
    finally:
        app.setOrganizationName(org)
        app.setOrganizationDomain(domain)
        app.setApplicationName(name)

    for key, value in carried.items():
        if not current.contains(key):
            current.setValue(key, value)
    current.setValue(MIGRATION_SENTINEL, True)
    current.sync()
