# Plan: Add "Template Reduction" Tab to lr_reduction Launcher

## Context

The `bvacaliuc/new_workflow_ui_plan` branch added `reduce_from_template()` — a new entry point for running reflectometry reduction from an XML template file. An example script (`example_nr_from_template.py`) demonstrates its usage but there is no GUI integration. The goal is to add a new tab to the existing Qt launcher that exposes all configuration options, runs reduction in a background thread with progress feedback, and supports saving/controlling plot output.

## Files to Modify

| File | Action | Purpose |
|------|--------|---------|
| `launcher/apps/template_reduce.py` | **Create** | New tab widget with all UI controls |
| `launcher/apps/reduction_worker.py` | **Create** | QThread worker for non-blocking reduction |
| `launcher/launcher.py` | **Modify** | Register the new tab (3 lines) |
| `src/lr_reduction/new_reduction_from_template.py` | **Modify** | Add progress_callback, save_plots, plot_dir params |
| `src/lr_reduction/nr_reduction_config.py` | **Modify** | Add `plot_save_dir` attribute |
| `src/lr_reduction/nr_reduction_calc.py` | **Modify** | Replace bare `plt.show()` with save-or-show helper |
| `tests/unit/lr_reduction/test_template_reduce_parsing.py` | **Create** | Test run number parsing logic |
| `tests/unit/lr_reduction/test_plot_save.py` | **Create** | Test plot save/show behavior |

## Implementation Steps

### Step 1: Backend — `nr_reduction_config.py`

Add one attribute after line 48 (`self.plotQ4 = False`):

```python
self.plot_save_dir = None  # Directory to save diagnostic plots (None = don't save)
```

### Step 2: Backend — `new_reduction_from_template.py`

**2a. Extend `reduce_from_template()` signature:**

```python
def reduce_from_template(runno, template_file, experiment_id, datapath=None,
                         template_path=None, override_params=None, plot=True,
                         progress_callback=None, save_plots=False, plot_dir=None):
```

All new params default to no-op values — existing callers unaffected.

**2b. Insert progress callbacks** between the 5 phases (read NEXUS, read template, reduce, assemble, save):

```python
def _report(phase, total, desc):
    if progress_callback:
        progress_callback(phase, total, desc)
```

Call `_report(0, 5, "Reading NEXUS metadata")` before h5py open, `_report(1, 5, "Configuring reduction")` before template read, etc.

**2c. Modify `plot_reflectivity()`** — add `save_path=None` and `show=True` params:

```python
def plot_reflectivity(data_array, RQ4=False, log_x=True, save_path=None, show=True):
    fig, ax = plt.subplots()
    # ... existing errorbar code ...
    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close(fig)
```

**2d. Modify `assemble_results()`** — add `save_plots=False, plot_dir=None, show_plots=True` params. Pass to `plot_reflectivity()`. Save as `REFL_{seq_id}_individual_settings.png`.

**2e. Update calls** in `reduce_from_template()` body:
- Pass `save_plots`, `plot_dir`, `show_plots=plot` through to `assemble_results()`
- Final reflectivity plot: save as `REFL_{seq_id}_combined.png` if `save_plots=True`

### Step 3: Backend — `nr_reduction_calc.py`

**3a. Add helper method to `NR_Reduction`:**

```python
def _show_or_save_plot(self, fig, name_hint):
    if getattr(self.config, 'plot_save_dir', None) is not None:
        save_path = Path(self.config.plot_save_dir) / f"{name_hint}.png"
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
    if self.config.plotON:
        plt.show()
    else:
        plt.close(fig)
```

**3b. Replace 5 `plt.show()` calls** (lines ~142, 420, 569, 644, 721) with `self._show_or_save_plot(fig, "descriptive_name")`. Ensure each plot block assigns to `fig` (most already do; the `reduce()` method at line ~128 uses bare `plt.errorbar` and needs wrapping in `fig, ax = plt.subplots()`).

### Step 4: Worker Thread — `launcher/apps/reduction_worker.py`

```python
class ReductionWorker(QThread):
    progress = Signal(int, int, str)  # overall_pct, total, description
    finished = Signal(list)           # list of result dicts
    error = Signal(str)               # error message
```

- Constructor receives all reduction parameters
- `run()` method:
  1. Set `matplotlib.use('Agg')` for thread safety
  2. Import `reduce_from_template`
  3. Loop over run numbers, calling `reduce_from_template()` with progress callback
  4. Progress callback maps per-run phases to overall percentage: `((run_idx * 5 + phase) / (total_runs * 5)) * 100`
  5. Emit `finished` signal with all results
  6. Catch exceptions, emit `error` signal
- `cancel()` method sets `_cancelled` flag checked between runs

### Step 5: Tab Widget — `launcher/apps/template_reduce.py`

**Layout (QVBoxLayout top-level):**

**Top section (fixed, QGridLayout):**

| Row | Col 1 | Col 2 | Col 3 |
|-----|-------|-------|-------|
| 0 | Run number(s) QLineEdit | | "Run number(s) — comma-separated or range" label |
| 1 | Template file QPushButton | | QLabel path (*.xml file browser) |
| 2 | Experiment ID QLineEdit | | "e.g. IPTS-36119" label |
| 3 | Output directory QPushButton | | QLabel path (folder browser) |
| 4 | Data path QPushButton | | QLabel path (folder browser, NEXUS dir) |
| 5 | DB path QPushButton | | QLabel path (folder browser) |
| 6 | DB file name QPushButton | | QLabel path (*.dat file browser) |
| 7 | Template save path QPushButton | | QLabel path (folder browser) |

**Middle section (QScrollArea with QGroupBox groups):**

- **Processing Flags**: Method (QComboBox: meanTheta/constantQ/constantTOF), Normalize (QCheckBox), AutoScale (QCheckBox), useCalcTheta (QCheckBox)
- **Plot Options**: Save plots (QCheckBox), Interact with plots (QCheckBox), Plot as RQ4 (QCheckBox)
- **Q-space Parameters**: qmin, qmax, dqbin, Qline_threshold (QLineEdit + QDoubleValidator)
- **Dead Time**: dead_time, dead_time_tof_step (QLineEdit + QDoubleValidator)
- **Detector/Peak**: DetResFn (QComboBox), DetSigma, peak_pad, peak_type (QComboBox)
- **Emission Time**: use_emission_time (QCheckBox)
- **Instrument Geometry**: IncidentTheta, plus optional fields for mmpix, dSampDet, ny, nx, dMod, xi_ref, dS1Samp (blank = use defaults from settings.json)

**Bottom section (fixed):**

| Widget | Description |
|--------|-------------|
| QPushButton("REDUCE") | Green background, triggers reduction |
| QProgressBar | 0-100, shows overall progress |
| QLabel | Status text ("Run 211029: Running reduction...") |

**Key methods following existing patterns:**
- `read_settings()` / `save_settings()` — QSettings with `"tmpl_"` prefix
- `check_inputs()` — validate run numbers parseable, template file exists, output dir exists, experiment ID non-empty
- `show_dialog()` — error message box (same as other tabs)
- `_parse_run_numbers(text)` — parse "211029,211030" or "211029-211031" into list[int]
- `_build_override_params()` — collect non-default widget values into dict
- `reduce()` — validate, save settings, disable button, start worker, connect signals
- `on_progress()` — update progress bar and status label
- `on_finished()` — re-enable button, optionally show interactive plot from main thread using returned results, show completion dialog
- `on_error()` — re-enable button, show error dialog

**Interactive plot handling:** Worker always uses Agg backend (no plt.show()). After worker finishes, if "Interact" is checked, the main thread calls `plot_reflectivity()` with the returned results dict to show interactive matplotlib windows.

### Step 6: Wire Up — `launcher/launcher.py`

Add import:
```python
from apps.template_reduce import TemplateReduce
```

Add tab after SLD calculator:
```python
tab_id += 1
self.template_tab = TemplateReduce()
self.addTab(self.template_tab, "Template reduction")
self.setTabText(tab_id, "Template reduction")
```

### Step 7: Tests (Red-Green TDD)

**`tests/unit/lr_reduction/test_template_reduce_parsing.py`:**
- Test `_parse_run_numbers()` with single, CSV, range, mixed, and invalid inputs

**`tests/unit/lr_reduction/test_plot_save.py`:**
- Test `plot_reflectivity()` saves file when `save_path` is given
- Test `plot_reflectivity()` doesn't call `plt.show()` when `show=False`
- Test `NRReductionConfig` has `plot_save_dir` attribute defaulting to `None`
- Test `reduce_from_template()` accepts `progress_callback` parameter (mock h5py/file ops)

## Verification

1. Run `pixi run test-reduction` to ensure existing tests still pass
2. Launch the UI with `python launcher/launcher.py` from the pixi environment
3. Verify the new "Template reduction" tab appears
4. Test file/folder browsers open with correct filters
5. Test settings persistence (close and reopen)
6. Run a reduction with test data:
   - Template: `/SNS/REF_L/shared/lr_reduction/new_workflow_test_outputs/test_template.xml`
   - Run: 211029
   - Experiment: IPTS-36119
   - Data path: `/SNS/REF_L/IPTS-30101/nexus`
   - DB path: `/SNS/REF_L/shared/Cd_DB_processing/DBs/`
   - DB file: `A1_air_div1_Cd_DB.dat`
7. Verify progress bar updates during reduction
8. Verify plots are saved to output directory when "Save plots" is checked
9. Verify interactive plots appear only when "Interact" is checked
10. Test cancellation during multi-run reduction
