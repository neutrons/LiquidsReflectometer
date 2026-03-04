"""Background worker thread for template-based reduction."""

import sys
import traceback

from qtpy.QtCore import QThread, Signal


class ReductionWorker(QThread):
    """Run reduce_from_template in a background thread with progress reporting."""

    progress = Signal(int, int, str)  # overall_pct, total, description
    finished = Signal(list)           # list of result dicts
    error = Signal(str)               # error message

    def __init__(self, run_numbers, template_file, experiment_id, datapath=None,
                 template_path=None, override_params=None, plot=False,
                 save_plots=False, plot_dir=None, log_to_stdout=False,
                 parent=None):
        super().__init__(parent)
        self.run_numbers = run_numbers
        self.template_file = template_file
        self.experiment_id = experiment_id
        self.datapath = datapath
        self.template_path = template_path
        self.override_params = override_params
        self.plot = plot
        self.save_plots = save_plots
        self.plot_dir = plot_dir
        self.log_to_stdout = log_to_stdout
        self._cancelled = False

    def cancel(self):
        self._cancelled = True

    def run(self):
        import matplotlib
        matplotlib.use('Agg')

        try:
            from lr_reduction.new_reduction_from_template import reduce_from_template

            total_runs = len(self.run_numbers)
            results = []

            for run_idx, runno in enumerate(self.run_numbers):
                if self._cancelled:
                    self.progress.emit(100, 100, "Cancelled")
                    break

                if self.log_to_stdout:
                    print(f"[Worker] Starting run {runno} ({run_idx + 1}/{total_runs})")

                def _progress_cb(phase, total, desc):
                    overall = int(((run_idx * 5 + phase) / (total_runs * 5)) * 100)
                    if self.log_to_stdout:
                        print(f"[Worker] Run {runno}: {desc} ({overall}%)")
                    self.progress.emit(overall, 100, f"Run {runno}: {desc}")

                result = reduce_from_template(
                    runno=runno,
                    template_file=self.template_file,
                    experiment_id=self.experiment_id,
                    datapath=self.datapath,
                    template_path=self.template_path,
                    override_params=self.override_params,
                    plot=False,  # never show interactive plots from thread
                    progress_callback=_progress_cb,
                    save_plots=self.save_plots,
                    plot_dir=self.plot_dir,
                )
                results.append(result)

            if not self._cancelled:
                self.progress.emit(100, 100, "Complete")
                if self.log_to_stdout:
                    print(f"[Worker] All {total_runs} run(s) completed successfully")
            else:
                if self.log_to_stdout:
                    print(f"[Worker] Cancelled after {len(results)}/{total_runs} run(s)")
            self.finished.emit(results)

        except Exception:
            tb = traceback.format_exc()
            if self.log_to_stdout:
                print(f"[Worker] ERROR:\n{tb}", file=sys.stderr)
            self.error.emit(tb)
