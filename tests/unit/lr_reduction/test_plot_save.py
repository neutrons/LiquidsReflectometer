"""Tests for plot save/show behavior in the new reduction workflow."""

import os
import tempfile
from unittest.mock import patch, MagicMock

import numpy as np
import pytest


class TestNRReductionConfigPlotSaveDir:
    """Test that NRReductionConfig has plot_save_dir attribute."""

    def test_plot_save_dir_defaults_to_none(self):
        from lr_reduction.nr_reduction_config import NRReductionConfig

        config = NRReductionConfig()
        assert config.plot_save_dir is None

    def test_plot_save_dir_can_be_set(self):
        from lr_reduction.nr_reduction_config import NRReductionConfig

        config = NRReductionConfig()
        config.plot_save_dir = "/tmp/test_plots"
        assert config.plot_save_dir == "/tmp/test_plots"


class TestPlotReflectivity:
    """Test plot_reflectivity save and show behavior."""

    def _make_data_array(self):
        """Create minimal test data."""
        q = np.linspace(0.01, 0.1, 10)
        r = np.ones(10) * 0.5
        dr = np.ones(10) * 0.01
        dq = np.ones(10) * 0.001
        return [{"Q": q, "R": r, "dR": dr, "dQ": dq}]

    def test_plot_reflectivity_saves_file(self):
        from lr_reduction.new_reduction_from_template import plot_reflectivity

        data = self._make_data_array()
        with tempfile.TemporaryDirectory() as tmpdir:
            save_path = os.path.join(tmpdir, "test_plot.png")
            plot_reflectivity(data, save_path=save_path, show=False)
            assert os.path.isfile(save_path)

    @patch("lr_reduction.new_reduction_from_template.plt")
    def test_plot_reflectivity_no_show_when_false(self, mock_plt):
        from lr_reduction.new_reduction_from_template import plot_reflectivity

        # Set up mock so subplots returns fig, ax mocks
        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        data = self._make_data_array()
        plot_reflectivity(data, show=False)
        mock_plt.show.assert_not_called()
        mock_plt.close.assert_called_once_with(mock_fig)

    @patch("lr_reduction.new_reduction_from_template.plt")
    def test_plot_reflectivity_shows_when_true(self, mock_plt):
        from lr_reduction.new_reduction_from_template import plot_reflectivity

        mock_fig = MagicMock()
        mock_ax = MagicMock()
        mock_plt.subplots.return_value = (mock_fig, mock_ax)

        data = self._make_data_array()
        plot_reflectivity(data, show=True)
        mock_plt.show.assert_called_once()


class TestReduceFromTemplateProgressCallback:
    """Test that reduce_from_template accepts progress_callback parameter."""

    def test_signature_accepts_progress_callback(self):
        """Verify the function signature accepts progress_callback without error."""
        import inspect
        from lr_reduction.new_reduction_from_template import reduce_from_template

        sig = inspect.signature(reduce_from_template)
        assert "progress_callback" in sig.parameters
        assert sig.parameters["progress_callback"].default is None

    def test_signature_accepts_save_plots(self):
        import inspect
        from lr_reduction.new_reduction_from_template import reduce_from_template

        sig = inspect.signature(reduce_from_template)
        assert "save_plots" in sig.parameters
        assert sig.parameters["save_plots"].default is False

    def test_signature_accepts_plot_dir(self):
        import inspect
        from lr_reduction.new_reduction_from_template import reduce_from_template

        sig = inspect.signature(reduce_from_template)
        assert "plot_dir" in sig.parameters
        assert sig.parameters["plot_dir"].default is None


class TestShowOrSavePlot:
    """Test _show_or_save_plot helper in NR_Reduction."""

    def _make_config(self):
        from lr_reduction.nr_reduction_config import NRReductionConfig

        config = NRReductionConfig()
        config.RBnum = [1]
        config.DBname = ["test.dat"]
        config.RB_Ymin = [50]
        config.RB_Ymax = [150]
        config.tof_min = [1000]
        config.tof_max = [50000]
        config.ThetaShift = [0]
        config.ScaleFactor = [1.0]
        config.useBS = [0]
        config.BkgROI = [[0, 0, 0, 0]]
        config.plotON = False
        return config

    def test_show_or_save_saves_file(self):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from lr_reduction.nr_reduction_calc import NR_Reduction

        config = self._make_config()
        with tempfile.TemporaryDirectory() as tmpdir:
            config.plot_save_dir = tmpdir
            calc = NR_Reduction(config)
            fig, ax = plt.subplots()
            ax.plot([1, 2, 3], [1, 2, 3])
            calc._show_or_save_plot(fig, "test_diagnostic")
            assert os.path.isfile(os.path.join(tmpdir, "test_diagnostic.png"))

    @patch("lr_reduction.nr_reduction_calc.plt")
    def test_show_or_save_no_show_when_plotON_false(self, mock_plt):
        from lr_reduction.nr_reduction_calc import NR_Reduction

        config = self._make_config()
        config.plotON = False
        config.plot_save_dir = None
        calc = NR_Reduction(config)
        mock_fig = MagicMock()
        calc._show_or_save_plot(mock_fig, "test")
        mock_plt.show.assert_not_called()
        mock_plt.close.assert_called_once_with(mock_fig)
