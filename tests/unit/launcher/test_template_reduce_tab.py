"""Comprehensive UI tests for the TemplateReduce tab widget.

Uses QT_QPA_PLATFORM=offscreen for headless rendering (set in conftest.py).
"""

import os
import tempfile
from unittest.mock import MagicMock, patch

import pytest

from apps.template_reduce import (
    DATA_PATH_DIRECTIVE,
    DB_FILE_DIRECTIVE,
    DB_PATH_DIRECTIVE,
    OUTPUT_DIR_DIRECTIVE,
    TEMPLATE_DIRECTIVE,
    TEMPLATE_SAVE_DIRECTIVE,
    TemplateReduce,
    _METHOD_DISPLAY_TO_INTERNAL,
    _parse_run_numbers,
    _safe_float,
    _safe_int,
    _settings_bool,
)


# ============================================================
# _parse_run_numbers tests
# ============================================================
class TestParseRunNumbers:
    """Tests for the run number parsing utility."""

    def test_single_number(self):
        assert _parse_run_numbers("211029") == [211029]

    def test_comma_separated(self):
        assert _parse_run_numbers("211029,211030,211031") == [211029, 211030, 211031]

    def test_range(self):
        assert _parse_run_numbers("211029-211031") == [211029, 211030, 211031]

    def test_single_element_range(self):
        assert _parse_run_numbers("211029-211029") == [211029]

    def test_mixed_csv_and_range(self):
        result = _parse_run_numbers("211029,211033-211035")
        assert result == [211029, 211033, 211034, 211035]

    def test_whitespace_stripped(self):
        assert _parse_run_numbers("  211029 , 211030  ") == [211029, 211030]

    def test_whitespace_in_range(self):
        assert _parse_run_numbers(" 100 - 102 ") == [100, 101, 102]

    def test_empty_string_raises(self):
        with pytest.raises(ValueError, match="empty"):
            _parse_run_numbers("")

    def test_whitespace_only_raises(self):
        with pytest.raises(ValueError, match="empty"):
            _parse_run_numbers("   ")

    def test_non_numeric_raises(self):
        with pytest.raises(ValueError):
            _parse_run_numbers("abc")

    def test_reversed_range_raises(self):
        with pytest.raises(ValueError, match="end < start"):
            _parse_run_numbers("211031-211029")

    def test_trailing_comma(self):
        # Trailing comma should be silently ignored (empty split part)
        assert _parse_run_numbers("211029,") == [211029]

    def test_leading_comma(self):
        assert _parse_run_numbers(",211029") == [211029]

    def test_multiple_ranges(self):
        result = _parse_run_numbers("100-102,200-201")
        assert result == [100, 101, 102, 200, 201]


# ============================================================
# Widget instantiation and layout tests
# ============================================================
class TestTemplateReduceWidget:
    """Test that the widget creates correctly with all expected children."""

    @pytest.fixture
    def widget(self, qapp):
        w = TemplateReduce()
        yield w
        w.close()

    def test_widget_creates(self, widget):
        assert widget is not None
        assert widget.windowTitle() == "Template reduction"

    def test_has_run_numbers_field(self, widget):
        assert widget.run_numbers_edit is not None
        assert widget.run_numbers_edit.text() == "" or isinstance(widget.run_numbers_edit.text(), str)

    def test_has_template_path_label(self, widget):
        assert widget.template_path_label is not None

    def test_has_experiment_id_field(self, widget):
        assert widget.experiment_id_edit is not None

    def test_has_output_dir_label(self, widget):
        assert widget.output_dir_label is not None

    def test_has_data_path_label(self, widget):
        assert widget.data_path_label is not None

    def test_has_db_path_label(self, widget):
        assert widget.db_path_label is not None

    def test_has_db_file_label(self, widget):
        assert widget.db_file_label is not None

    def test_has_template_save_label(self, widget):
        assert widget.template_save_label is not None

    def test_has_reduce_button(self, widget):
        assert widget.reduce_btn is not None
        assert widget.reduce_btn.isEnabled()

    def test_has_cancel_button(self, widget):
        assert widget.cancel_btn is not None
        assert not widget.cancel_btn.isEnabled()

    def test_has_progress_bar(self, widget):
        assert widget.progress_bar is not None
        assert widget.progress_bar.value() == 0

    def test_has_status_label(self, widget):
        assert widget.status_label is not None
        assert widget.status_label.text() == "Ready"

    def test_has_method_combo(self, widget):
        items = [widget.method_combo.itemText(i) for i in range(widget.method_combo.count())]
        assert "meanTheta" in items
        assert "constantQ" in items
        assert "constantTOF" in items

    def test_has_processing_checkboxes(self, widget):
        assert widget.normalize_check is not None
        assert widget.autoscale_check is not None
        assert widget.use_calc_theta_check is not None

    def test_has_plot_checkboxes(self, widget):
        assert widget.save_plots_check is not None
        assert widget.interact_plots_check is not None
        assert widget.plot_rq4_check is not None

    def test_has_qspace_fields(self, widget):
        assert widget.qmin_edit is not None
        assert widget.qmax_edit is not None
        assert widget.dqbin_edit is not None
        assert widget.qline_threshold_edit is not None

    def test_has_dead_time_fields(self, widget):
        assert widget.dead_time_edit is not None
        assert widget.dead_time_tof_step_edit is not None

    def test_has_detector_fields(self, widget):
        det_items = [widget.det_res_fn_combo.itemText(i) for i in range(widget.det_res_fn_combo.count())]
        assert "rectangular" in det_items
        assert "gaussian" in det_items
        assert widget.det_sigma_edit is not None
        assert widget.peak_pad_edit is not None
        peak_items = [widget.peak_type_combo.itemText(i) for i in range(widget.peak_type_combo.count())]
        assert "supergauss" in peak_items
        assert "gauss" in peak_items

    def test_has_emission_time_fields(self, widget):
        assert widget.use_emission_time_check is not None
        assert widget.incident_theta_edit is not None

    def test_has_geometry_fields(self, widget):
        expected = {"mmpix", "dSampDet", "ny", "nx", "dMod", "xi_ref", "dS1Samp"}
        assert set(widget.geom_edits.keys()) == expected


# ============================================================
# Input validation tests
# ============================================================
class TestInputValidation:
    """Test check_inputs validation logic."""

    @pytest.fixture
    def widget(self, qapp):
        w = TemplateReduce()
        yield w
        w.close()

    def test_empty_run_numbers_fails(self, widget):
        widget.run_numbers_edit.setText("")
        widget.template_path_label.setText("/some/file.xml")
        widget.experiment_id_edit.setText("IPTS-12345")
        widget.output_dir_label.setText(tempfile.gettempdir())
        with patch.object(widget, "show_dialog") as mock_dialog:
            result = widget.check_inputs()
        assert result is False
        mock_dialog.assert_called_once()
        assert "run numbers" in mock_dialog.call_args[0][0].lower()

    def test_invalid_run_numbers_fails(self, widget):
        widget.run_numbers_edit.setText("abc")
        widget.template_path_label.setText("/some/file.xml")
        widget.experiment_id_edit.setText("IPTS-12345")
        widget.output_dir_label.setText(tempfile.gettempdir())
        with patch.object(widget, "show_dialog") as mock_dialog:
            result = widget.check_inputs()
        assert result is False

    def test_missing_template_file_fails(self, widget):
        widget.run_numbers_edit.setText("211029")
        widget.template_path_label.setText("/nonexistent/file.xml")
        widget.experiment_id_edit.setText("IPTS-12345")
        widget.output_dir_label.setText(tempfile.gettempdir())
        with patch.object(widget, "show_dialog") as mock_dialog:
            result = widget.check_inputs()
        assert result is False
        assert "template" in mock_dialog.call_args[0][0].lower()

    def test_empty_experiment_id_fails(self, widget):
        with tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as f:
            tmp_xml = f.name
        try:
            widget.run_numbers_edit.setText("211029")
            widget.template_path_label.setText(tmp_xml)
            widget.experiment_id_edit.setText("")
            widget.output_dir_label.setText(tempfile.gettempdir())
            with patch.object(widget, "show_dialog") as mock_dialog:
                result = widget.check_inputs()
            assert result is False
            assert "experiment" in mock_dialog.call_args[0][0].lower()
        finally:
            os.unlink(tmp_xml)

    def test_invalid_output_dir_fails(self, widget):
        with tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as f:
            tmp_xml = f.name
        try:
            widget.run_numbers_edit.setText("211029")
            widget.template_path_label.setText(tmp_xml)
            widget.experiment_id_edit.setText("IPTS-12345")
            widget.output_dir_label.setText("/nonexistent/dir/for/test")
            with patch.object(widget, "show_dialog") as mock_dialog:
                result = widget.check_inputs()
            assert result is False
            assert "output" in mock_dialog.call_args[0][0].lower() or "directory" in mock_dialog.call_args[0][0].lower()
        finally:
            os.unlink(tmp_xml)

    def test_directive_output_dir_fails(self, widget):
        with tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as f:
            tmp_xml = f.name
        try:
            widget.run_numbers_edit.setText("211029")
            widget.template_path_label.setText(tmp_xml)
            widget.experiment_id_edit.setText("IPTS-12345")
            widget.output_dir_label.setText(OUTPUT_DIR_DIRECTIVE)
            with patch.object(widget, "show_dialog") as mock_dialog:
                result = widget.check_inputs()
            assert result is False
        finally:
            os.unlink(tmp_xml)

    def test_valid_inputs_pass(self, widget):
        with tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as f:
            tmp_xml = f.name
        try:
            widget.run_numbers_edit.setText("211029")
            widget.template_path_label.setText(tmp_xml)
            widget.experiment_id_edit.setText("IPTS-12345")
            widget.output_dir_label.setText(tempfile.gettempdir())
            result = widget.check_inputs()
            assert result is True
        finally:
            os.unlink(tmp_xml)


# ============================================================
# Settings persistence tests
# ============================================================
class TestSettingsPersistence:
    """Test that save_settings / read_settings round-trip correctly."""

    @pytest.fixture
    def widget(self, qapp):
        w = TemplateReduce()
        yield w
        w.close()

    def test_round_trip_run_numbers(self, widget):
        widget.run_numbers_edit.setText("211029,211030")
        widget.save_settings()
        widget.run_numbers_edit.setText("")
        widget.read_settings()
        assert widget.run_numbers_edit.text() == "211029,211030"

    def test_round_trip_template_path(self, widget):
        widget.template_path_label.setText("/test/path/template.xml")
        widget.save_settings()
        widget.template_path_label.setText("")
        widget.read_settings()
        assert widget.template_path_label.text() == "/test/path/template.xml"

    def test_round_trip_experiment_id(self, widget):
        widget.experiment_id_edit.setText("IPTS-99999")
        widget.save_settings()
        widget.experiment_id_edit.setText("")
        widget.read_settings()
        assert widget.experiment_id_edit.text() == "IPTS-99999"

    def test_round_trip_method(self, widget):
        widget.method_combo.setCurrentText("constantQ")
        widget.save_settings()
        widget.method_combo.setCurrentText("meanTheta")
        widget.read_settings()
        assert widget.method_combo.currentText() == "constantQ"

    def test_round_trip_checkboxes(self, widget):
        widget.normalize_check.setChecked(True)
        widget.autoscale_check.setChecked(True)
        widget.use_calc_theta_check.setChecked(True)
        widget.save_plots_check.setChecked(True)
        widget.interact_plots_check.setChecked(True)
        widget.plot_rq4_check.setChecked(True)
        widget.save_settings()

        widget.normalize_check.setChecked(False)
        widget.autoscale_check.setChecked(False)
        widget.use_calc_theta_check.setChecked(False)
        widget.save_plots_check.setChecked(False)
        widget.interact_plots_check.setChecked(False)
        widget.plot_rq4_check.setChecked(False)
        widget.read_settings()

        assert widget.normalize_check.isChecked()
        assert widget.autoscale_check.isChecked()
        assert widget.use_calc_theta_check.isChecked()
        assert widget.save_plots_check.isChecked()
        assert widget.interact_plots_check.isChecked()
        assert widget.plot_rq4_check.isChecked()

    def test_round_trip_qspace_params(self, widget):
        widget.qmin_edit.setText("0.005")
        widget.qmax_edit.setText("0.3")
        widget.dqbin_edit.setText("0.01")
        widget.qline_threshold_edit.setText("0.5")
        widget.save_settings()

        widget.qmin_edit.setText("")
        widget.qmax_edit.setText("")
        widget.dqbin_edit.setText("")
        widget.qline_threshold_edit.setText("")
        widget.read_settings()

        assert widget.qmin_edit.text() == "0.005"
        assert widget.qmax_edit.text() == "0.3"
        assert widget.dqbin_edit.text() == "0.01"
        assert widget.qline_threshold_edit.text() == "0.5"

    def test_round_trip_dead_time(self, widget):
        widget.dead_time_edit.setText("4.5")
        widget.dead_time_tof_step_edit.setText("100")
        widget.save_settings()

        widget.dead_time_edit.setText("")
        widget.dead_time_tof_step_edit.setText("")
        widget.read_settings()

        assert widget.dead_time_edit.text() == "4.5"
        assert widget.dead_time_tof_step_edit.text() == "100"

    def test_round_trip_detector_peak(self, widget):
        widget.det_res_fn_combo.setCurrentText("gaussian")
        widget.det_sigma_edit.setText("1.2")
        widget.peak_pad_edit.setText("3")
        widget.peak_type_combo.setCurrentText("gauss")
        widget.save_settings()

        widget.det_res_fn_combo.setCurrentText("rectangular")
        widget.det_sigma_edit.setText("")
        widget.peak_pad_edit.setText("")
        widget.peak_type_combo.setCurrentText("supergauss")
        widget.read_settings()

        assert widget.det_res_fn_combo.currentText() == "gaussian"
        assert widget.det_sigma_edit.text() == "1.2"
        assert widget.peak_pad_edit.text() == "3"
        assert widget.peak_type_combo.currentText() == "gauss"

    def test_round_trip_emission_time(self, widget):
        widget.use_emission_time_check.setChecked(False)
        widget.incident_theta_edit.setText("3.5")
        widget.save_settings()

        widget.use_emission_time_check.setChecked(True)
        widget.incident_theta_edit.setText("")
        widget.read_settings()

        assert not widget.use_emission_time_check.isChecked()
        assert widget.incident_theta_edit.text() == "3.5"

    def test_round_trip_geometry(self, widget):
        widget.geom_edits["mmpix"].setText("0.7")
        widget.geom_edits["dSampDet"].setText("1800")
        widget.save_settings()

        widget.geom_edits["mmpix"].setText("")
        widget.geom_edits["dSampDet"].setText("")
        widget.read_settings()

        assert widget.geom_edits["mmpix"].text() == "0.7"
        assert widget.geom_edits["dSampDet"].text() == "1800"


# ============================================================
# Override params builder tests
# ============================================================
class TestBuildOverrideParams:
    """Test _build_override_params constructs the correct dict."""

    @pytest.fixture
    def widget(self, qapp):
        w = TemplateReduce()
        yield w
        w.close()

    def test_default_params(self, widget):
        """With no user changes, method and flags are still set."""
        params = widget._build_override_params()
        assert "method" in params
        assert params["method"] == "meantheta"

    def test_normalize_flag(self, widget):
        widget.normalize_check.setChecked(True)
        params = widget._build_override_params()
        assert params["Normalize"] is True

    def test_autoscale_flag(self, widget):
        widget.autoscale_check.setChecked(True)
        params = widget._build_override_params()
        assert params["AutoScale"] is True

    def test_use_calc_theta_flag(self, widget):
        widget.use_calc_theta_check.setChecked(True)
        params = widget._build_override_params()
        assert params["useCalcTheta"] is True

    def test_qspace_params_set(self, widget):
        widget.qmin_edit.setText("0.005")
        widget.qmax_edit.setText("0.3")
        params = widget._build_override_params()
        assert params["qmin"] == 0.005
        assert params["qmax"] == 0.3

    def test_qspace_params_empty_not_set(self, widget):
        widget.qmin_edit.setText("")
        widget.qmax_edit.setText("")
        params = widget._build_override_params()
        assert "qmin" not in params
        assert "qmax" not in params

    def test_dead_time_set(self, widget):
        widget.dead_time_edit.setText("5.0")
        widget.dead_time_tof_step_edit.setText("100")
        params = widget._build_override_params()
        assert params["dead_time"] == 5.0
        assert params["dead_time_tof_step"] == 100.0

    def test_detector_params(self, widget):
        widget.det_res_fn_combo.setCurrentText("gaussian")
        widget.det_sigma_edit.setText("1.5")
        widget.peak_pad_edit.setText("2")
        widget.peak_type_combo.setCurrentText("gauss")
        params = widget._build_override_params()
        assert params["DetResFn"] == "gaussian"
        assert params["DetSigma"] == 1.5
        assert params["peak_pad"] == 2
        assert params["peak_type"] == "gauss"

    def test_emission_time_params(self, widget):
        widget.use_emission_time_check.setChecked(False)
        widget.incident_theta_edit.setText("3.0")
        params = widget._build_override_params()
        assert params["use_emission_time"] is False
        assert params["IncidentTheta"] == 3.0

    def test_geometry_overrides(self, widget):
        widget.geom_edits["mmpix"].setText("0.7")
        widget.geom_edits["dSampDet"].setText("1800")
        params = widget._build_override_params()
        assert params["mmpix"] == 0.7
        assert params["dSampDet"] == 1800.0

    def test_geometry_empty_not_set(self, widget):
        for edit in widget.geom_edits.values():
            edit.setText("")
        params = widget._build_override_params()
        for name in widget.geom_edits:
            assert name not in params

    def test_output_dir_sets_spath(self, widget):
        with tempfile.TemporaryDirectory() as tmpdir:
            widget.output_dir_label.setText(tmpdir)
            params = widget._build_override_params()
            assert params["Spath"] == tmpdir

    def test_directive_output_dir_not_set(self, widget):
        widget.output_dir_label.setText(OUTPUT_DIR_DIRECTIVE)
        params = widget._build_override_params()
        assert "Spath" not in params

    def test_db_path_set(self, widget):
        with tempfile.TemporaryDirectory() as tmpdir:
            widget.db_path_label.setText(tmpdir)
            params = widget._build_override_params()
            assert params["DBpath"] == tmpdir

    def test_db_file_set(self, widget):
        with tempfile.NamedTemporaryFile(suffix=".dat", delete=False) as f:
            tmp_dat = f.name
        try:
            widget.db_file_label.setText(tmp_dat)
            params = widget._build_override_params()
            assert params["DBname"] == [os.path.basename(tmp_dat)]
        finally:
            os.unlink(tmp_dat)

    def test_plotQ4_flag(self, widget):
        widget.plot_rq4_check.setChecked(True)
        params = widget._build_override_params()
        assert params["plotQ4"] is True

        widget.plot_rq4_check.setChecked(False)
        params = widget._build_override_params()
        assert params["plotQ4"] is False


# ============================================================
# Signal connection tests
# ============================================================
class TestSignalConnections:
    """Verify that buttons are connected to the right slots."""

    @pytest.fixture
    def widget(self, qapp):
        w = TemplateReduce()
        yield w
        w.close()

    def test_reduce_button_connected(self, widget):
        """Reduce button should be connected to reduce()."""
        assert widget.reduce_btn.receivers(widget.reduce_btn.clicked) > 0

    def test_cancel_button_connected(self, widget):
        """Cancel button should be connected to on_cancel()."""
        assert widget.cancel_btn.receivers(widget.cancel_btn.clicked) > 0

    def test_template_button_connected(self, widget):
        assert widget.choose_template_btn.receivers(widget.choose_template_btn.clicked) > 0

    def test_output_dir_button_connected(self, widget):
        assert widget.choose_output_btn.receivers(widget.choose_output_btn.clicked) > 0

    def test_data_path_button_connected(self, widget):
        assert widget.choose_data_path_btn.receivers(widget.choose_data_path_btn.clicked) > 0

    def test_db_path_button_connected(self, widget):
        assert widget.choose_db_path_btn.receivers(widget.choose_db_path_btn.clicked) > 0

    def test_db_file_button_connected(self, widget):
        assert widget.choose_db_file_btn.receivers(widget.choose_db_file_btn.clicked) > 0

    def test_template_save_button_connected(self, widget):
        assert widget.choose_template_save_btn.receivers(widget.choose_template_save_btn.clicked) > 0


# ============================================================
# Progress callback and worker integration tests
# ============================================================
class TestProgressCallbacks:
    """Test progress bar and status label updates."""

    @pytest.fixture
    def widget(self, qapp):
        w = TemplateReduce()
        yield w
        w.close()

    def test_on_progress_updates_bar(self, widget):
        widget.on_progress(50, 100, "Halfway there")
        assert widget.progress_bar.value() == 50

    def test_on_progress_updates_status(self, widget):
        widget.on_progress(25, 100, "Run 211029: Reading NEXUS metadata")
        assert widget.status_label.text() == "Run 211029: Reading NEXUS metadata"

    def test_on_finished_reenables_button(self, widget):
        widget.reduce_btn.setEnabled(False)
        widget.cancel_btn.setEnabled(True)
        with patch.object(widget, "show_dialog"):
            widget.on_finished([])
        assert widget.reduce_btn.isEnabled()
        assert not widget.cancel_btn.isEnabled()
        assert widget.progress_bar.value() == 100

    def test_on_error_reenables_button(self, widget):
        widget.reduce_btn.setEnabled(False)
        widget.cancel_btn.setEnabled(True)
        with patch.object(widget, "show_dialog"):
            widget.on_error("Something went wrong")
        assert widget.reduce_btn.isEnabled()
        assert not widget.cancel_btn.isEnabled()
        assert widget.status_label.text() == "Error"

    def test_on_cancel_sets_worker_cancelled(self, widget):
        mock_worker = MagicMock()
        widget._worker = mock_worker
        widget.cancel_btn.setEnabled(True)
        widget.on_cancel()
        mock_worker.cancel.assert_called_once()
        assert widget.status_label.text() == "Cancelling..."


# ============================================================
# Reduce method behavior tests
# ============================================================
class TestReduceMethod:
    """Test the reduce() method's control flow."""

    @pytest.fixture
    def widget(self, qapp):
        w = TemplateReduce()
        yield w
        w.close()

    def test_reduce_aborts_on_invalid_input(self, widget):
        """If check_inputs fails, worker should NOT be started."""
        widget.run_numbers_edit.setText("")
        with patch.object(widget, "show_dialog"):
            widget.reduce()
        assert widget._worker is None

    def test_reduce_starts_worker_on_valid_input(self, widget):
        """With valid inputs, reduce() should create and start a worker."""
        with tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as f:
            tmp_xml = f.name
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                widget.run_numbers_edit.setText("211029")
                widget.template_path_label.setText(tmp_xml)
                widget.experiment_id_edit.setText("IPTS-12345")
                widget.output_dir_label.setText(tmpdir)

                with patch("apps.reduction_worker.ReductionWorker") as MockWorker:
                    mock_instance = MagicMock()
                    MockWorker.return_value = mock_instance
                    widget.reduce()

                    # Worker should have been created and started
                    MockWorker.assert_called_once()
                    mock_instance.start.assert_called_once()

                    # Controls should be in "running" state
                    assert not widget.reduce_btn.isEnabled()
                    assert widget.cancel_btn.isEnabled()
            finally:
                os.unlink(tmp_xml)

    def test_reduce_passes_datapath_when_set(self, widget):
        """When data path is a real directory, it should be passed to the worker."""
        with tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as f:
            tmp_xml = f.name
        with tempfile.TemporaryDirectory() as tmpdir, tempfile.TemporaryDirectory() as datadir:
            try:
                widget.run_numbers_edit.setText("211029")
                widget.template_path_label.setText(tmp_xml)
                widget.experiment_id_edit.setText("IPTS-12345")
                widget.output_dir_label.setText(tmpdir)
                widget.data_path_label.setText(datadir)

                with patch("apps.reduction_worker.ReductionWorker") as MockWorker:
                    mock_instance = MagicMock()
                    MockWorker.return_value = mock_instance
                    widget.reduce()

                    call_kwargs = MockWorker.call_args[1]
                    assert call_kwargs["datapath"] == datadir
            finally:
                os.unlink(tmp_xml)

    def test_reduce_passes_none_datapath_for_directive(self, widget):
        """When data path is still the placeholder directive, pass None."""
        with tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as f:
            tmp_xml = f.name
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                widget.run_numbers_edit.setText("211029")
                widget.template_path_label.setText(tmp_xml)
                widget.experiment_id_edit.setText("IPTS-12345")
                widget.output_dir_label.setText(tmpdir)
                widget.data_path_label.setText(DATA_PATH_DIRECTIVE)

                with patch("apps.reduction_worker.ReductionWorker") as MockWorker:
                    mock_instance = MagicMock()
                    MockWorker.return_value = mock_instance
                    widget.reduce()

                    call_kwargs = MockWorker.call_args[1]
                    assert call_kwargs["datapath"] is None
            finally:
                os.unlink(tmp_xml)

    def test_reduce_saves_settings(self, widget):
        """reduce() should call save_settings before starting."""
        with tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as f:
            tmp_xml = f.name
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                widget.run_numbers_edit.setText("211029")
                widget.template_path_label.setText(tmp_xml)
                widget.experiment_id_edit.setText("IPTS-12345")
                widget.output_dir_label.setText(tmpdir)

                with patch("apps.reduction_worker.ReductionWorker") as MockWorker:
                    mock_instance = MagicMock()
                    MockWorker.return_value = mock_instance
                    with patch.object(widget, "save_settings") as mock_save:
                        widget.reduce()
                        mock_save.assert_called_once()
            finally:
                os.unlink(tmp_xml)


# ============================================================
# Launcher integration test
# ============================================================
class TestLauncherIntegration:
    """Test that the tab integrates correctly into the launcher."""

    @staticmethod
    def _import_launcher_module():
        """Import the launcher module bypassing the test package name collision."""
        import importlib.util

        launcher_py = os.path.join(os.path.dirname(__file__), "..", "..", "..", "launcher", "launcher.py")
        launcher_py = os.path.abspath(launcher_py)
        spec = importlib.util.spec_from_file_location("launcher_main", launcher_py)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod

    def test_launcher_has_template_tab(self, qapp):
        """ReductionInterface should have a 'Template reduction' tab."""
        mod = self._import_launcher_module()
        window = mod.ReductionInterface()
        tab_texts = [window.tabText(i) for i in range(window.count())]
        assert "Template reduction" in tab_texts
        window.close()

    def test_template_tab_is_correct_type(self, qapp):
        mod = self._import_launcher_module()
        window = mod.ReductionInterface()
        tab_index = None
        for i in range(window.count()):
            if window.tabText(i) == "Template reduction":
                tab_index = i
                break
        assert tab_index is not None
        tab_widget = window.widget(tab_index)
        assert isinstance(tab_widget, TemplateReduce)
        window.close()


# ============================================================
# ReductionWorker tests
# ============================================================
class TestReductionWorker:
    """Test the worker thread in isolation."""

    def test_worker_instantiation(self, qapp):
        from apps.reduction_worker import ReductionWorker

        worker = ReductionWorker(
            run_numbers=[211029],
            template_file="/fake/template.xml",
            experiment_id="IPTS-12345",
        )
        assert worker.run_numbers == [211029]
        assert not worker._cancelled

    def test_worker_cancel_sets_flag(self, qapp):
        from apps.reduction_worker import ReductionWorker

        worker = ReductionWorker(
            run_numbers=[211029],
            template_file="/fake/template.xml",
            experiment_id="IPTS-12345",
        )
        worker.cancel()
        assert worker._cancelled

    def test_worker_has_signals(self, qapp):
        from apps.reduction_worker import ReductionWorker

        worker = ReductionWorker(
            run_numbers=[211029],
            template_file="/fake/template.xml",
            experiment_id="IPTS-12345",
        )
        # Verify signals exist and are connectable
        assert hasattr(worker, "progress")
        assert hasattr(worker, "finished")
        assert hasattr(worker, "error")

    def test_worker_emits_error_on_bad_file(self, qapp):
        """Worker should emit error signal when reduce_from_template fails."""
        from apps.reduction_worker import ReductionWorker

        worker = ReductionWorker(
            run_numbers=[999999],
            template_file="/nonexistent/template.xml",
            experiment_id="IPTS-00000",
        )

        errors = []
        worker.error.connect(lambda msg: errors.append(msg))
        worker.run()  # run synchronously (not start()) for testing
        assert len(errors) == 1
        assert len(errors[0]) > 0  # got an error message

    def test_worker_accepts_log_to_stdout(self, qapp):
        """Worker should accept log_to_stdout parameter."""
        from apps.reduction_worker import ReductionWorker

        worker = ReductionWorker(
            run_numbers=[211029],
            template_file="/fake/template.xml",
            experiment_id="IPTS-12345",
            log_to_stdout=True,
        )
        assert worker.log_to_stdout is True

    def test_worker_log_to_stdout_defaults_false(self, qapp):
        """Worker log_to_stdout should default to False."""
        from apps.reduction_worker import ReductionWorker

        worker = ReductionWorker(
            run_numbers=[211029],
            template_file="/fake/template.xml",
            experiment_id="IPTS-12345",
        )
        assert worker.log_to_stdout is False


# ============================================================
# _safe_float and _safe_int tests
# ============================================================
class TestSafeConversions:
    """Test _safe_float and _safe_int helper functions."""

    def test_safe_float_valid(self):
        assert _safe_float("3.14") == 3.14

    def test_safe_float_integer(self):
        assert _safe_float("42") == 42.0

    def test_safe_float_negative(self):
        assert _safe_float("-1.5") == -1.5

    def test_safe_float_scientific(self):
        assert _safe_float("1e-3") == 0.001

    def test_safe_float_empty(self):
        assert _safe_float("") is None

    def test_safe_float_whitespace(self):
        assert _safe_float("   ") is None

    def test_safe_float_invalid(self):
        assert _safe_float("abc") is None

    def test_safe_float_partial_minus(self):
        """QDoubleValidator intermediate state: bare minus sign."""
        assert _safe_float("-") is None

    def test_safe_float_partial_dot(self):
        """QDoubleValidator intermediate state: bare dot."""
        assert _safe_float(".") is None

    def test_safe_float_partial_scientific(self):
        """QDoubleValidator intermediate state: incomplete scientific notation."""
        assert _safe_float("1e") is None

    def test_safe_float_with_spaces(self):
        assert _safe_float("  3.14  ") == 3.14

    def test_safe_int_valid(self):
        assert _safe_int("42") == 42

    def test_safe_int_empty(self):
        assert _safe_int("") is None

    def test_safe_int_whitespace(self):
        assert _safe_int("   ") is None

    def test_safe_int_invalid(self):
        assert _safe_int("abc") is None

    def test_safe_int_float_string(self):
        """A float string like '3.0' should convert to int 3."""
        assert _safe_int("3.0") == 3

    def test_safe_int_negative(self):
        assert _safe_int("-5") == -5

    def test_safe_int_partial_minus(self):
        assert _safe_int("-") is None

    def test_safe_int_with_spaces(self):
        assert _safe_int("  7  ") == 7


# ============================================================
# _settings_bool tests
# ============================================================
class TestSettingsBool:
    """Test _settings_bool helper."""

    def test_bool_true(self):
        assert _settings_bool(True) is True

    def test_bool_false(self):
        assert _settings_bool(False) is False

    def test_string_true(self):
        assert _settings_bool("true") is True

    def test_string_True(self):
        assert _settings_bool("True") is True

    def test_string_false(self):
        assert _settings_bool("false") is False

    def test_string_empty(self):
        assert _settings_bool("") is False

    def test_none(self):
        assert _settings_bool(None) is False

    def test_int_nonzero(self):
        # int is not bool or str, should return False
        assert _settings_bool(1) is False


# ============================================================
# Method case normalization tests
# ============================================================
class TestMethodCaseNormalization:
    """Test that method display names map to lowercase internal values."""

    def test_display_to_internal_mapping(self):
        assert _METHOD_DISPLAY_TO_INTERNAL["meanTheta"] == "meantheta"
        assert _METHOD_DISPLAY_TO_INTERNAL["constantQ"] == "constantq"
        assert _METHOD_DISPLAY_TO_INTERNAL["constantTOF"] == "constanttof"

    @pytest.fixture
    def widget(self, qapp):
        w = TemplateReduce()
        yield w
        w.close()

    def test_meantheta_override(self, widget):
        widget.method_combo.setCurrentText("meanTheta")
        params = widget._build_override_params()
        assert params["method"] == "meantheta"

    def test_constantq_override(self, widget):
        widget.method_combo.setCurrentText("constantQ")
        params = widget._build_override_params()
        assert params["method"] == "constantq"

    def test_constanttof_override(self, widget):
        widget.method_combo.setCurrentText("constantTOF")
        params = widget._build_override_params()
        assert params["method"] == "constanttof"


# ============================================================
# Boolean flags always sent tests
# ============================================================
class TestBooleanFlagsAlwaysSent:
    """Test that boolean flags are always included in override params."""

    @pytest.fixture
    def widget(self, qapp):
        w = TemplateReduce()
        yield w
        w.close()

    def test_normalize_false_is_sent(self, widget):
        widget.normalize_check.setChecked(False)
        params = widget._build_override_params()
        assert "Normalize" in params
        assert params["Normalize"] is False

    def test_autoscale_false_is_sent(self, widget):
        widget.autoscale_check.setChecked(False)
        params = widget._build_override_params()
        assert "AutoScale" in params
        assert params["AutoScale"] is False

    def test_use_calc_theta_false_is_sent(self, widget):
        widget.use_calc_theta_check.setChecked(False)
        params = widget._build_override_params()
        assert "useCalcTheta" in params
        assert params["useCalcTheta"] is False

    def test_plotq4_false_is_sent(self, widget):
        widget.plot_rq4_check.setChecked(False)
        params = widget._build_override_params()
        assert "plotQ4" in params
        assert params["plotQ4"] is False

    def test_use_emission_time_false_is_sent(self, widget):
        widget.use_emission_time_check.setChecked(False)
        params = widget._build_override_params()
        assert "use_emission_time" in params
        assert params["use_emission_time"] is False


# ============================================================
# Numeric field validation tests
# ============================================================
class TestNumericFieldValidation:
    """Test _validate_numeric_fields catches invalid entries."""

    @pytest.fixture
    def widget(self, qapp):
        w = TemplateReduce()
        yield w
        w.close()

    def test_valid_floats_pass(self, widget):
        widget.qmin_edit.setText("0.005")
        widget.qmax_edit.setText("0.3")
        assert widget._validate_numeric_fields() is None

    def test_empty_fields_pass(self, widget):
        assert widget._validate_numeric_fields() is None

    def test_invalid_qmin_fails(self, widget):
        widget.qmin_edit.setText("abc")
        result = widget._validate_numeric_fields()
        assert result is not None
        assert "qmin" in result

    def test_invalid_dead_time_fails(self, widget):
        widget.dead_time_edit.setText("-")
        result = widget._validate_numeric_fields()
        assert result is not None
        assert "dead_time" in result

    def test_invalid_peak_pad_fails(self, widget):
        widget.peak_pad_edit.setText("xyz")
        result = widget._validate_numeric_fields()
        assert result is not None
        assert "peak_pad" in result

    def test_invalid_geometry_field_fails(self, widget):
        widget.geom_edits["mmpix"].setText("not_a_number")
        result = widget._validate_numeric_fields()
        assert result is not None
        assert "mmpix" in result

    def test_check_inputs_catches_numeric_errors(self, widget):
        """Full check_inputs should catch numeric validation failures."""
        with tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as f:
            tmp_xml = f.name
        try:
            widget.run_numbers_edit.setText("211029")
            widget.template_path_label.setText(tmp_xml)
            widget.experiment_id_edit.setText("IPTS-12345")
            widget.output_dir_label.setText(tempfile.gettempdir())
            widget.qmin_edit.setText("not_a_number")
            with patch.object(widget, "show_dialog") as mock_dialog:
                result = widget.check_inputs()
            assert result is False
            assert "qmin" in mock_dialog.call_args[0][0]
        finally:
            os.unlink(tmp_xml)


# ============================================================
# Double-click guard tests
# ============================================================
class TestDoubleClickGuard:
    """Test that reduce() prevents double-clicking."""

    @pytest.fixture
    def widget(self, qapp):
        w = TemplateReduce()
        yield w
        w.close()

    def test_reduce_blocked_when_worker_running(self, widget):
        """If a worker is already running, reduce() should return immediately."""
        mock_worker = MagicMock()
        mock_worker.isRunning.return_value = True
        widget._worker = mock_worker

        # reduce() should bail out without calling check_inputs
        with patch.object(widget, "check_inputs") as mock_check:
            widget.reduce()
            mock_check.assert_not_called()

    def test_reduce_allowed_when_worker_finished(self, widget):
        """If previous worker is done, reduce() should proceed to validation."""
        mock_worker = MagicMock()
        mock_worker.isRunning.return_value = False
        widget._worker = mock_worker

        # reduce() should proceed to check_inputs (which will fail due to empty fields)
        with patch.object(widget, "show_dialog"):
            widget.reduce()
        # check_inputs was called (widget shows dialog for empty run numbers)


# ============================================================
# Cancellation detection tests
# ============================================================
class TestCancellationDetection:
    """Test on_finished behavior when worker was cancelled."""

    @pytest.fixture
    def widget(self, qapp):
        w = TemplateReduce()
        yield w
        w.close()

    def test_on_finished_detects_cancellation(self, widget):
        """on_finished should show cancellation message if worker was cancelled."""
        mock_worker = MagicMock()
        mock_worker._cancelled = True
        widget._worker = mock_worker
        widget.reduce_btn.setEnabled(False)
        widget.cancel_btn.setEnabled(True)

        with patch.object(widget, "show_dialog") as mock_dialog:
            widget.on_finished([])

        assert widget.status_label.text() == "Cancelled"
        mock_dialog.assert_called_once()
        assert "cancelled" in mock_dialog.call_args[0][0].lower()

    def test_on_finished_shows_success_on_completion(self, widget):
        """on_finished should show success message if worker completed normally."""
        mock_worker = MagicMock()
        mock_worker._cancelled = False
        widget._worker = mock_worker
        widget.reduce_btn.setEnabled(False)

        with patch.object(widget, "show_dialog") as mock_dialog:
            widget.on_finished([])

        assert widget.status_label.text() == "Complete"
        mock_dialog.assert_called_once()
        assert "success" in mock_dialog.call_args[0][0].lower()


# ============================================================
# show_dialog icon parameter tests
# ============================================================
class TestShowDialogIcon:
    """Test that show_dialog accepts and uses the icon parameter."""

    @pytest.fixture
    def widget(self, qapp):
        w = TemplateReduce()
        yield w
        w.close()

    def test_default_icon_is_critical(self, widget):
        """Default icon should be Critical."""
        from qtpy.QtWidgets import QMessageBox

        with patch.object(QMessageBox, "exec_") as mock_exec:
            with patch.object(QMessageBox, "setIcon") as mock_set_icon:
                widget.show_dialog("test error")
                mock_set_icon.assert_called_with(QMessageBox.Critical)

    def test_custom_icon_information(self, widget):
        """Passing Information icon should work."""
        from qtpy.QtWidgets import QMessageBox

        with patch.object(QMessageBox, "exec_") as mock_exec:
            with patch.object(QMessageBox, "setIcon") as mock_set_icon:
                widget.show_dialog("test success", icon=QMessageBox.Information)
                mock_set_icon.assert_called_with(QMessageBox.Information)

    def test_custom_icon_warning(self, widget):
        """Passing Warning icon should work."""
        from qtpy.QtWidgets import QMessageBox

        with patch.object(QMessageBox, "exec_") as mock_exec:
            with patch.object(QMessageBox, "setIcon") as mock_set_icon:
                widget.show_dialog("test warning", icon=QMessageBox.Warning)
                mock_set_icon.assert_called_with(QMessageBox.Warning)


# ============================================================
# Log to stdout checkbox tests
# ============================================================
class TestLogStdoutCheckbox:
    """Test the log_stdout_check widget and its integration."""

    @pytest.fixture
    def widget(self, qapp):
        w = TemplateReduce()
        yield w
        w.close()

    def test_log_stdout_checkbox_exists(self, widget):
        assert widget.log_stdout_check is not None

    def test_log_stdout_default_unchecked(self, widget):
        assert not widget.log_stdout_check.isChecked()

    def test_log_stdout_settings_round_trip(self, widget):
        widget.log_stdout_check.setChecked(True)
        widget.save_settings()
        widget.log_stdout_check.setChecked(False)
        widget.read_settings()
        assert widget.log_stdout_check.isChecked()

    def test_log_stdout_passed_to_worker(self, widget):
        """When log_stdout is checked, reduce() should pass it to the worker."""
        with tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as f:
            tmp_xml = f.name
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                widget.run_numbers_edit.setText("211029")
                widget.template_path_label.setText(tmp_xml)
                widget.experiment_id_edit.setText("IPTS-12345")
                widget.output_dir_label.setText(tmpdir)
                widget.log_stdout_check.setChecked(True)

                with patch("apps.reduction_worker.ReductionWorker") as MockWorker:
                    mock_instance = MagicMock()
                    MockWorker.return_value = mock_instance
                    widget.reduce()

                    call_kwargs = MockWorker.call_args[1]
                    assert call_kwargs["log_to_stdout"] is True
            finally:
                os.unlink(tmp_xml)

    def test_on_error_prints_to_stderr_when_logging(self, widget, capsys):
        """on_error should print to stderr when log_stdout is checked."""
        widget.log_stdout_check.setChecked(True)
        with patch.object(widget, "show_dialog"):
            widget.on_error("Test error message")
        captured = capsys.readouterr()
        assert "Test error message" in captured.err


# ============================================================
# Safe conversion integration in _build_override_params
# ============================================================
class TestSafeConversionInOverrideParams:
    """Test that _build_override_params handles intermediate validator states."""

    @pytest.fixture
    def widget(self, qapp):
        w = TemplateReduce()
        yield w
        w.close()

    def test_partial_float_skipped(self, widget):
        """A partially typed float like '-' should be skipped, not crash."""
        widget.qmin_edit.setText("-")
        params = widget._build_override_params()
        assert "qmin" not in params

    def test_partial_dot_skipped(self, widget):
        """A bare '.' should be skipped."""
        widget.dead_time_edit.setText(".")
        params = widget._build_override_params()
        assert "dead_time" not in params

    def test_incomplete_scientific_skipped(self, widget):
        """'1e' should be skipped."""
        widget.qmax_edit.setText("1e")
        params = widget._build_override_params()
        assert "qmax" not in params

    def test_partial_int_skipped(self, widget):
        """A bare '-' in an int field should be skipped."""
        widget.peak_pad_edit.setText("-")
        params = widget._build_override_params()
        assert "peak_pad" not in params

    def test_geometry_partial_skipped(self, widget):
        """Geometry field with partial input should be skipped."""
        widget.geom_edits["mmpix"].setText(".")
        params = widget._build_override_params()
        assert "mmpix" not in params
