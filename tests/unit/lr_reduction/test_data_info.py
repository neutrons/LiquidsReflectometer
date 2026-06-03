from unittest.mock import Mock, patch

from lr_reduction.data_info import DataType


def _set_mock_logs(*mock_classes, logs):
    for mock_class in mock_classes:
        mock_class.return_value = logs


@patch("lr_reduction.theta_selection.SampleLogValues")
@patch("lr_reduction.data_info.SampleLogValues")
def test_from_workspace_direct_beam_earth_centered_coordinates(
    mock_data_info_sample_logs, mock_theta_selection_sample_logs
):
    """Test direct beam detection with earth-centered coordinates"""
    mock_logs = {
        "BL4B:CS:Mode:Coordinates": 0,
        "thi": 0.800,
        "tthd": 0.805,
    }
    _set_mock_logs(mock_data_info_sample_logs, mock_theta_selection_sample_logs, logs=mock_logs)

    result = DataType.from_workspace(Mock())
    assert result == DataType.DIRECT_BEAM


@patch("lr_reduction.theta_selection.SampleLogValues")
@patch("lr_reduction.data_info.SampleLogValues")
def test_from_workspace_direct_beam_free_liquid_mode(mock_data_info_sample_logs, mock_theta_selection_sample_logs):
    """Test direct beam detection with free liquid mode"""
    mock_logs = {
        "BL4B:CS:ExpPl:OperatingMode": "Free Liquid",
        "thi": 0.210,
        "tthd": 0.212,
    }
    _set_mock_logs(mock_data_info_sample_logs, mock_theta_selection_sample_logs, logs=mock_logs)

    result = DataType.from_workspace(Mock())
    assert result == DataType.DIRECT_BEAM


@patch("lr_reduction.theta_selection.SampleLogValues")
@patch("lr_reduction.data_info.SampleLogValues")
def test_from_workspace_coordinate_mode_overrides_legacy_operating_mode(
    mock_data_info_sample_logs, mock_theta_selection_sample_logs
):
    """Beam-centered coordinate mode should win over the legacy Free Liquid fallback."""
    mock_logs = {
        "BL4B:CS:Mode:Coordinates": 1,
        "BL4B:CS:ExpPl:OperatingMode": "Free Liquid",
        "thi": 0.210,
        "ths": 0.0006,
        "tthd": 0.0008,
    }
    _set_mock_logs(mock_data_info_sample_logs, mock_theta_selection_sample_logs, logs=mock_logs)

    result = DataType.from_workspace(Mock())
    assert result == DataType.DIRECT_BEAM


@patch("lr_reduction.theta_selection.SampleLogValues")
@patch("lr_reduction.data_info.SampleLogValues")
def test_from_workspace_reflected_beam_earth_centered(mock_data_info_sample_logs, mock_theta_selection_sample_logs):
    """Test reflected beam detection with earth-centered coordinates"""
    mock_logs = {
        "BL4B:CS:Mode:Coordinates": 0,
        "thi": 0.5,
        "tthd": 1.0,
    }
    _set_mock_logs(mock_data_info_sample_logs, mock_theta_selection_sample_logs, logs=mock_logs)

    result = DataType.from_workspace(Mock())
    assert result == DataType.REFLECTED_BEAM


@patch("lr_reduction.theta_selection.SampleLogValues")
@patch("lr_reduction.data_info.SampleLogValues")
def test_from_workspace_reflected_beam_beam_centered(
    mock_data_info_sample_logs, mock_theta_selection_sample_logs
):
    """Test reflected beam detection with beam-centered coordinates"""
    mock_logs = {
        "BL4B:CS:Mode:Coordinates": 1,
        "BL4B:CS:ExpPl:OperatingMode": "Other",
        "ths": 0.5,
        "tthd": 1.0,
    }
    _set_mock_logs(mock_data_info_sample_logs, mock_theta_selection_sample_logs, logs=mock_logs)

    result = DataType.from_workspace(Mock())
    assert result == DataType.REFLECTED_BEAM


@patch("lr_reduction.theta_selection.SampleLogValues")
@patch("lr_reduction.data_info.SampleLogValues")
def test_from_workspace_direct_beam_beam_centered(mock_data_info_sample_logs, mock_theta_selection_sample_logs):
    """Test direct beam detection with beam-centered coordinates"""
    mock_logs = {
        "BL4B:CS:Mode:Coordinates": 1,
        "BL4B:CS:ExpPl:OperatingMode": "Other",
        "ths": 0.0006,
        "tthd": 0.0008,
    }
    _set_mock_logs(mock_data_info_sample_logs, mock_theta_selection_sample_logs, logs=mock_logs)

    result = DataType.from_workspace(Mock())
    assert result == DataType.DIRECT_BEAM


@patch("lr_reduction.theta_selection.SampleLogValues")
@patch("lr_reduction.data_info.SampleLogValues")
def test_from_workspace_missing_logs(mock_data_info_sample_logs, mock_theta_selection_sample_logs):
    """Test that missing logs default to reflected beam"""
    _set_mock_logs(mock_data_info_sample_logs, mock_theta_selection_sample_logs, logs={})

    result = DataType.from_workspace(Mock())
    assert result == DataType.REFLECTED_BEAM
