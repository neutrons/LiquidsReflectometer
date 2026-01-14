from unittest.mock import Mock, patch

from lr_reduction.data_info import DataType


@patch("lr_reduction.data_info.SampleLogValues")
def test_from_workspace_direct_beam_earth_centered_coordinates(mock_sample_logs_class):
    """Test direct beam detection with earth-centered coordinates"""
    mock_logs = {
        "BL4B:CS:Mode:Coordinates": 0,
        "thi": 0.005,
        "tthd": 0.005,
    }
    mock_sample_logs_class.return_value = mock_logs

    result = DataType.from_workspace(Mock())
    assert result == DataType.DIRECT_BEAM


@patch("lr_reduction.data_info.SampleLogValues")
def test_from_workspace_direct_beam_free_liquid_mode(mock_sample_logs_class):
    """Test direct beam detection with free liquid mode"""
    mock_logs = {
        "BL4B:CS:Mode:Coordinates": 1,
        "BL4B:CS:ExpPl:OperatingMode": "Free Liquid",
        "thi": 0.0,
        "tthd": 0.0,
    }
    mock_sample_logs_class.return_value = mock_logs

    result = DataType.from_workspace(Mock())
    assert result == DataType.DIRECT_BEAM


@patch("lr_reduction.data_info.SampleLogValues")
def test_from_workspace_reflected_beam_earth_centered(mock_sample_logs_class):
    """Test reflected beam detection with earth-centered coordinates"""
    mock_logs = {
        "BL4B:CS:Mode:Coordinates": 0,
        "thi": 0.5,
        "tthd": 1.0,
    }
    mock_sample_logs_class.return_value = mock_logs

    result = DataType.from_workspace(Mock())
    assert result == DataType.REFLECTED_BEAM


@patch("lr_reduction.data_info.SampleLogValues")
def test_from_workspace_reflected_beam_beam_centered(mock_sample_logs_class):
    """Test reflected beam detection with beam-centered coordinates"""
    mock_logs = {
        "BL4B:CS:Mode:Coordinates": 1,
        "BL4B:CS:ExpPl:OperatingMode": "Other",
        "ths": 0.5,
        "tthd": 1.0,
    }
    mock_sample_logs_class.return_value = mock_logs

    result = DataType.from_workspace(Mock())
    assert result == DataType.REFLECTED_BEAM


@patch("lr_reduction.data_info.SampleLogValues")
def test_from_workspace_direct_beam_beam_centered(mock_sample_logs_class):
    """Test direct beam detection with beam-centered coordinates"""
    mock_logs = {
        "BL4B:CS:Mode:Coordinates": 1,
        "BL4B:CS:ExpPl:OperatingMode": "Other",
        "ths": 0.0,
        "tthd": 0.0,
    }
    mock_sample_logs_class.return_value = mock_logs

    result = DataType.from_workspace(Mock())
    assert result == DataType.DIRECT_BEAM


@patch("lr_reduction.data_info.SampleLogValues")
def test_from_workspace_missing_logs(mock_sample_logs_class):
    """Test that missing logs default to reflected beam"""
    mock_sample_logs_class.side_effect = KeyError("Missing log")

    result = DataType.from_workspace(Mock())
    assert result == DataType.REFLECTED_BEAM


@patch("lr_reduction.data_info.SampleLogValues")
def test_from_workspace_exception_handling(mock_sample_logs_class):
    """Test that exceptions default to reflected beam"""
    mock_sample_logs_class.side_effect = Exception("Unexpected error")

    result = DataType.from_workspace(Mock())
    assert result == DataType.REFLECTED_BEAM
