from lr_reduction.data_info import DataType


class _FakeProperty:
    def __init__(self, value):
        self.value = value


class _FakeRun:
    def __init__(self, logs):
        self._logs = {name: _FakeProperty(value) for name, value in logs.items()}

    def getProperty(self, key):  # noqa: N802
        return self._logs[key]

    def hasProperty(self, key):  # noqa: N802
        return key in self._logs


class _FakeWorkspace:
    def __init__(self, logs):
        self._run = _FakeRun(logs)

    def getRun(self):  # noqa: N802
        return self._run


def test_from_workspace_direct_beam_earth_centered_coordinates():
    """Test direct beam detection with earth-centered coordinates"""
    workspace = _FakeWorkspace(
        {
            "BL4B:CS:Mode:Coordinates": 0,
            "thi": 0.800,
            "tthd": 0.805,
        }
    )

    result = DataType.from_workspace(workspace)
    assert result == DataType.DIRECT_BEAM


def test_from_workspace_direct_beam_free_liquid_mode():
    """Test direct beam detection with free liquid mode"""
    workspace = _FakeWorkspace(
        {
            "BL4B:CS:ExpPl:OperatingMode": "Free Liquid",
            "thi": 0.210,
            "tthd": 0.212,
        }
    )

    result = DataType.from_workspace(workspace)
    assert result == DataType.DIRECT_BEAM


def test_from_workspace_coordinate_mode_overrides_legacy_operating_mode():
    """Beam-centered coordinate mode should win over the legacy Free Liquid fallback."""
    workspace = _FakeWorkspace(
        {
            "BL4B:CS:Mode:Coordinates": 1,
            "BL4B:CS:ExpPl:OperatingMode": "Free Liquid",
            "thi": 0.210,
            "ths": 0.0006,
            "tthd": 0.0008,
        }
    )

    result = DataType.from_workspace(workspace)
    assert result == DataType.DIRECT_BEAM


def test_from_workspace_reflected_beam_earth_centered():
    """Test reflected beam detection with earth-centered coordinates"""
    workspace = _FakeWorkspace(
        {
            "BL4B:CS:Mode:Coordinates": 0,
            "thi": 0.5,
            "tthd": 1.0,
        }
    )

    result = DataType.from_workspace(workspace)
    assert result == DataType.REFLECTED_BEAM


def test_from_workspace_reflected_beam_beam_centered():
    """Test reflected beam detection with beam-centered coordinates"""
    workspace = _FakeWorkspace(
        {
            "BL4B:CS:Mode:Coordinates": 1,
            "BL4B:CS:ExpPl:OperatingMode": "Other",
            "ths": 0.5,
            "tthd": 1.0,
        }
    )

    result = DataType.from_workspace(workspace)
    assert result == DataType.REFLECTED_BEAM


def test_from_workspace_direct_beam_beam_centered():
    """Test direct beam detection with beam-centered coordinates"""
    workspace = _FakeWorkspace(
        {
            "BL4B:CS:Mode:Coordinates": 1,
            "BL4B:CS:ExpPl:OperatingMode": "Other",
            "ths": 0.0006,
            "tthd": 0.0008,
        }
    )

    result = DataType.from_workspace(workspace)
    assert result == DataType.DIRECT_BEAM


def test_from_workspace_missing_logs():
    """Test that missing logs default to reflected beam"""
    workspace = _FakeWorkspace({})

    result = DataType.from_workspace(workspace)
    assert result == DataType.REFLECTED_BEAM
