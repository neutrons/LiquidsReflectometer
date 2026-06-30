import pytest

from lr_reduction.legacy.theta_selection import is_earth_centered_geometry, theta_log_name


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


@pytest.mark.parametrize(
    ("logs", "expected"),
    [
        ({"BL4B:CS:Mode:Coordinates": 0}, True),
        ({"BL4B:CS:Mode:Coordinates": 1, "BL4B:CS:ExpPl:OperatingMode": "Free Liquid"}, False),
        ({"BL4B:CS:Mode:Coordinates": 1, "BL4B:CS:ExpPl:OperatingMode": "Bound Liquid"}, False),
        ({"BL4B:CS:ExpPl:OperatingMode": "Free Liquid"}, True),
        ({}, False),
    ],
)
def test_is_earth_centered_geometry(logs, expected):
    assert is_earth_centered_geometry(_FakeWorkspace(logs)) is expected


@pytest.mark.parametrize(
    ("logs", "expected"),
    [
        ({"BL4B:CS:Mode:Coordinates": 0}, "thi"),
        ({"BL4B:CS:Mode:Coordinates": 1, "BL4B:CS:ExpPl:OperatingMode": "Free Liquid"}, "ths"),
        ({"BL4B:CS:ExpPl:OperatingMode": "Free Liquid"}, "thi"),
        ({"BL4B:CS:Mode:Coordinates": 1, "BL4B:CS:ExpPl:OperatingMode": "Bound Liquid"}, "ths"),
    ],
)
def test_theta_log_name(logs, expected):
    assert theta_log_name(_FakeWorkspace(logs)) == expected
