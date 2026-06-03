import math

import numpy as np
import pytest

from lr_reduction.reduction_template_reader import ReductionParameters
from lr_reduction.template import get_default_template_file, process_from_template_ws


def test_get_default_template_file_up_down_and_default(tmp_path):
    out = tmp_path

    # No template -> ValueError
    with pytest.raises(ValueError):
        get_default_template_file(str(out), 0)

    # Case: fallback default template
    default_file = out / "template.xml"
    default_file.touch()
    # tthd doesn't matter when only default exists
    assert get_default_template_file(str(out), 0) == str(default_file)

    # Case: up geometry
    up_file = out / "template_up.xml"
    up_file.touch()
    assert get_default_template_file(str(out), 1) == str(up_file)

    # Case: down geometry
    down_file = out / "template_down.xml"
    down_file.touch()
    assert get_default_template_file(str(out), -1) == str(down_file)


class _FakeProperty:
    def __init__(self, value):
        self.value = [value]
        self.size = 1


class _FakeRun:
    def __init__(self, logs):
        self._logs = {name: _FakeProperty(value) for name, value in logs.items()}

    def __contains__(self, key):
        return key in self._logs

    def __getitem__(self, key):
        return self._logs[key]

    def getProperty(self, key):  # noqa: N802
        return self._logs[key]

    def hasProperty(self, key):  # noqa: N802
        return key in self._logs


class _FakeWorkspace:
    def __init__(self, logs):
        self._run = _FakeRun(logs)

    def getRun(self):  # noqa: N802
        return self._run


class _DummyEventReflectivity:
    INSTRUMENT_4B = "REF_L"
    captured_theta = None

    def __init__(self, *_args, **kwargs):
        self.theta = kwargs["theta"]
        type(self).captured_theta = kwargs["theta"]

    def specular(self, **_kwargs):
        return np.array([0.1, 0.2, 0.3]), np.array([1.0, 2.0]), np.array([0.1, 0.2])

    def to_dict(self):
        return {}


@pytest.mark.parametrize(
    ("logs", "expected_theta"),
    [
        (
            {
                "LambdaRequest": 5.0,
                "thi": 1.5,
                "ths": 0.2,
                "BL4B:CS:Mode:Coordinates": 0,
                "BL4B:CS:ExpPl:OperatingMode": "Bound Liquid",
            },
            1.5,
        ),
        (
            {
                "LambdaRequest": 5.0,
                "thi": 1.4,
                "ths": 0.3,
                "BL4B:CS:Mode:Coordinates": 1,
                "BL4B:CS:ExpPl:OperatingMode": "Free Liquid",
            },
            0.3,
        ),
        (
            {
                "LambdaRequest": 5.0,
                "thi": 1.3,
                "ths": 0.4,
                "BL4B:CS:Mode:Coordinates": 1,
                "BL4B:CS:ExpPl:OperatingMode": "Bound Liquid",
            },
            0.4,
        ),
    ],
)
def test_process_from_template_ws_selects_theta_log(monkeypatch, logs, expected_theta):
    template_data = ReductionParameters()
    monkeypatch.setattr("lr_reduction.template.event_reduction.EventReflectivity", _DummyEventReflectivity)

    process_from_template_ws(_FakeWorkspace(logs), template_data, normalize=False)

    assert _DummyEventReflectivity.captured_theta == pytest.approx(math.radians(expected_theta))
