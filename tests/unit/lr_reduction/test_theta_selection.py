import pytest

from lr_reduction.theta_selection import is_earth_centered_geometry, theta_log_name


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
    assert is_earth_centered_geometry(logs) is expected


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
    assert theta_log_name(logs) == expected
