import subprocess
import sys
from unittest.mock import patch

import pytest

from lr_autoreduce.reduce_REF_L import confirm_data_availability, parse_command_arguments, str_to_bool


def test_parse_args_minimal():
    argv = ["prog", "events.nxs.h5", "output"]
    with patch.object(sys, "argv", argv):
        args = parse_command_arguments()

    assert args.events_file == "events.nxs.h5"
    assert args.output_dir == "output"

    # Optional positionals should fall back to defaults
    assert args.old_version_flag is None
    assert args.template_file is None
    assert args.avg_overlap == "false"
    assert args.const_q == "false"
    assert args.fit_first_peak == "false"
    assert args.theta_offset == "0"

    # Optional flag
    assert args.no_publish is False


def test_parse_args_all_positionals():
    argv = [
        "prog",
        "events.nxs.h5",
        "output",
        "1",  # old_version_flag
        "template.xml",  # template_file
        "true",  # avg_overlap
        "true",  # const_q
        "true",  # fit_first_peak
        "5.0",  # theta_offset
    ]
    with patch.object(sys, "argv", argv):
        args = parse_command_arguments()

    assert args.events_file == "events.nxs.h5"
    assert args.output_dir == "output"
    assert args.old_version_flag == "1"
    assert args.template_file == "template.xml"
    assert args.avg_overlap == "true"
    assert args.const_q == "true"
    assert args.fit_first_peak == "true"
    assert args.theta_offset == "5.0"


def test_parse_args_no_publish():
    argv = [
        "prog",
        "events.nxs.h5",
        "output",
        "--no_publish",
    ]
    with patch.object(sys, "argv", argv):
        args = parse_command_arguments()

    assert args.no_publish is True


def test_str_to_bool():
    assert str_to_bool("true") is True
    assert str_to_bool("True") is True
    assert str_to_bool("false") is False
    assert str_to_bool("FALSE") is False


@pytest.fixture
def sample_logs():
    # Minimal stand-in for SampleLogValues
    return {"experiment_identifier": "IPTS-123456"}


def test_confirm_data_exception_logs_notice(sample_logs):
    with patch("lr_autoreduce.reduce_REF_L.subprocess.run", side_effect=subprocess.CalledProcessError(1, "cmd")):
        with patch("lr_autoreduce.reduce_REF_L.logger") as mock_logger:
            confirm_data_availability(sample_logs)

    mock_logger.notice.assert_called_once()
    assert "Could not set data availability" in mock_logger.notice.call_args[0][0]
