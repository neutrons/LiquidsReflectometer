"""Unit tests for resolving the emission-time coefficients of a run."""

import pytest
from mantid.simpleapi import CreateSampleWorkspace, DeleteWorkspace, mtd

from lr_reduction.exceptions import CorrectionError, EmissionTimeCorrectionError, LrValidationError
from lr_reduction.models.config import EmissionTimeConfig
from lr_reduction.operations.emission_time import (
    CHOPPER_MODERATOR_LOG,
    CHOPPER_MULTIPLIER_LOG,
    CHOPPER_OFFSET_LOG,
    resolve_emission_coefficients,
)
from lr_reduction.utils.sample_logs import SampleLogs


@pytest.fixture
def make_run():
    """Factory for a workspace carrying the given chopper logs; deletes each on teardown."""
    names = []

    def make(**logs):
        name = mtd.unique_hidden_name()
        CreateSampleWorkspace(WorkspaceType="Event", NumBanks=1, BankPixelWidth=1, NumEvents=1, OutputWorkspace=name)
        for log_name, value in logs.items():
            # The DAS records chopper PVs as series; a vector log stands in for one.
            SampleLogs(name).insert(log_name, [value])
        names.append(name)
        return name

    yield make
    for name in names:
        if mtd.doesExist(name):
            DeleteWorkspace(name)


def chopper_logs(flag=1, offset=12.5, multiplier=3.25):
    return {CHOPPER_MODERATOR_LOG: flag, CHOPPER_OFFSET_LOG: offset, CHOPPER_MULTIPLIER_LOG: multiplier}


def test_disabled_correction_is_zero_even_with_logs(make_run):
    run = make_run(**chopper_logs())
    assert resolve_emission_coefficients(EmissionTimeConfig(apply=False), run) == (0.0, 0.0)


def test_disabled_correction_needs_no_logs(make_run):
    assert resolve_emission_coefficients(EmissionTimeConfig(apply=False), make_run()) == (0.0, 0.0)


def test_configured_coefficients_win_over_logs(make_run):
    run = make_run(**chopper_logs())
    config = EmissionTimeConfig(coefficients=[1.0, 2.0])
    assert resolve_emission_coefficients(config, run) == (1.0, 2.0)


def test_configured_coefficients_need_no_logs(make_run):
    config = EmissionTimeConfig(coefficients=[1.0, 2.0])
    assert resolve_emission_coefficients(config, make_run()) == (1.0, 2.0)


def test_chopper_logs_used_when_flag_set(make_run):
    run = make_run(**chopper_logs(flag=1))
    t0, t1 = resolve_emission_coefficients(EmissionTimeConfig(), run)
    assert (t0, t1) == (12.5, 3.25)
    assert isinstance(t0, float) and isinstance(t1, float)


def test_flag_clear_means_no_correction(make_run):
    run = make_run(**chopper_logs(flag=0))
    assert resolve_emission_coefficients(EmissionTimeConfig(), run) == (0.0, 0.0)


def test_missing_flag_log_raises(make_run):
    with pytest.raises(EmissionTimeCorrectionError, match=CHOPPER_MODERATOR_LOG) as raised:
        resolve_emission_coefficients(EmissionTimeConfig(), make_run())
    # The message says how to proceed.
    assert "coefficients" in str(raised.value) and "apply: false" in str(raised.value)


@pytest.mark.parametrize("missing", [CHOPPER_OFFSET_LOG, CHOPPER_MULTIPLIER_LOG])
def test_flag_set_but_coefficient_log_missing_raises(make_run, missing):
    logs = chopper_logs(flag=1)
    del logs[missing]
    with pytest.raises(EmissionTimeCorrectionError, match=missing):
        resolve_emission_coefficients(EmissionTimeConfig(), make_run(**logs))


def test_error_is_catchable_by_domain_and_kind():
    assert issubclass(EmissionTimeCorrectionError, CorrectionError)
    assert issubclass(EmissionTimeCorrectionError, LrValidationError)
