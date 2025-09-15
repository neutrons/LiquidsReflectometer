import numpy as np
import pytest
from mantid.simpleapi import AddSampleLog, CreateSingleValuedWorkspace, LoadEmptyInstrument, mtd
from numpy.testing import assert_allclose, assert_almost_equal

from lr_reduction.gravity_correction import GravityDirection, _log_value, _theta_in, _theta_sample, gravity_correction
from lr_reduction.utils import workspace_handle


class TestGravityDirection:
    """Test suite for GravityDirection enum"""

    def test_enum_values(self):
        """Test that enum values are correctly defined"""
        assert GravityDirection.DOWN == -1
        assert GravityDirection.OFF == 0
        assert GravityDirection.UP == 1

    def test_enum_creation(self):
        """Test creating GravityDirection from integer values"""
        assert GravityDirection(-1) == GravityDirection.DOWN
        assert GravityDirection(0) == GravityDirection.OFF
        assert GravityDirection(1) == GravityDirection.UP
        with pytest.raises(ValueError):
            GravityDirection(2)

    def test_from_value_with(self):
        """Test from_value class method with valid integer inputs"""
        assert GravityDirection.from_value(-1) == GravityDirection.DOWN
        assert GravityDirection.from_value(0) == GravityDirection.OFF
        assert GravityDirection.from_value(1) == GravityDirection.UP
        assert GravityDirection.from_value(None) is None
        with pytest.raises(ValueError):
            GravityDirection.from_value(2)

    def test_int_conversion(self):
        """Test converting enum values back to integers"""
        assert int(GravityDirection.DOWN) == -1
        assert int(GravityDirection.OFF) == 0
        assert int(GravityDirection.UP) == 1

    def test_string_representation(self):
        """Test string representation of enum values"""
        assert str(GravityDirection.DOWN) == "-1"
        assert str(GravityDirection.OFF) == "0"
        assert str(GravityDirection.UP) == "1"


def test_log_value():
    workspace = mtd.unique_hidden_name()
    CreateSingleValuedWorkspace(OutputWorkspace=workspace)
    run = workspace_handle(workspace).getRun()


    with pytest.raises(RuntimeError):
        _log_value(run, "NonExistent Log")
    assert _log_value(run, "NonExistent Log", default=100) == 100

    AddSampleLog(Workspace=workspace, LogName="test_log", LogText="42.0", LogType="Number")
    assert _log_value(run, "test_log") == 42.0
    assert _log_value(run, "test_log", default=0) == 42.0

    AddSampleLog(Workspace=workspace, LogName="test_log", LogText="42.0", LogType="Number Series")
    assert _log_value(run, "test_log") == 42.0
    assert _log_value(run, "test_log", default=0) == 42.0

    AddSampleLog(Workspace=workspace, LogName="test_log", LogText="42.0 43.0", LogType="Number Series")
    assert _log_value(run, "test_log") == 42.0
    assert _log_value(run, "test_log", default=0) == 42.0

    AddSampleLog(Workspace=workspace, LogName="test_log", LogText="one string", LogType="String")
    assert _log_value(run, "test_log") == "one string"


def test_theta_in():
    workspace = mtd.unique_hidden_name()
    CreateSingleValuedWorkspace(OutputWorkspace=workspace)
    AddSampleLog(Workspace=workspace, LogName="BL4B:Mot:thi.RBV", LogText="42.0", LogType="Number Series")

    # test last `else`
    assert _theta_in(workspace) == 38.0

    # test `elif`
    AddSampleLog(Workspace=workspace, LogName="BL4B:CS:ExpPl:OperatingMode", LogText="ON", LogType="String")
    AddSampleLog(Workspace=workspace, LogName="BL4B:CS:ExpPl:OperatingMode", LogText="Bound Liquid", LogType="String")
    assert _theta_in(workspace) == 38.0
    AddSampleLog(Workspace=workspace, LogName="BL4B:CS:ExpPl:OperatingMode", LogText="Free Liquid", LogType="String")
    assert _theta_in(workspace) == 42.0

    # test first `if`
    AddSampleLog(Workspace=workspace, LogName="BL4B:CS:Mode:Coordinates", LogText="0", LogType="Number Series")
    assert _theta_in(workspace) == 42.0
    AddSampleLog(Workspace=workspace, LogName="BL4B:CS:Mode:Coordinates", LogText="1", LogType="Number Series")
    assert _theta_in(workspace) == 38.0


def test_theta_sample():
    """Test _theta_sample function with various workspace configurations"""
    workspace = mtd.unique_hidden_name()
    CreateSingleValuedWorkspace(OutputWorkspace=workspace)

    # Test with default instrument values
    angles = _theta_sample(workspace, wavelengths=np.array([3.0, 6.0, 9.0]), theta_in=2.0)
    assert_almost_equal(angles, [2.00026, 2.00104, 2.00235], decimal=4)

def test_gravity_correction():
    workspace = mtd.unique_hidden_name()
    LoadEmptyInstrument(InstrumentName="REF_L", OutputWorkspace=workspace)

    # Test with default instrument values
    # no correction (no "BL4B:Mot:ths.RBV" log value)
    angles = gravity_correction(workspace, wavelengths=np.array([3.0, 6.0, 9.0]), theta_in=2.0)
    assert_almost_equal(angles, [0.00000, 0.00000, 0.00000], decimal=4)

    AddSampleLog(Workspace=workspace, LogName="BL4B:Mot:ths.RBV", LogText="1.0", LogType="Number")
    angles = gravity_correction(workspace, wavelengths=np.array([3.0, 6.0, 9.0]), theta_in=2.0)
    assert_allclose(angles, [4.5e-06, 18.2e-06, 41.0e-06], atol=1e-06)

    AddSampleLog(Workspace=workspace, LogName="BL4B:Mot:ths.RBV", LogText="-1.0", LogType="Number")
    angles = gravity_correction(workspace, wavelengths=np.array([3.0, 6.0, 9.0]), theta_in=2.0)
    assert_allclose(angles, [-4.5e-06, -18.2e-06, -41.0e-06], atol=1e-06)
