import pytest
from mantid.simpleapi import AddSampleLog, CreateSingleValuedWorkspace, mtd

from lr_reduction.data_info import CoordinateSystem, DataType


class TestDataType:
    """Test suite for DataType enum"""

    def test_from_workspace_direct_beam_earth_centered_coordinates(self):
        """Test direct beam detection with earth-centered coordinates"""
        workspace = mtd.unique_hidden_name()
        CreateSingleValuedWorkspace(OutputWorkspace=workspace)
        AddSampleLog(Workspace=workspace, LogName="BL4B:CS:Mode:Coordinates", LogText="0", LogType="Number Series")
        AddSampleLog(Workspace=workspace, LogName="thi", LogText="0.800", LogType="Number")
        AddSampleLog(Workspace=workspace, LogName="tthd", LogText="0.805", LogType="Number")

        result = DataType.from_workspace(workspace)
        assert result == DataType.DIRECT_BEAM


    def test_from_workspace_direct_beam_free_liquid_mode(self):
        """Test direct beam detection with free liquid mode"""
        workspace = mtd.unique_hidden_name()
        CreateSingleValuedWorkspace(OutputWorkspace=workspace)
        AddSampleLog(Workspace=workspace, LogName="BL4B:CS:ExpPl:OperatingMode", LogText="Free Liquid", LogType="String")
        AddSampleLog(Workspace=workspace, LogName="thi", LogText="0.210", LogType="Number")
        AddSampleLog(Workspace=workspace, LogName="tthd", LogText="0.212", LogType="Number")

        result = DataType.from_workspace(workspace)
        assert result == DataType.DIRECT_BEAM


    def test_from_workspace_reflected_beam_earth_centered(self):
        """Test reflected beam detection with earth-centered coordinates"""
        workspace = mtd.unique_hidden_name()
        CreateSingleValuedWorkspace(OutputWorkspace=workspace)
        AddSampleLog(Workspace=workspace, LogName="BL4B:CS:Mode:Coordinates", LogText="0", LogType="Number Series")
        AddSampleLog(Workspace=workspace, LogName="thi", LogText="0.5", LogType="Number")
        AddSampleLog(Workspace=workspace, LogName="tthd", LogText="1.0", LogType="Number")

        result = DataType.from_workspace(workspace)
        assert result == DataType.REFLECTED_BEAM


    def test_from_workspace_reflected_beam_beam_centered(self):
        """Test reflected beam detection with beam-centered coordinates"""
        workspace = mtd.unique_hidden_name()
        CreateSingleValuedWorkspace(OutputWorkspace=workspace)
        AddSampleLog(Workspace=workspace, LogName="BL4B:CS:Mode:Coordinates", LogText="1", LogType="Number Series")
        AddSampleLog(Workspace=workspace, LogName="BL4B:CS:ExpPl:OperatingMode", LogText="Other", LogType="String")
        AddSampleLog(Workspace=workspace, LogName="thi", LogText="0.5", LogType="Number")
        AddSampleLog(Workspace=workspace, LogName="tthd", LogText="1.0", LogType="Number")

        result = DataType.from_workspace(workspace)
        assert result == DataType.REFLECTED_BEAM


    def test_from_workspace_direct_beam_beam_centered(self):
        """Test direct beam detection with beam-centered coordinates"""
        workspace = mtd.unique_hidden_name()
        CreateSingleValuedWorkspace(OutputWorkspace=workspace)
        AddSampleLog(Workspace=workspace, LogName="BL4B:CS:Mode:Coordinates", LogText="1", LogType="Number Series")
        AddSampleLog(Workspace=workspace, LogName="BL4B:CS:ExpPl:OperatingMode", LogText="Other", LogType="String")
        AddSampleLog(Workspace=workspace, LogName="ths", LogText="0.0006", LogType="Number")
        AddSampleLog(Workspace=workspace, LogName="tthd", LogText="0.0008", LogType="Number")

        result = DataType.from_workspace(workspace)
        assert result == DataType.DIRECT_BEAM


    def test_from_workspace_missing_logs(self):
        """Test that missing logs default to reflected beam"""
        workspace = mtd.unique_hidden_name()
        CreateSingleValuedWorkspace(OutputWorkspace=workspace)

        result = DataType.from_workspace(workspace)
        assert result == DataType.REFLECTED_BEAM


class TestCoordinateSystem:
    """Test suite for CoordinateSystem enum"""

    @pytest.mark.parametrize("log_coordinates,expected_coordinate_system", [
        ("0", CoordinateSystem.EARTH_CENTERED),
        ("1", CoordinateSystem.BEAM_CENTERED),
    ])
    def test_from_workspace_coordinates_only(self, log_coordinates, expected_coordinate_system):
        """Test coordinate system detection via coordinates log only"""
        workspace = mtd.unique_hidden_name()
        CreateSingleValuedWorkspace(OutputWorkspace=workspace)
        AddSampleLog(Workspace=workspace, LogName="BL4B:CS:Mode:Coordinates", LogText=log_coordinates, LogType="Number Series")

        result = CoordinateSystem.from_workspace(workspace)
        assert result == expected_coordinate_system

    @pytest.mark.parametrize("log_operatingmode,expected_coordinate_system", [
        ("Free Liquid", CoordinateSystem.EARTH_CENTERED),
        ("Other", CoordinateSystem.BEAM_CENTERED),
    ])
    def test_from_workspace_operatingmode_only(self, log_operatingmode, expected_coordinate_system):
        """Test coordinate system detection via operating mode log only"""
        workspace = mtd.unique_hidden_name()
        CreateSingleValuedWorkspace(OutputWorkspace=workspace)
        AddSampleLog(Workspace=workspace, LogName="BL4B:CS:ExpPl:OperatingMode", LogText=log_operatingmode, LogType="String")

        result = CoordinateSystem.from_workspace(workspace)
        assert result == expected_coordinate_system

    @pytest.mark.parametrize("log_coordinates,log_operatingmode,expected_coordinate_system", [
        ("0", "Other", CoordinateSystem.EARTH_CENTERED),
        ("0", "Free Liquid", CoordinateSystem.EARTH_CENTERED),
        ("1", "Other", CoordinateSystem.BEAM_CENTERED),
        ("1", "Free Liquid", CoordinateSystem.BEAM_CENTERED),
    ])
    def test_from_workspace_both_logs_present(self, log_coordinates, log_operatingmode, expected_coordinate_system):
        """Test coordinate system detection when both logs are present

        Note: the coordinates log takes precedence over operating mode.
        """
        workspace = mtd.unique_hidden_name()
        CreateSingleValuedWorkspace(OutputWorkspace=workspace)
        AddSampleLog(Workspace=workspace, LogName="BL4B:CS:Mode:Coordinates", LogText=log_coordinates, LogType="Number Series")
        AddSampleLog(Workspace=workspace, LogName="BL4B:CS:ExpPl:OperatingMode", LogText=log_operatingmode, LogType="String")

        result = CoordinateSystem.from_workspace(workspace)
        assert result == expected_coordinate_system
