import pytest
from mantid.simpleapi import AddSampleLog, CreateWorkspace, DeleteWorkspace, mtd

from lr_reduction.mantid_utils import SampleLogValues


class TestSampleLogValues:
    @classmethod
    def setup_class(cls):
        """
        Setup method runs once for the whole test suite.

        This method creates a workspace with sample data and adds several sample logs
        of different types (integer, float, string, and number series) to the workspace.
        """
        ws = CreateWorkspace(DataX=[0, 1], DataY=[0, 10], OutputWorkspace=mtd.unique_hidden_name())
        AddSampleLog(ws, LogName="property_int", LogText="42", LogType="Number", NumberType="Int")
        AddSampleLog(ws, LogName="property_float", LogText="3.14", LogType="Number", NumberType="Double")
        AddSampleLog(ws, LogName="property_str", LogText="hello", LogType="String")
        AddSampleLog(ws, LogName="property_series", LogText="42", LogType="Number Series")
        cls.workspace = str(ws)

    @classmethod
    def teardown_class(cls):
        """Delete the workspace after all tests in the suite have run."""
        DeleteWorkspace(cls.workspace)

    def test_initialization(self):
        assert SampleLogValues(self.workspace)

    def test_contains_property(self):
        sample_logs = SampleLogValues(self.workspace)
        for property_name in ["property_int", "property_float", "property_str", "property_series"]:
            assert property_name in sample_logs
        assert sample_logs["property_int"] == 42
        assert sample_logs["property_float"] == 3.14
        assert sample_logs["property_str"] == "hello"
        assert sample_logs["property_series"] == 42

    def test_does_not_contain_property(self):
        sample_logs = SampleLogValues(self.workspace)
        assert "property_unknown" not in sample_logs

    def test_property(self):
        sample_logs = SampleLogValues(self.workspace)
        assert sample_logs.property("property_int").value == 42
        assert sample_logs.property("property_float").value == 3.14
        assert sample_logs.property("property_str").value == "hello"
        assert sample_logs.property("property_series").firstValue() == 42

    def test_mean(self):
        sample_logs = SampleLogValues(self.workspace)
        assert sample_logs.mean("property_series") == 42


if __name__ == "__main__":
    pytest.main([__file__])
