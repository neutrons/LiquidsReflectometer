from lr_reduction.typing import MantidWorkspace
from lr_reduction.utils import workspace_handle


class SampleLogs:
    """
    Wrapper around Mantid's run object so that `SampleLogs(workspace)[property_name]`
    returns the first value of the property if it is a vector, or the value if it is a scalar.
    Mantid's default run object `workspace.getRun()[property_name]` returns the property object, not its value.
    Usually, we're interested in the property's value, not the object itself.
    With this wrapper, we would write:
        sample_logs = SampleLogs(workspace)
        value = sample_logs[property_name]  # value if scalar, first value if vector
    instead of:
        sample_logs = workspace.getRun()
        value = sample_logs.getProperty(property_name).value  # if scalar
        value = sample_logs.getProperty(property_name).firstValue()  # if vector
    """

    def __init__(self, input_workspace: MantidWorkspace):
        self._run = workspace_handle(input_workspace).getRun()

    def __contains__(self, property_name):
        return self._run.hasProperty(property_name)

    def __getitem__(self, property_name):
        value = self._run.getProperty(property_name).value
        if isinstance(value, (int, float, str)):  # scalar sample logs can only be one of these three types
            return value
        else:
            return value[0]  # return the first value

    def property(self, property_name: str):
        """property object for the given property name"""
        return self._run.getProperty(property_name)

    def mean(self, property_name) -> float:
        """mean value of the property"""
        return self._run.getStatistics(property_name).mean
