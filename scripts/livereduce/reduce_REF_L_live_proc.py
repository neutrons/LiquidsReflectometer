"""
    Accumulation workspace

    To test: nc -l  31466 < REF_M_25640.adara
"""
from mantid import simpleapi

simpleapi.logger.notice("Starting proc")
try:
    simpleapi.CloneWorkspace(InputWorkspace=input, OutputWorkspace=output)
except:  # noqa E722
    return input  # noqa F706
