.. _dead_time_correction:

SingleReadoutDeadTimeCorrection
===============================

Dead time is the time after an event that a detector is not able to detect another event.
For a paralyzable detector, an event that happens during the dead time restarts the dead time. For
a non-paralyzable detector, the event is simply lost and does not cause additional dead time.

Dead-time correction corrects for detector dead time by weighing the events according to:

.. math:: N = M \frac{1}{(1-\mathrm{rate} \cdot (\frac{t_{\mathrm{dead}}}{t_{\mathrm{bin}}}))}

for non-paralyzable detectors and

.. math:: N = M \frac{\mathrm{Re} (\mathrm{W}(-\mathrm{rate} \cdot (\frac{t_{\mathrm{dead}}}{t_{\mathrm{bin}}})) )}{\frac{\mathrm{rate}}{t_{\mathrm{bin}}}}

for paralyzable detectors, where

| :math:`N` = true count
| :math:`M` = measured count
| :math:`t_{\mathrm{dead}}` = dead time
| :math:`t_{\mathrm{bin}}` = TOF bin width
| :math:`\mathrm{rate}` = measured count rate
| :math:`\mathrm{W}` = `Lambert W function <https://en.wikipedia.org/wiki/Lambert_W_function>`_

The class ``SingleReadoutDeadTimeCorrection`` is a Mantid-style algorithm for computing the
dead-time correction for an event workspace. One can optionally include error events in the
dead-time computation.

Properties
----------

.. list-table::
   :widths: 20 20 20 20 20
   :header-rows: 1

   * - Name
     - Direction
     - Type
     - Default
     - Description
   * - InputWorkspace
     - Input
     - EventWorkspace
     - Mandatory
     - Input workspace used to compute dead-time correction
   * - InputErrorEventsWorkspace
     - Input
     - EventWorkspace
     -
     - Input workspace with error events used to compute dead-time correction
   * - DeadTime
     - Input
     - number
     - 4.2
     - Dead time in microseconds
   * - UseDeadTimeThreshold
     - Input
     - boolean
     - False
     - If True, use a correction of 0 for TOF bins requiring corrections greater than ``DeadTimeThreshold``
   * - DeadTimeThreshold
     - Input, Optional
     - number
     - 1.5
     - If ``UseDeadTimeThreshold`` is True, this is the upper limit for dead-time correction ratios
   * - TOFStep
     - Input
     - number
     - 100.0
     - TOF bins to compute dead-time correction, in microseconds
   * - Paralyzable
     - Input
     - boolean
     - False
     - If True, paralyzable correction will be applied, non-paralyzable otherwise
   * - TOFRange
     - Input
     - dbl list
     - [0.0, 0.0]
     - TOF range to use to compute dead-time correction
   * - OutputWorkspace
     - Output
     - MatrixWorkspace
     - Mandatory
     - Output workspace containing the dead-time correction factor for each TOF bin

Usage
-----
Example using ``SingleReadoutDeadTimeCorrection``

.. testcode::

    import mantid.simpleapi as mtd_api
    from lr_reduction import template
    from lr_reduction.DeadTimeCorrection import SingleReadoutDeadTimeCorrection
    from lr_reduction.utils import amend_config

    mtd_api.config["default.facility"] = "SNS"
    mtd_api.config["default.instrument"] = "REF_L"


    def test_deadtime(nexus_dir):
    """
    Test the time-resolved reduction that uses a measured reference.
    It is generally used at 30 Hz but it also works at 60 Hz.
    """
    with amend_config(data_dir=nexus_dir):
        ws = mtd_api.Load("REF_L_198409")

    algo = SingleReadoutDeadTimeCorrection()
    algo.PyInit()
    algo.setProperty("InputWorkspace", ws)
    algo.setProperty("OutputWorkspace", "dead_time_corr")
    algo.setProperty("UseDeadTimeThreshold", True)
    algo.setProperty("DeadTimeThreshold", 1.1)

    algo.PyExec()
    corr_ws = algo.getProperty("OutputWorkspace").value
    corr = corr_ws.readY(0)
    for c in corr:
        assert c <= 1.1
        
