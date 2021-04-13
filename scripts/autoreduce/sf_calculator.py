"""
    Scaling factor calculation for automated reduction
"""
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np

from mantid.api import *
from mantid.kernel import *
import functools

from scipy.signal import find_peaks, peak_widths
from typing import List, Tuple


THI_TOLERANCE = 0.002


class CompareTwoNXSDataForSFcalculator(object):
    """
        will return -1, 0 or 1 according to the position of the nexusToPosition in relation to the
        nexusToCompareWith based on the following criteria
        #1: number of attenuators (ascending order)
        #2: lambda requested (descending order)
        #3: S2W (ascending order)
        #4: S2H (descending order)
        #5 if everything up to this point is identical, return 0
    """
    nexusToCompareWithRun = None
    nexusToPositionRun = None
    resultComparison = 0

    def __init__(self, nxsdataToCompareWith, nxsdataToPosition):
        self.nexusToCompareWithRun = nxsdataToCompareWith.getRun()
        self.nexusToPositionRun = nxsdataToPosition.getRun()

        compare = self.compareParameter('LambdaRequest', 'descending')
        if compare != 0:
            self.resultComparison = compare
            return

        compare = self.compareParameter('thi', 'descending', tolerance=THI_TOLERANCE)
        if compare != 0:
            self.resultComparison = compare
            return

        compare = self.compareParameter('vAtt', 'ascending')
        if compare != 0:
            self.resultComparison = compare
            return

        pcharge1 = self.nexusToCompareWithRun.getProperty('gd_prtn_chrg').value/nxsdataToCompareWith.getNEvents()
        pcharge2 = self.nexusToPositionRun.getProperty('gd_prtn_chrg').value/nxsdataToPosition.getNEvents()

        self.resultComparison = -1 if pcharge1 < pcharge2 else 1

    def compareParameter(self, param, order, tolerance=0.0):
        """
            Compare parameters for the two runs
            :param string param: name of the parameter to compare
            :param string order: ascending or descending
            :param float tolerance: tolerance to apply to the comparison [optional]
        """
        _nexusToCompareWithRun = self.nexusToCompareWithRun
        _nexusToPositionRun = self.nexusToPositionRun

        _paramNexusToCompareWith = float(_nexusToCompareWithRun.getProperty(param).value[0])
        _paramNexusToPosition = float(_nexusToPositionRun.getProperty(param).value[0])

        if abs(_paramNexusToPosition - _paramNexusToCompareWith) <= tolerance:
            return 0

        if order == 'ascending':
            resultLessThan = -1
            resultMoreThan = 1
        else:
            resultLessThan = 1
            resultMoreThan = -1

        if _paramNexusToPosition < _paramNexusToCompareWith:
            return resultLessThan
        elif _paramNexusToPosition > _paramNexusToCompareWith:
            return resultMoreThan
        else:
            return 0

    def result(self):
        return self.resultComparison


def sorter_function(r1, r2):
    """
        Sorter function used by with the 'sorted' call to sort the direct beams.
    """
    return CompareTwoNXSDataForSFcalculator(r2, r1).result()


class ScalingFactor(object):

    def __init__(self, run_list, sort_by_runs=True, sf_file='test.cfg', tof_step=200,
                 medium='Si', slit_tolerance=0.06, dry_run=False):
        self._run_list = run_list
        self._sort_by_runs = sort_by_runs
        self._sf_file = sf_file
        self._tof_steps = tof_step
        self._use_low_res_cut = False
        self._incident_medium = medium
        self._slit_tolerance = slit_tolerance
        self._dry_run = dry_run

    def execute(self):
        lr_data = []
        run_list = self._run_list

        # Load the data
        for run in run_list:
            workspace = LoadEventNexus(Filename="REF_L_%s" % run, OutputWorkspace="__data_file_%s" % run)
            lr_data.append(workspace)

        sort_by_runs = True
        if sort_by_runs is True:
            lr_data_sorted = sorted(lr_data, key=lambda r: r.getRunNumber())
        else:
            lr_data_sorted = sorted(lr_data, key=functools.cmp_to_key(sorter_function))

        # Set the output properties
        run_numbers = [r.getRunNumber() for r in lr_data_sorted]
        ws_names = [str(r) for r in lr_data_sorted]

        # Compute the scaling factors if requested
        self._compute_scaling_factors(lr_data_sorted)

    def find_peak_new(self, ws, crop=25, factor=1.) -> Tuple[List[int], List[int]]:
        """Find peak by Mantid FindPeaks

        - factor: float
            peak width multiplication factor to peak range

        Assumption
        1. The maximum count is belonged to a READ peak
        2. 
        """
        # Sum detector counts into 1D
        # x = ws.extractX()
        y = ws.extractY()
        y = np.reshape(y, (256, 304, y.shape[1]))
        p_vs_t = np.sum(y, axis=0)
        signal = np.sum(p_vs_t, axis=1)

        # Max index as the "observed" peak center
        max_index = np.argmax(signal)

        # Fit peak by Gaussian
        # create workspace
        import datetime
        now = datetime.datetime.now()
        ws_name = f'REL{now.hour:02}{now.minute:02}{now.second:02}{now.microsecond:04}.dat'
        CreateWorkspace(DataX=np.arange(len(signal)), DataY=signal, DataE=np.sqrt(signal), OutputWorkspace=ws_name)

        # prepare fitting
        model_ws_name = f'{ws_name}_model'
        param_ws_name = f'{ws_name}_parameter'
        peak_ws_name = f'{ws_name}_peaks'

        FitPeaks(InputWorkspace=ws_name,
                 OutputWorkspace=peak_ws_name,
                 PeakCenters=f'{max_index}',
                 FitWindowBoundaryList=f'{crop},{signal.shape[0]-crop}',
                 HighBackground=False,
                 ConstrainPeakPositions=False,
                 FittedPeaksWorkspace=model_ws_name,
                 OutputPeakParametersWorkspace=param_ws_name,
                 RawPeakParameters=False)

        # Retrieve value
        peak_width = mtd[param_ws_name].cell(0, 3)
        peak_center = mtd[param_ws_name].cell(0, 2)

        print(f'[INFO FIT]{ws}: Max = {max_index}, Peak center = {peak_center}, Width = {peak_width}')

        # Form output
        peak = [peak_center - factor * peak_width, peak_center + factor * peak_width]

        # Export fitting result
        model_ws = mtd[model_ws_name]
        data_ws = mtd[ws_name]
        f_name = f'{str(ws_name).split(".")[0]}.png'
        plt.figure(figsize=(10, 8))
        plt.plot(data_ws.readX(0), data_ws.readY(0), color='black', label='Counts')
        plt.plot(model_ws.readX(0), model_ws.readY(0), color='red', linestyle='-.', label='Fitted Gaussian')
        plt.plot([peak_center, peak_center], [0, np.max(signal)], color='green')
        plt.plot([peak_center - 0.5 * peak_width, peak_center + 0.5 * peak_width],
                 [0.5 * np.max(signal), 0.5 * np.max(signal)],
                 color='green')
        plt.plot([peak[0], peak[1]], [0., 0., ], color='green', linewidth=2)
        plt.xlim(120, 160)
        plt.title(f'{ws_name}: center @ {peak_center} width = {peak_width}')
        plt.legend()
        # plt.show()
        plt.savefig(f_name)
        plt.close()
        import time
        time.sleep(0.1)

        # TODO Delete workspaces
        for ws_name in [peak_ws_name, model_ws_name, param_ws_name]:
            DeleteWorkspace(ws_name)

        return peak, [0, 255]

    def find_peak(self, ws):
        """
            Find the peak in y
            TODO: find peak in x
        """
        # SaveNexusProcessed(InputWorkspace=ws, Filename=f'/tmp/{ws}.nxs')
        y=ws.extractY()
        x=ws.extractX()
        y = np.reshape(y, (256, 304, y.shape[1]))

        p_vs_t = np.sum(y, axis=0)
        counts = np.sum(p_vs_t, axis=1)

        avg = np.average(counts)

        # Crop pixels on each side where background can create a peak
        _crop=25
        peaks, props = find_peaks(counts[_crop:-_crop],
                                  threshold=None,
                                  width=3,
                                  prominence=0.5*avg)
        width = peak_widths(counts[_crop:-_crop], peaks, rel_height=0.05)

        _peak_index = 0
        _peak_max = 0
        if len(peaks)>0:
            for i in range(len(peaks)):
                if counts[peaks[i]+_crop] > _peak_max:
                    _peak_index = i
                    _peak_max = counts[peaks[i]+_crop]
                    print(f'[VZ DEBUG] _peak_max = {_peak_max} from {peaks[i]} + {_crop} =  {peaks[i]+_crop}')

        try:
            print(f'[VZ DEBUG] _peak_index = {_peak_index}, peak position = {peaks[_peak_index]}, width = {width[0][_peak_index]}')

            peak = [np.int(np.floor(peaks[_peak_index]+_crop-2.0*width[0][_peak_index])),
                    np.int(np.floor(peaks[_peak_index]+_crop+2.0*width[0][_peak_index]))]
        except:
            print(counts)
            print(avg)
            print(peaks)
            print(props)
            raise

        # FIXME DEBUG Output data set to fit
        import datetime

        content = ''
        for i in range(len(counts)):
            content += f'{i:-5}    {counts[i]:-10}\n'
        now = datetime.datetime.now()
        f_name = f'{now.hour:02}{now.minute:02}{now.second:02}{now.microsecond:04}.dat'
        with open(f_name, 'w') as dat:
            dat.write(content)

        return peak, [0, 255]
        
    def _compute_scaling_factors(self, lr_data_sorted):
        """
            If we need to compute the scaling factors, group the runs by their wavelength request
            @param lr_data_sorted: ordered list of workspaces
        """
        group_list = []
        current_group = []
        _current_wl = None
        _current_thi = None
        for r in lr_data_sorted:
            wl_ = r.getRun().getProperty('LambdaRequest').value[0]
            thi = r.getRun().getProperty('thi').value[0]

            if _current_thi is None or abs(thi-_current_thi)>THI_TOLERANCE or not _current_wl == wl_:
                # New group
                _current_wl = wl_
                _current_thi = thi
                if len(current_group)>0:
                    group_list.append(current_group)
                current_group = []

            current_group.append(r)

        # Add in the last group
        group_list.append(current_group)

        summary = ""
        for g in group_list:
            if len(g) == 0:
                continue

            direct_beam_runs = []
            peak_ranges = []
            x_ranges = []
            bck_ranges = []

            for run in g:
                print("processing: %g" % run.getRunNumber())
                peak, low_res = self.find_peak(run)
                
                peak2, _ = self.find_peak_new(run)
                print(f'[............ OUTPUT] run {run}: peak = {peak2}')

                att = run.getRun().getProperty('vAtt').value[0]-1
                wl = run.getRun().getProperty('LambdaRequest').value[0]
                thi = run.getRun().getProperty('thi').value[0]
                direct_beam_runs.append(run.getRunNumber())
                peak_ranges.append(int(peak[0]))
                peak_ranges.append(int(peak[1]))
                x_ranges.append(int(low_res[0]))
                x_ranges.append(int(low_res[1]))
                bck_ranges.append(int(peak[0])-3)
                bck_ranges.append(int(peak[1])+3)

                summary += "%10s wl=%5s thi=%5s att=%s %5s,%5s %5s,%5s ... ... (%5s, %5s)\n" % \
                    (run.getRunNumber(), wl, thi, att, peak[0], peak[1], low_res[0], low_res[1], round(peak2[0]), round(peak2[1]))

            # Determine TOF range from first file
            sample = g[0].getInstrument().getSample()
            source = g[0].getInstrument().getSource()
            source_sample_distance = sample.getDistance(source)
            detector = g[0].getDetector(0)
            sample_detector_distance = detector.getPos().getZ()
            source_detector_distance = source_sample_distance + sample_detector_distance
            h = 6.626e-34  # m^2 kg s^-1
            m = 1.675e-27  # kg
            wl = g[0].getRun().getProperty('LambdaRequest').value[0]
            chopper_speed = g[0].getRun().getProperty('SpeedRequest1').value[0]
            wl_offset = 0.0
            tof_min = source_detector_distance / h * m * (wl + wl_offset*60.0/chopper_speed - 1.7*60.0/chopper_speed) * 1e-4
            tof_max = source_detector_distance / h * m * (wl + wl_offset*60.0/chopper_speed + 1.7*60.0/chopper_speed) * 1e-4
            tof_range = [tof_min, tof_max]

            summary += "      TOF: %s\n\n" % tof_range

            # Compute the scaling factors
            if not self._dry_run:
                LRScalingFactors(DirectBeamRuns=direct_beam_runs,
                                 TOFRange=tof_range, TOFSteps=self._tof_steps,
                                 SignalPeakPixelRange=peak_ranges,
                                 SignalBackgroundPixelRange=bck_ranges,
                                 LowResolutionPixelRange=x_ranges,
                                 IncidentMedium=self._incident_medium,
                                 SlitTolerance=self._slit_tolerance,
                                 ScalingFactorFile=self._sf_file)

            print(summary)
