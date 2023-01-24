#!/usr/bin/python3
import sys
import os
import time
import threading

import matplotlib.pyplot as plt

if ("MANTIDPATH" in os.environ):
    del os.environ["MANTIDPATH"]

sys.path.insert(0,"/opt/mantid63/bin")
sys.path.insert(1,"/opt/mantid63/lib")

import mantid
mantid.kernel.config.setLogLevel(3)
from mantid.simpleapi import *
import numpy as np



print(mantid.__version__)
LIVE_DATA_WS = 'live_data'
UPDATE_TIME = 30 # seconds



def thread_function():
    is_alive = True
    #fig, ax = plt.subplots(figsize=(6,6))
    plt.subplots(2, 1, dpi=100, figsize=(6,9), sharex=True)
    plt.subplots_adjust(hspace=0.5)
    ax = plt.subplot(2, 1, 1)
    ax2 = plt.subplot(2, 1, 2)
    previous_data = None
    previous_delta = None

    while is_alive:
        if LIVE_DATA_WS in mtd:
            try:
                n_events = mtd[LIVE_DATA_WS].getNumberEvents()
                #ipts = mtd[LIVE_DATA_WS].getRun()["experiment_identifier"].value
                #run_number = mtd[LIVE_DATA_WS].getRun()["run_number"].value
                print("Events: %g" % n_events)

                if n_events == 0:
                    continue

                ws = SumSpectra(LIVE_DATA_WS)
                tof = Rebin(ws, [ws.getTofMin(), 300, ws.getTofMax()], OutputWorkspace='tof_')
                x = tof.readX(0)
                x = (x[1:]+x[:-1])/2.0
                y = tof.readY(0)
                
                ax.clear()
                total = np.sum(y)
                ax.step(x, y/total, where='mid', label='Total')
                
                if previous_data:
                    tof_previous_data = Rebin(previous_data, [ws.getTofMin(), 300, ws.getTofMax()],
                                              OutputWorkspace='tof_previous_data_')

                    x_prev = tof_previous_data.readX(0)
                    x_prev = (x_prev[1:]+x_prev[:-1])/2.0
                    y_prev = tof_previous_data.readY(0)
               
                    total = np.sum(y_prev)
                    ax.step(x_prev, y_prev/total, where='mid')

                ax.set_title('%g events | %s' % (n_events, time.ctime()))
                ax.set_xlabel('TOF')
                ax.set_ylabel('Events')
                ax.legend(['Total', 'Previous'])
                #ax.set_yscale('log')
                #ax.set_xscale('log')
                
                # Difference between this chunk and the previous
                
                if previous_data:
                    n_previous_events = previous_data.getNumberEvents()
                    if n_events > n_previous_events:              
                        ax2.clear()
                        y_delta = y-y_prev
                        total_d = np.sum(y_delta)
                        #ax2.step(x_prev, y_prev/total, where='mid', label='Total')
                        ax2.step(x_prev, y_delta/total_d, where='mid', label='Difference')
                                            
                        if previous_delta is not None:
                            ax2.step(x_prev, previous_delta, where='mid', label='Previous Diff')
                        previous_delta = y_delta/total_d
                        
                        ax2.set_title('Difference from previous [%g -> %g events]' % (n_previous_events, n_events))
                        ax2.set_xlabel('TOF')
                        ax2.set_ylabel('Events')
                        ax2.legend(['Difference', 'Previous diff'])
                else:
                    ax2.clear()
                    ax2.set_title('No data yet')

                plt.pause(10)
                #plt.savefig('/SNS/REF_L/shared/livedata.png')
        
                DeleteWorkspace('tof_')
                previous_data = CloneWorkspace(ws)
            except:
                print(sys.exc_info()[1])
                print('Stopping')
                AlgorithmManager.cancelAll()
                is_alive = False



        time.sleep(int(UPDATE_TIME/2.0))


x = threading.Thread(target=thread_function)
x.start()

try:
    StartLiveData(Instrument='REF_L',
                  FromNow=False,
                  FromStartOfRun=True,
                  UpdateEvery=UPDATE_TIME,
                  Listener='SNSLiveEventDataListener',
                  Address='bl4b-daq1.sns.gov:31415',
                  #PostProcessingScript='output=input',
                  AccumulationWorkspace='acc',
                  PreserveEvents=True,
                  OutputWorkspace=LIVE_DATA_WS)

    # If we were to have only the StartLiveData call, we'd have to 
    # keep this process alive:
    #time.sleep(2000)
except:
    print("STOPPING")
    AlgorithmManager.cancelAll()
    time.sleep(1)
        
        
        
                  
