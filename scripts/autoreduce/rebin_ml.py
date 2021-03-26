# Rebin and smooth the reflectivity curve to prepare it for ML
import sys
from mantid.simpleapi import *

def smooth_and_rebin(file_path, reference, output_file):
    ws = Load(file_path)
    ws_ref = Load(reference)
    
    ws_ref = CropWorkspace(InputWorkspace=ws_ref, XMin=0.0085, XMax=0.13)
    
    ws_smoothed = SmoothData(InputWorkspace=ws, NPoints=8)
    
    ws_ref = ConvertToHistogram(ws_ref)
    ws_smoothed = ConvertToHistogram(ws_smoothed)
    
    
    ws_final = RebinToWorkspace(WorkspaceToRebin=ws_smoothed,
                                WorkspaceToMatch=ws_ref)
    
    ws_final = ConvertToPointData(ws_final)
    
    SaveAscii(InputWorkspace=ws_final,
              Filename=output_file,
              WriteXError=False,
              WriteSpectrumID=False,
              Separator="Space")

if __name__ == "__main__":
    smooth_and_rebin(sys.argv[1], sys.argv[2], sys.argv[3])