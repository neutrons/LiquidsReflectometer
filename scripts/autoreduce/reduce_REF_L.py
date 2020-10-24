"""
    Auto-reduction script for the Liquids Reflectometer
    For reference:
        Type 0: Normal sample data
        Type 1: Direct beams for scaling factors
        Type 2: Zero-attenuator direct beams
        Type 3: Data that we don't need to treat
"""
import sys
import os
import json
import warnings
warnings.simplefilter('ignore', RuntimeWarning)


if ("MANTIDPATH" in os.environ):
    del os.environ["MANTIDPATH"]
#sys.path.insert(0,"/opt/mantid50/bin")
#sys.path.insert(1,"/opt/mantid50/lib")
sys.path.insert(0,"/opt/mantidnightly/bin")
sys.path.insert(1,"/opt/mantidnightly/lib")

import mantid
from mantid.simpleapi import *

try:
    from sf_calculator import ScalingFactor
    DIRECT_BEAM_CALC_AVAILABLE = True
except:
    import scipy
    logger.notice("Scaling factor calculation upgrade not available: scipy=%s" % scipy.__version__)
    DIRECT_BEAM_CALC_AVAILABLE = False

event_file_path=sys.argv[1]
output_dir=sys.argv[2]

event_file = os.path.split(event_file_path)[-1]
# The legacy format is REF_L_xyz_event.nxs
# The new format is REF_L_xyz.nxs.h5
run_number = event_file.split('_')[2]
run_number = run_number.replace('.nxs.h5', '')


# Reduction options
#-------------------------------------------------------------------------
# Wavelength below which we don't need the absolute normalization
WL_CUTOFF = 10.0  

# Default primary fraction range to be used if it is not defined in the template
PRIMARY_FRACTION_RANGE = [5, 290]

NORMALIZE_TO_UNITY = False

# Allowed values: ['DirectBeam', 'WithReference']
NORMALIZATION_TYPE = "WithReference"
#NORMALIZATION_TYPE = "DirectBeam"

# Allowed values: dict or ""
# D2O
REFL1D_PARS = json.dumps(dict(back_sld=6.4,
                              back_roughness=2.7,
                              front_sld=0,
                              layers=[],
                              scale=1.0,
                              background=0.0))
# Quartz
#REFL1D_PARS = json.dumps(dict(back_sld=4.09,
#                              back_roughness=4.28,
#                              front_sld=0,
#                              layers=[],
#                              scale=0.9169,
#                              background=3.753e-07))
#-------------------------------------------------------------------------


# Locate the template file
# If no template file is available, the automated reduction will generate one
template_file = ""
if os.path.isfile("template.xml"):
    template_file = "template.xml"
elif os.path.isfile(os.path.join(output_dir, "template.xml")):
    template_file = os.path.join(output_dir, "template.xml")
elif os.path.isfile("/SNS/REF_L/shared/autoreduce/template.xml"):
    template_file = "/SNS/REF_L/shared/autoreduce/template.xml"

print("Using template: %s" % template_file)
# Run the auto-reduction
ws = LoadEventNexus(Filename=event_file_path)

# Check the measurement geometry
if ws.getRun().getProperty('BL4B:CS:ExpPl:OperatingMode').value[0] == 'Free Liquid':
    NORMALIZATION_TYPE = "WithReference"
else:
    NORMALIZATION_TYPE = "DirectBeam"

# Determine whether this is data or whether we need to compute scaling factors
data_type = ws.getRun().getProperty("data_type").value[0]

if data_type == 1 and DIRECT_BEAM_CALC_AVAILABLE:
    sequence_number = ws.getRun().getProperty("sequence_number").value[0]
    first_run_of_set = ws.getRun().getProperty("sequence_id").value[0]
    incident_medium = ws.getRun().getProperty("incident_medium").value[0]

    _fpath = os.path.join(output_dir, "sf_%s_%s_auto.cfg" % (first_run_of_set, incident_medium))

    sf = ScalingFactor(run_list=range(first_run_of_set, first_run_of_set + sequence_number),
                       sort_by_runs=True,
                       tof_step=200,
                       medium=incident_medium,
                       slit_tolerance=0.06,
                       sf_file=_fpath)
else:
    output = LRAutoReduction(#Filename=event_file_path,
                             InputWorkspace=ws,
                             ScaleToUnity=NORMALIZE_TO_UNITY,
                             ScalingWavelengthCutoff=WL_CUTOFF,
                             PrimaryFractionRange=PRIMARY_FRACTION_RANGE,
                             OutputDirectory=output_dir,
                             SlitTolerance=0.06,
                             ReadSequenceFromFile=True,
                             OrderDirectBeamsByRunNumber=True,
                             TemplateFile=template_file, FindPeaks=False,
                             NormalizationType=NORMALIZATION_TYPE,
                             Refl1DModelParameters=REFL1D_PARS)
    first_run_of_set=int(output[1])


#-------------------------------------------------------------------------
# Produce plot for the web monitor
default_file_name = 'REFL_%s_combined_data_auto.txt' % first_run_of_set
if os.path.isfile(default_file_name):
    print("Loading %s" % os.path.join(output_dir, default_file_name))
    reflectivity = LoadAscii(Filename=os.path.join(output_dir, default_file_name), Unit="MomentumTransfer")

    try:
        from postprocessing.publish_plot import plot1d
    except ImportError:
        from finddata.publish_plot import plot1d
    x = reflectivity.readX(0)
    y = reflectivity.readY(0)
    dy = reflectivity.readE(0)
    dx = reflectivity.readDx(0)
    
    if int(run_number) - first_run_of_set < 10:
        for r in range(0, 10):
            if os.path.isfile('REFL_%s_%s_%s_auto.nxs' % (first_run_of_set, r+1, first_run_of_set+r)):
                plot1d(first_run_of_set+r, [[x, y, dy, dx]], instrument='REF_L', 
                       x_title=u"Q (1/A)", x_log=True,
                       y_title="Reflectivity", y_log=True, show_dx=False)
    else:
        plot1d(run_number, [[x, y, dy, dx]], instrument='REF_L', 
               x_title=u"Q (1/A)", x_log=True,
               y_title="Reflectivity", y_log=True, show_dx=False)


