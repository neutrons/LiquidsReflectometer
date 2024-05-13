import os
import epics
from epics import PV
 
 
# Set the EPICS_CA_ADDR_LIST environment variable
os.environ['EPICS_CA_ADDR_LIST'] = 'bl4b-dassrv1:5066'
 
name='S1'
X=PV('BL4B:Mot:s1:X:Gap.RBV')
Xcen=PV('BL4B:Mot:s1:X:Center.RBV')
Y=PV('BL4B:Mot:s1:Y:Gap.RBV')
Ycen=PV('BL4B:Mot:s1:Y:Center.RBV')
 
 
print(X.get())