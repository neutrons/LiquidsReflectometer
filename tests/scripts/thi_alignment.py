"""
Created on Mon Oct  3 09:36:53 2022
 
@author: esw
"""
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
 
 
#-------------------------------------------------
 
def get_ypix_v_i(fname,nx,ny):
    f = h5py.File(fname, 'r')
   
 
    eid=np.array(f['entry/bank1_events/event_id'])               
    ect=np.array(f['entry/bank1_events/total_counts'])
    ect=ect[0]   
 
    ipix=np.zeros(ny)
 
    #convert the event id into and x and y pixel value
    #x= eid // ny
    y= eid % ny
   
    for i in range(ect):
        ipix[y[i]]=ipix[y[i]]+1
   
    return ipix
#-----------------------------------------------
 
 
def get_zd(fname):
    f = h5py.File(fname, 'r')
    return np.array(f['entry/bank1_events/event_id'])               
 
#-----------------------------------------------
 
#-----------------------------------------------
 
def gauss(x, a, b, c):
    return a*np.exp(-(x-b)**2/(2*c**2))            
 
#-----------------------------------------------
 
nx=256
ny=304
 
path='C:/Users/esw/Documents/NR/4B/INSTRUMENT/water/'
 
#zi scans (0.70)
runs=[197310,197327]
runs=[204912,204921]
 
zstep=1
 
th=[-0.2,-0.4,-0.6,-0.8,-1.0]
 
 
runs=[]
th=[-0.25,-0.45,-0.65,-0.85,-1.05]
 
#-------------------------------
 
nn=runs[1]-runs[0]
ypix=np.arange(ny)
 
plt.cla()
 
mm=np.zeros(nn+1)
pixCOM=np.zeros(nn+1)
pixGauss=np.zeros(nn+1)
ypix=np.arange(ny)
 
for i in range(nn+1):
    name='REF_L_'+str(runs[0]+i)+'.nxs.h5'
    #print(name)
 
    ipix=get_ypix_v_i(path+name,nx,ny)
    plt.plot(ypix,ipix)
    cen=0
    mm[i]=i*zstep
 
    plt.plot(ipix)
    plt.xlim(0,170)
   
    #get the center of mass of the distribution
    com=sum(ypix*ipix)/sum(ipix)
    #fit a gaussian to the distribution
    par,cov = curve_fit(gauss,ypix,ipix, p0=(max(ipix),com,1))   
 
    yy=np.arange(min(ypix),max(ypix),0.01)
    fit = gauss(yy, par[0],par[1],par[2])   
    
    pixCOM[i]=com
    pixGauss[i]=par[1]
 
d1=np.zeros(int(len(pixCOM)/2))
d2=np.zeros(int(len(pixCOM)/2))
 
for i in range(0,len(pixCOM),2):
    v=int(i/2)
    print(i,v)
    print((pixCOM[i]-pixCOM[i+1]),(pixGauss[i]-pixGauss[i+1]))
    d1[v]=(pixGauss[i]-pixGauss[i+1])
    d2[v]=(pixCOM[i]-pixCOM[i+1])
 
plt.cla()
plt.plot(d1,th,'o')
plt.plot(d2,th,'o')  
 
#now fit a line
a, b = np.polyfit(d1,th, 1)
print(a,b)
x=np.arange(0,100.0,0.01)
y=a*x+b
plt.plot(x,y)
 
#now fit a line
a, b = np.polyfit(d2,th, 1)
print(a,b)
x=np.arange(0,100.0,0.01)
y=a*x+b
plt.plot(x,y)
 
