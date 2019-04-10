# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 09:29:53 2018

MacCormack Solver for St Venant Equations. 

Simple Channel Flow

@author: Joe Shuttleworth @ Cardiff University 

Contact: shuttleworthjd@cardiff.ac.uk

"""
import numpy as np
import matplotlib.pyplot as plt

# Fixed variables
dt = 0.02
dx = 10.0
L = 500
g = 9.807
height = 0.5
Qin = 0.1
Hin = 0.5
nx = 10
k= dt/dx
TimeSim = 200
t = 0
x = list(range(nx))
# Channel and Water Elevation, Discharge and chanel width
z =np.zeros(nx)
h = np.ones(nx)*height
q = np.ones(nx)*0.01
flux=np.zeros(nx)

# Water Depth
depth = h-z
print(depth)

hpredictor = h
hcorrector = h
qpredictor = q
qcorrector = q

def MacCormack(h,q,depth,g,k,z):

    flux = fluxcalc(q, depth)
    hpredictor, qpredictor, depth = predictor(h, q, flux, g, depth, k, z)
    flux = fluxcalc(qpredictor, depth)
    hcorrector, qcorrector = corrector(h, q, flux, g, depth, hpredictor, qpredictor, k)
    return hcorrector, qcorrector, flux

def fluxcalc (q, depth):
    for i in range(0,nx):
        flux[i] = q[i]**2/depth[i]
    return flux

def predictor(h,q,flux,g,depth,k, z):
    for i in range(1,nx-1):
        hpredictor[i] = h[i] - k * (q[i+1] - q[i])
        qpredictor[i] = q[i]- k*(flux[i+1]-flux[i])-0.5*k*g*(depth[i]+depth[i+1])*(h[i+1]-h[i])
        depth[i] = hpredictor[i] - z[i]
    return hpredictor, qpredictor, depth

def corrector(h, q, flux, g, depth, hpredictor, qpredictor, k):
    for i in range(1,nx-1):
        hcorrector[i] = h[i] - k* (qpredictor[i] - qpredictor[i-1])
        qcorrector[i] = q[i] - k* (flux[i] - flux[i-1]) - 0.5*k*g * (depth[i]+depth[i-1]) * (h[i]-h[i-1])
        # Averaging Step - MacCormack Step        
        hcorrector[i] = 0.5*(hpredictor[i] + hcorrector[i])
        qcorrector[i] = 0.5*(qpredictor[i] + qcorrector[i])
    return hcorrector, qcorrector

def resetvalues(hcorrector,qcorrector,h,q):
    h = hcorrector
    q = qcorrector
    depth = h - z
    return h, q, depth

def bound(q,h,depth,qpredictor,hpredictor):
    q[0]=Qin
    h[0]=Hin
    depth[0]=0.5
    qpredictor[0] = Qin
    hpredictor[0] = Hin
    return q, h, depth, qpredictor, hpredictor
nt = 0

while t < TimeSim:
    hcorrector, qcorrector, flux = MacCormack(h,q,depth,g,k,z)
    h,q,depth = resetvalues(hcorrector, qcorrector, h, q)
    t = t + dt
    nt = nt + 1
    q,h,depth,qpredictor, hpredictor = bound(q,h,depth, qpredictor, hpredictor)
    
    
plt.plot(x,h)
plt.plot(x,q)
print(nt)
    