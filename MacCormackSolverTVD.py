# -*- coding: utf-8 -*-
"""
Created on Fri Apr 13 09:29:53 2018

MacCormack Solver for St Venant Equations. 

Simple Channel Flow - With Shock Capturing Capabilities

@author: Joe Shuttleworth @ Cardiff University 

Contact: shuttleworthjd@cardiff.ac.uk

"""
import numpy as np
import matplotlib.pyplot as plt

# Fixed variables
dt = 1
dx = 10.0
L = 500
g = 9.807
height = 0.5
Qin = 0.1
Hin = 0.5
nx = 50
k= dt/dx
TimeSim = 200
TVDChoice = 1
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

TVDcont = np.zeros(nx)
TVDmom = np.zeros(nx)
cr = np.zeros(nx)
Grplus = np.zeros(nx)
Grminus = np.zeros(nx)
hnew = np.ones(nx)*Hin
qnew = np.ones(nx)*Qin


def MacCormack(h,q,depth,g,k,z):

    flux = fluxcalc(q, depth)
    hpredictor, qpredictor, depth = predictor(h, q, flux, g, depth, k, z)
    flux = fluxcalc(qpredictor, depth)
    hcorrector, qcorrector = corrector(h, q, flux, g, depth, hpredictor, qpredictor, k)
    hnew, qnew = calc(hcorrector, qcorrector, hpredictor, qpredictor)
    if TVDChoice == 1:
        TVDcont, TVDmom = TVD(hnew,qnew,depth,dt,dx)
        for i in range(1,nx-1):
            hnew[i] = hnew[i] + TVDcont[i]
            qnew[i] = qnew[i] + TVDmom[i]    
    return hnew, qnew, flux

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
    return hcorrector, qcorrector

def calc(hcorrector, qcorrector, hpredictor, qpredictor):
    for i in range(1, nx-1):
        hnew[i] = 0.5*(hpredictor[i] + hcorrector[i])
        qnew[i] = 0.5*(qpredictor[i] + qcorrector[i])
        depth[i] = hnew[i] - z[i]
    return hnew, qnew

def resetvalues(hnew,qnew,h,q):
    h = hnew
    q = qnew
    depth = h - z
    return h, q, depth

def bound(q,h,depth,qpredictor,hpredictor):
    q[0]=Qin
    h[0]=Hin
    depth[0]=0.5
    qpredictor[0] = Qin
    hpredictor[0] = Hin
    return q, h, depth, qpredictor, hpredictor

def TVD(h,q,depth,dt,dx):
    for i in range(1,nx-1):
        hminus = h[i] - h[i-1]
        hplus = h[i+1] - h[i]
        qminus = q[i] - q[i-1]
        qplus = q[i+1] - q[i]
        rplus = (hminus*hplus + qminus*qplus)/(hplus*hplus + qplus*qplus)
        rminus = (hminus*hplus + qminus*qplus)/(hminus*hminus +qminus*qminus)
        sigmaplus = max(0,min(2*rplus,1))
        sigmaminus = max(0,min(2*rminus,1))
        cr[i] = (abs(q[i]/depth[i]) + np.sqrt(g*depth[i]*dt/dx))
        if cr[1] <= 0.5:
            c = cr[i] * (1-cr[i])
        else:
            c = 0.25
        Grplus[i] = 0.5*c*(1-sigmaplus)
        Grminus[i] = 0.5*c*(1-sigmaminus)
    for i in range(1,nx-1):
        TVDcont[i] = (Grplus[i] + Grminus[i+1])*(h[i+1]-h[i]) - (Grplus[i-1]+Grminus[i])*(h[i]-h[i+1])
        TVDmom[i] = (Grplus[i] + Grminus[i+1])*(q[i+1]-q[i]) - (Grplus[i-1]+Grminus[i])*(q[i]-q[i+1])
    return TVDcont, TVDmom

nt = 0

while t < TimeSim:
    hcorrector, qcorrector, flux = MacCormack(h,q,depth,g,k,z)
    h,q,depth = resetvalues(hnew, qnew, h, q)
    t = t + dt
    nt = nt + 1
    q,h,depth,qpredictor, hpredictor = bound(q,h,depth, qpredictor, hpredictor)
    
    
plt.plot(x,h)
plt.plot(x,q)
print(nt)
    