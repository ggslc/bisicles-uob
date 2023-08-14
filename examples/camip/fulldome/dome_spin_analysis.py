#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 10:08:48 2023

@author: ggslc
"""

import sys
import os
import glob
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.getcwd() + '/../python')
from bisiclesIO import BisiclesData

def km(xm):
    return xm * 1.0e-3

def circ(r):
    th = np.linspace(0,2.*np.pi,100)
    x = r * np.cos(th)
    y = r * np.sin(th)
    return x,y    

def plot_state(bike, ax):


    pc = ax.pcolormesh(km(bike.x), km(bike.y), np.ma.masked_array(bike.thk, bike.thk < 1e-5), 
                       vmin=0, vmax=1600, cmap='viridis_r')
    ax.contour(km(bike.x), km(bike.y), bike.ice_frac, [0.1, 0.5, 0.9], colors=['r','k','r'], 
               linewidths=[0.5,0.5,0.5])
    #plt.colorbar(pc, label='h',shrink=0.5,orientation='horizontal')
    ax.contour(km(bike.x), km(bike.y), bike.hab,[0])
    
    w = int(len(bike.x)/16)
    xs = slice(int(w/2), len(bike.x), w)
    ys = slice(int(w/2), len(bike.y), w)

    
    ax.quiver(km(bike.x[xs]), km(bike.y[ys]), bike.xvel[ys, xs], bike.yvel[ys, xs], scale=5000)
    
    ax.plot(*circ(750), color='magenta', ls='--', lw=1)
    #ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    
    dx = km(bike.x[1] - bike.x[0])
    ax.set_title(f't = {bike.time} a, dx_min = {dx} km')





def areas(bike):
    
    dx = km(bike.x[1] - bike.x[0])
    af = np.sum(bike.ice_frac)*dx**2
    ag = np.sum(np.where(bike.beta*bike.thk > 1.e-10, 1.0, 0.0))*dx**2
    
    time = bike.time
    return time, af, ag

def plot_ts(filespec, level):

    files = sorted(glob.glob(filespec))   
    fig = plt.figure(figsize=(8,8), dpi=200)
    
    axi = fig.add_subplot(2,2,1, aspect='equal') 
    initial = BisiclesData(files[0], level=level,iord=0)
    plot_state(initial, axi)
    
    axf = fig.add_subplot(2,2,2, aspect='equal') 
    final = BisiclesData(files[-1], level=level,iord=0)
    plot_state(final, axf)
    
    
    axa = fig.add_subplot(2,1,2) 
    step = max(1,int(len(files)/32))
    a = np.array([areas(BisiclesData(f,level=level,iord=0)) for f in files[::step] ])
      
    axa.plot(a[:,0],a[:,1],'o-',label='all')
    axa.plot(a[:,0],a[:,2],'o-',label='grounded')
    axa.set_xlabel('time (a)')
    axa.set_ylabel('area (km^2)')
    axa.set_xlim(0,10000)
    axa.legend()
    axa.axhline(np.pi*750**2,color='k',lw=1,ls='--')
    axa.set_title(filespec)
    
    return a[:,0],a[:,1],a[:,2]
      
t0, a0, b0 = plot_ts('plot.camip_fulldome_spin_tc-1.??????.2d.hdf5', level=0)    
t1, a1, b1 = plot_ts('plot.camip_fulldome_spin_tc0.??????.2d.hdf5', level=1)   
t2, a2, b2 = plot_ts('plot.camip_fulldome_spin_tc1.??????.2d.hdf5', level=2)   
t3, a3, b3 = plot_ts('plot.camip_fulldome_spin_tc2.??????.2d.hdf5', level=3)  
#t4, a4, b4 = plot_ts('plot.camip_fulldome_spin_tc3.??????.2d.hdf5', level=4) 
# %%
fig = plt.figure(figsize=(4,4), dpi=200)
#plt.plot(ot0,ob0**0.5,'o--',label='0 (old)')
#plt.plot(ot1,ob1**0.5,'o--',label='1 (old)')
plt.plot(t0,a0**0.5,'s-',label='0',ms=2)
plt.plot(t1,a1**0.5,'s-',label='1',ms=2)
plt.plot(t2,a2**0.5,'s-',label='2',ms=2)
plt.plot(t3,a3**0.5,'s-',label='3',ms=2)
plt.axhline((np.pi*750**2)**0.5)
#plt.plot(t4,b4**0.5,'s-',label='4',ms=2)
plt.xlim(0,10000)
plt.ylim(1328,1332)
plt.legend(title='n ref')
plt.xlabel('time')
plt.ylabel(' sqrt(grounded area) (km)')
plt.grid()

