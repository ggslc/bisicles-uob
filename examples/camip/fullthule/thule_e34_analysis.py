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
    th = np.linspace(0.,2.*np.pi,100)
    x = r * np.cos(th)
    y = r * np.sin(th)
    return x,y    


def r_analytic(time):
    dr =  np.where(time > 1000, 750.0/(2*np.pi) * (np.cos(time*2.0*np.pi/1000) - 1.0),0)
    return 750.0 + dr

def plot_state(bike, ax):


    pc = ax.pcolormesh(km(bike.x), km(bike.y), np.ma.masked_array(bike.thk, bike.ice_frac < 0.1), 
                       vmin=0, vmax=1400, cmap='pink_r')
    ax.contour(km(bike.x), km(bike.y), bike.ice_frac, [0.05, 0.1, 0.5, 0.9, 0.95], colors=['r','r','k','r','r'], 
               linewidths=[0.5,0.5,0.5])
    #plt.colorbar(pc, label='h',shrink=0.5,orientation='horizontal')
    ax.contour(km(bike.x), km(bike.y), bike.hab,[0])
    
    w = int(len(bike.x)/16)
    xs = slice(int(w/2), len(bike.x), w)
    ys = slice(int(w/2), len(bike.y), w)

    dt = 100.0
    ax.quiver(km(bike.x[xs]), km(bike.y[ys]), bike.xvel[ys, xs], bike.yvel[ys, xs], scale=1600e+3/dt)
    
    #ax.plot(*circ(750), color='magenta', ls='--', lw=1)
    
    
    #ax.plot(*circ(r_analytic(bike.time)), color='magenta', ls='--', lw=2, label='analytic front')
    #ax.set_xlabel('x (km)')
    ax.set_ylabel('y (km)')
    
    dx = km(bike.x[1] - bike.x[0])
    ax.set_title(f't = {bike.time} a, dx_min = {dx} km')
    

def areas(bike):
    
    dx = km(bike.x[1] - bike.x[0])
    af = np.sum(bike.ice_frac)*dx**2
    #af = np.sum(np.where(bike.ice_frac > 0.5, 1, 0))*dx**2
    ag = np.sum(np.where(bike.hab > 1.e-10, 1.0, 0.0))*dx**2
    
    time = bike.time
    return time, af, ag

def plot_ts(filespec, level):

    files = sorted(glob.glob(filespec))   
    fig = plt.figure(figsize=(8,8), dpi=200)
    
    axi = fig.add_subplot(2,2,1, aspect='equal') 
    initial = BisiclesData(files[0], level=level,iord=0)
    plot_state(initial, axi)
    axi.legend(loc='lower left')
    
    n = len(files)
    axm = fig.add_subplot(2,2,2, aspect='equal') 
    mid = BisiclesData(files[int(3*n/4)], level=level,iord=0)
    plot_state(mid, axm)

    
    axf = fig.add_subplot(2,2,3, aspect='equal') 
    final = BisiclesData(files[-1], level=level,iord=0)
    plot_state(final, axf)
    
    
    axa = fig.add_subplot(2,2,4) 
    step = max(1,int(len(files)/32))
    a = np.array([areas(BisiclesData(f,level=level,iord=1)) for f in files[::step] ])
      
    Rsq = 1600**2
    t = a[:,0]
    axa.plot(t,a[:,1]/Rsq,'o-',label='ice (numerical)',ms=3)
    
    #analytic_area = np.pi*r_analytic(t)**2 
    #axa.plot(t, analytic_area/ Rsq,'--', label='ice (analytic)' )
    #axa.plot(t,(a[:,1] - analytic_area)/analytic_area * 10, '-', label='ice (10 * rel error)' )    


    axa.plot(t,a[:,2]/Rsq,'o-',label='grounded',ms=3)
    axa.set_xlabel('time (a)')
    axa.set_ylabel('area / domain area ')
    #axa.set_xlim(-100,1000)
    #axa.set_ylim(0.,1.0)
    axa.legend(fontsize='xx-small',loc='upper right')
    axa.axhline((np.pi*r_analytic(0)**2) / Rsq,color='k',lw=1,ls='--')
    axa.axhline((np.pi*r_analytic(500)**2) / Rsq,color='k',lw=1,ls='--')
    #axa.set_title(filespec)
    #axa.grid()
    fig.subplots_adjust(wspace=0.15)
    
    return t,a[:,1],a[:,2]
    
    
t0, a0, b0 = plot_ts('plot.camip_fullthule_expt34_tc-1.??????.2d.hdf5', level=0)    
t1, a1, b1 = plot_ts('plot.camip_fullthule_expt34_tc0.??????.2d.hdf5', level=1)   
t2, a2, b2 = plot_ts('plot.camip_fullthule_expt34_tc1.??????.2d.hdf5', level=2)
#t3, a3, b3 = plot_ts('plot.camip_fulldome_expt2_tc2.??????.2d.hdf5', level=3)    
#
