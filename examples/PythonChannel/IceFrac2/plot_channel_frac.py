#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import sys
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.getcwd() + '/')
from bisiclesIO import BisiclesData


def extract(pattern, level):

    files = sorted(glob(pattern))
    nt = len(files)
    area = np.empty(nt)
    cz_area = np.empty(nt)
    time = np.empty(nt)
    snap_time = [8,16,32,64,128,256]    
    k_snap = 0
    snap = []
    
    for it, f in enumerate(files):
        bike = BisiclesData(f,level=level,iord=0) 
        time[it] = bike.time
        area[it] = np.sum(bike.frac)*bike.dx**2       
        tol = 1.0e-2
        cz_area[it] = np.sum(np.where(np.abs(bike.frac-0.5) < 0.5-tol,1.0,0.0))*bike.dx**2
        
        if len(snap) < k_snap + 1:
            snap.append(bike)
        
        if np.abs(bike.time - snap_time[k_snap]) < np.abs(snap[k_snap].time - snap_time[k_snap]):
            snap[k_snap] = bike

        if bike.time > 0.5 + snap_time[k_snap]:
            k_snap += 1


    return time, area, cz_area, snap


def snap_plot(bike, ax):
    
    km = 1.0e-3
    tol = 1.0e-3
    ax.pcolormesh(bike.x*km,bike.y*km,np.where(bike.thk > tol, bike.speed, np.nan),
                  vmin=0,vmax=1024,cmap='gist_ncar_r')
    ax.contour(bike.x*km,bike.y*km,bike.frac,[0.1,0.5,0.9],colors=['c','k','m'],linewidths=1)
    ax.text(0,8,f't={bike.time}')
    ax.set_yticks([])
    

#%%
# read results from plot files for all expts/refinements
expts = ['rn','ri','rp']
tcs = [-1,0,1,2,3]
#tcs = [-1]
result = {}

for j, expt in enumerate(expts):
    result[expt] = {}
    for tc in tcs:
        level = tc + 1
        result[expt][tc] =  extract(f'plot.channel_frac.{expt}.tc{tc}.??????.2d.hdf5', level)
 
#%%
# area & calving zone length plots
fig, axs = plt.subplots(3,2,figsize=(8,9),sharex=True)
plt.subplots_adjust(wspace=0.25)        
for j, expt in enumerate(expts):
    for tc in tcs: 
        time, area, cz_area, snap = result[expt][tc]
        axs[j,0].plot(time, area*1.0e-6,label=f'{level}')   
        axs[j,0].set_ylabel(r'frac area (km^2)')

        axs[j,1].plot(time, cz_area/16.0e+3/1e3,label=f'{level}')    
        axs[j,1].set_ylabel(r'cf zone length (km) (->0)')

    axs[j,0].legend(title='AMR')  
    axs[j,0].set_title(f'{expt}')
    axs[j,1].set_title(f'{expt}')
      
axs[j,0].set_xlabel('time')
axs[j,1].set_xlabel('time')

#%% snap shot plots
for expt in expts:
    for tc in tcs: 
        time, area, cz_area, snap = result[expt][tc]    
        fig, axs = plt.subplots(len(snap),1,figsize=(12,9),sharex=True)
        fig.suptitle(f'expt:{expt} AMR:{tc+1}')
        for j, bike in enumerate(snap):
            snap_plot(bike, axs.flat[j])




