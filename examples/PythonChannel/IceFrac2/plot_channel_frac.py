#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import sys
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.gridspec as gridspec
sys.path.append(os.getcwd() + '/')
from bisiclesIO import BisiclesData

#%%

def extract(pattern, level):

    files = sorted(glob(pattern))
    nt = len(files)
    area = np.empty(nt)
    cz_area = np.empty(nt)
    time = np.empty(nt)
    snap_time = [8,16,32,64,128,192]    
    k_snap = 0
    snap = []
    
    for it, f in enumerate(files):
        bike = BisiclesData(f,level=level,iord=0) 
        time[it] = bike.time
        area[it] = np.sum(bike.frac)*bike.dx**2       
        tol = 1.0e-2
        cz_area[it] = np.sum(np.where(np.abs(bike.frac-0.5) < 0.5-tol,1.0,0.0))*bike.dx**2
        
        if k_snap < len(snap_time):
            if len(snap) < k_snap + 1:
                snap.append(bike)
            
            if np.abs(bike.time - snap_time[k_snap]) < np.abs(snap[k_snap].time - snap_time[k_snap]):
                snap[k_snap] = bike
    
            if bike.time > 0.5 + snap_time[k_snap]:
                k_snap += 1



    return time, area, cz_area, snap

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

def snap_plot(bike, ax):
    
    km = 1.0e-3
    tol = 1.0e-3
    pc = ax.pcolormesh(bike.x*km,bike.y*km,np.where(bike.thk > tol, bike.speed*km, np.nan),
                  vmin=0,vmax=1,cmap='gist_ncar_r')
    cr = ax.contour(bike.x*km,bike.y*km,bike.frac,[0.1,0.5,0.9],colors=['c','k','m'],linewidths=2)
    ax.text(4,8,f't={bike.time}')
    ax.set_yticks([])
    ax.set_aspect('equal')
    
    return pc,cr
    
def snap_plot_expt(result, expt, level):
    
    time, area, cz_area, snap = result   
    
    
    fig = plt.figure(figsize=(12,9))

    
    #fig, axs = plt.subplots(len(snap),1,sharex=True)
    #fig.suptitle(f'expt:{expt} AMR:{tc+1}')
    
    #column of long axes +  column of short axes
    gs = gridspec.GridSpec(max(6,len(snap)), 2, width_ratios=[5,1], figure=fig)
        
    
    for j, bike in enumerate(snap):
        pc, cr = snap_plot(bike, fig.add_subplot(gs[j*2]))
        
    #velocity color bar  
    ax = fig.add_subplot(gs[3])
    cax = inset_axes(ax,width=1,height=0.5,loc='center')
    ax.set_axis_off()
    fig.colorbar(pc, cax=cax, label=r'$|u|$ (km/a)', orientation='horizontal',)
    
    #frac contour key
    ax = fig.add_subplot(gs[7])
    xx,yy = np.meshgrid(np.linspace(0,1,10),np.linspace(-0.25,1.25,10))
    cs = ax.contour(yy,[0.1,0.5,0.9],colors=['c','k','m'],linewidths=2)
    ax.clabel(cs)
    ax.set_axis_off()
    ax.set_title(r'$f$')
    fig.suptitle(f'expt {expt} / AMR {level}')
    fig.savefig(f'snapshots_{expt}_AMR{level}.png',dpi=200)
    


#%%
# area & calving zone length plots


fig, axs = plt.subplots(3,2,figsize=(8,9),sharex=True)
plt.subplots_adjust(wspace=0.25)        
for j, expt in enumerate(expts):
    for tc in tcs: 
        level = tc + 1
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
for expt in expts[0:]:
    for tc in tcs[0:]: 
        snap_plot_expt(result[expt][tc], expt, tc+1)



