#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

convert plot files to the CalvingMIP submission format.

1. 321 x 321 5km node cenetred mesh. Corresponds to *nodes* on our
 coarset level, so we must interpolate from cell centers

2. Profile A-H

Profile A (along the line x=0 in the positive y direction)
Profile B (along the line y=x in the positive x and y direction)
Profile C (along the line y=0 in the positive x direction)
Profile D (along the line y=-x in the positive x and negative y direction)
Profile E (along the line x=0 in the negative y direction)
Profile F (along the line y=x in the negative x and y direction)
Profile G (along the line y=0 in the negative x direction)
Profile H (along the line y=-x in the negative x and positive y direction)


"""

import sys
import os
import glob
import numpy as np
from scipy.interpolate import RectBivariateSpline, NearestNDInterpolator
from netCDF4 import Dataset


sys.path.append(os.getcwd() + '/../python')
from bisiclesIO import BisiclesData

def km(xm):
    return xm * 1.0e-3

    
def grid_321():
    x = np.linspace(-0.8e+6, 0.8e+6, 321)
    y = x.copy()
    return x,y

   
def transect(x, y, z, x0, y0, x1, y1):
    #straight transect from point (x0, y0) to (x1, y1)
    
   dx = x[1]-x[0]
   
   R = ((x1-x0)**2 + (y1 - y0)**2)**0.5
   rt = np.arange(0, R, dx)
   cos_theta = (x1-x0)/R
   sin_theta = (y1-y0)/R
   xt = x0 + rt * cos_theta
   yt = y0 + rt * sin_theta
   xx, yy = np.meshgrid(x, y)
   zt = RectBivariateSpline(y,x,z,  kx =3, ky = 3)(yt, xt, grid=False)
   #zt = NearestNDInterpolator(list(zip(yy.flatten(),xx.flatten())),z.flatten())(yt, xt)
   return rt, zt
   

    
    

def bike_to_321_nc(plot_file_list, level, nc_file, preview=False, time_offset=100.0):

    #Jim wants Time100, Visit likes Time   
    time_100_name = 'Time' if preview else 'Time100' 

    #read some data to work out dimension
    bike = BisiclesData(plot_file_list[0],level=level,iord=0)    
    xf0, dummy = transect(bike.x, bike.y, bike.thk, 0.0, 0.0, 800e+3, 0.0)
            

    ncout = Dataset(nc_file,'w')
    
    t100 = np.linspace(0.0,1000.0,11)
    t1 = np.linspace(0.0,1000.0,1001)
    x, y = grid_321()
    
    #dimensions
    
    t1dim = ncout.createDimension('Time1',size=len(t1))
    t100dim = ncout.createDimension(time_100_name,size=len(t100))
    xdim = ncout.createDimension('X',size=len(x))
    ydim = ncout.createDimension('Y',size=len(y))
    theta_list = np.pi*np.arange(0.0, 2.0, 1./4.)
    letter_list = ['A','B','C','D','E','F','G','H']
    
    pdim = [ncout.createDimension(f'Profile {letter}', size=len(xf0)) for letter in letter_list]
    
    
    #var defs
    t1v = ncout.createVariable('Time1','f8','Time1', fill_value=0.0)
    t1v.setncattr('units','a')
    t100v = ncout.createVariable(time_100_name,'f8',(time_100_name) ,fill_value=np.nan)
    t100v.setncattr('units','a')
    xv = ncout.createVariable('X','f8',('X'))
    xv.setncattr('units','m')
    yv = ncout.createVariable('Y','f8',('Y'))
    yv.setncattr('units','m')


    def create2D(name, units):
        v = ncout.createVariable(name,'f8',('X','Y',time_100_name), fill_value=np.nan)
        v.setncattr('units',units)
        return v
    
    def create0D(name, units):
        v = ncout.createVariable(name,'f8',('Time1'), fill_value=np.nan)
        v.setncattr('units',units)
        return v
        
    def create1D(name, letter,  units):
       v = ncout.createVariable(f'{name}{letter}','f8',(f'Profile {letter}','Time1'), fill_value=np.nan)
       v.setncattr('units',units)
       return v
   
    
    #1D, 1 year fields (transects)
    lithkP = [create1D('lithk', letter, 'm') for letter in letter_list]
    topgP = [create1D('topg', letter, 'm') for letter in letter_list]
    usrfP = [create1D('usrf', letter, 'm') for letter in letter_list]
    sP = [create1D('s', letter, 'm') for letter in letter_list]
    xvelmeanP = [create1D('xvelmean', letter, 'm/a') for letter in letter_list]
    yvelmeanP = [create1D('yvelmean', letter, 'm/a') for letter in letter_list]
    maskP = [create1D('mask', letter, '') for letter in letter_list]  
    
    #0D, 1 year fields
    tendlicalvfv = create0D( 'tendlicalvf', 'kg/a')
    tendligroundfv = create0D( 'tendligroundfv', 'kg/a')
    iareaflv = create0D('iareafl', 'm^2')
    iareagrv = create0D('iareagr', 'm^2')
    limv = create0D('lim', 'kg')
    limnswv = create0D('limnsw', 'kg')
    
    
    #2D, 100 year fields
    lithkv = create2D('lithk','m') 
    topgv = create2D('topg','m') 
    xvelmeanv =  create2D('xvelmean','m/a')
    yvelmeanv = create2D('yvelmean','m/a')
    maskv = create2D('mask','')
    

    xv[:] = x
    yv[:] = y
    
    tendlicalvfv[:] = 0 
    
    
    it100 = 0
    it1 = 0
    for it, plot_file in enumerate(plot_file_list):
        print(it, plot_file)
        bike = BisiclesData(plot_file,level=level,iord=0)
        time = bike.time #- time_offset
        dx = bike.x[1] - bike.x[0]
        dx2 = dx*dx
        
        if time in t1:
            #totals
            t1v[it1] = time
            kgm3 = 917.0
            tendlicalvfv[it1] = np.sum(bike.cflux)*dx2 * kgm3
            iareagrv[it1] = np.sum(np.where(bike.hab > 0, 1, 0))*dx2
            iareaflv[it1] = np.sum(np.where(bike.ice_frac > 0.01, bike.ice_frac, 0))*dx2 - iareagrv[it1] 
            limv[it1] = np.sum(bike.thk)*dx2* kgm3
            limnswv[it1] = np.sum(np.where(bike.hab > 0, bike.hab, 0))*dx2* kgm3
            
            #transects
            for k, letter in enumerate(letter_list):
                X = 800.0e3 * np.sin(theta_list[k])
                Y = 800.0e3 * np.cos(theta_list[k])
                xf, f = transect(bike.x, bike.y, bike.ice_frac, 0.0, 0.0, X, Y)
                def maskf(x,z):
                    zm = np.where(f > 0.5, z, np.nan)[0:len(xf0)]
                    xm = x[0:len(xf0)]
                    return xm,zm
                sP[k][:,it1], lithkP[k][:,it1]= maskf(*transect(bike.x, bike.y, bike.thk, 0.0, 0.0, X, Y))
                __, usrfP[k][:,it1] = maskf(*transect(bike.x, bike.y, bike.usrf, 0.0, 0.0, X, Y))
                __, topgP[k][:,it1] = maskf(*transect(bike.x, bike.y, bike.topg, 0.0, 0.0, X, Y))
                __, xvelmeanP[k][:,it1] = maskf(*transect(bike.x, bike.y, bike.xvel, 0.0, 0.0, X, Y))
                __, yvelmeanP[k][:,it1] = maskf(*transect(bike.x, bike.y, bike.yvel, 0.0, 0.0, X, Y))
                __, hab = transect(bike.x, bike.y, bike.hab, 0.0, 0.0, X, Y)
                maskP[k][:,it1] = np.where(hab > 0, 1, np.where(f > 0.5, 2, 3))[0:len(xf0)]
                
                
            #update index
            it1 += 1
        
        if time in t100: 
            print(f'2D, t =  {bike.time - 100}')
            
            def bilin(z):
                return RectBivariateSpline(bike.x, bike.y, z,  kx =1, ky = 1)(x, y)


            mask = np.where(bike.hab > 0, 1, np.where(bike.ice_frac > 0.5, 2, 3))
            mask321 = np.round(bilin(mask),decimals=0)

            def maskf(z):
                return np.where(mask321 < 3, z, np.nan)
            
            #plt.imshow(bilin(bike.cflux))
            t100v[it100] = time
            lithkv[:,:,it100] = maskf(bilin(bike.thk)).transpose()
            xvelmeanv[:,:,it100] = maskf(bilin(bike.xvel)).transpose()
            yvelmeanv[:,:,it100] = maskf(bilin(bike.yvel)).transpose()
            maskv[:,:,it100] = mask321.transpose()
            topgv[:,:,it100] =  bilin(bike.topg).transpose()
            it100 += 1
    
    ncout.close()

files = sorted(glob.glob('plot.camip_fulldome_expt2_tc0.??????.2d.hdf5'))
bike_to_321_nc(files, 1, 'test.nc', preview=True)
