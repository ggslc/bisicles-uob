#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 12:37:14 2019

post process inverse problem data to produce
basal friction coefs, etc, for spin up / forward runs

differs from 95x by explictly making the raised (~200 m) ice free
region in front of the NEGIS stickier

@author: stephen
"""

from netCDF4 import Dataset
import numpy as np
from  scipy.interpolate import RectBivariateSpline
from scipy import ndimage
from osgeo import ogr, osr
import os

FILETOOLS_DIR='~/Development/bisicles-uob/code/filetools/'
EXCFG='2d.Linux.64.g++.gfortran.DEBUG.OPT.ex'

def filetool(tool, args):
    return FILETOOLS_DIR + tool + EXCFG + ' ' + args


os.system(filetool('extract', 'ctrl.bm_gris_rb_m3.03lev.000000000250.2d.hdf5 post_inverse_tmp.2d.hdf5 Cwshelf muCoef xVelb yVelb Z_base thickness'))

#meshes levels: 4.8km (0) 2.4km (1) 1.2km(2) 0.6km(3)
os.system(filetool('flatten', 'post_inverse_tmp.2d.hdf5 post_inverse_tmp.nc 3'))

ncin =  Dataset('post_inverse_tmp.nc', 'r')

x = ncin.variables["x"][:] -653000
y = ncin.variables["y"][:] -3384500
h =  ncin.variables["thickness"][:,:]
topg =  ncin.variables["Z_base"][:,:]
mucoef =  ncin.variables["muCoef"][:,:]
cthird =  ncin.variables["Cwshelf"][:,:]
u = ncin.variables["xVelb"][:,:]
v = ncin.variables["yVelb"][:,:]
uu = (u**2 + v**2)**0.5
rhoo = 1027.0
rhoi = 917.0
hf = np.where(topg < 0.0, - rhoo/rhoi * topg, 0.0)
hab = h - hf

ny = len(y)
nx = len(x)


#THESE ARE ABRITRARY VALUES THAT MAY BE IMPORTANT IN
#AREAS THAT ARE ICE-FREE IN THE INVERSE PROBLEM. 
cthird_flat_bed = 3.0e+3
delta_cthird_raised_bed = 2.0e+4

raised = topg > -200.0
cthird_ice_free = cthird_flat_bed + np.where(raised, delta_cthird_raised_bed  , 0.0)
cthird = np.where(np.logical_and(hab > 0, h > 10.0), cthird, cthird_ice_free)



# Joughin regularized coeffs. Need a default value of |u| in ice free areas
def jreg(cthird, uu, uf):
    uureg = np.where( h > 1.0, uu, uf)
    return cthird * ( uureg/uf + 1.0 )**(1.0/3.0)

cthird_jreg_300 = jreg(cthird, uu, 300.0)

#ouput netcdf
print ('writing ...')
out='gris_bedmachine5_inverse_ps2'
ncout = Dataset(out + '.nc','w')
#dimensions
xdim = ncout.createDimension('x',size=len(x))
ydim = ncout.createDimension('y',size=len(y))
#var defs

xv = ncout.createVariable('x','f8',('x'))
yv = ncout.createVariable('y','f8',('y'))

def create2D(name):
    v = ncout.createVariable(name,'f8',('y','x'))
    return v


cthirdv = create2D('cthird')
cthird300v = create2D('cthird_jreg_300')
muv = create2D('mucoef')
hv = create2D('thk')
bv = create2D('topg')
uuv = create2D('speed')

#data
xv[:] = x
yv[:] = y
cthirdv[:,:] = cthird
cthird300v[:,:] = cthird_jreg_300
muv[:,:] = mucoef
hv[:,:] = h
bv[:,:] = topg
uuv[:,:] =  uu


ncout.close()

os.system(filetool('nctoamr',out + '.nc ' + out + '.2d.hdf5' 
                   + ' thk topg speed mucoef cthird cthird_jreg_300'))
