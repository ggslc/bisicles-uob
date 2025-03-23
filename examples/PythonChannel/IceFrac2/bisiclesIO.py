#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 10:49:03 2020

@author: scornford
"""
from amrfile import io as amrio
amrio.freeAll()
import numpy as np
from netCDF4 import Dataset

import os

def _system(cmd):
    print(cmd)
    os.system(cmd)
   

def nctoamr(nc_file, hdf5_file, var_string):
    _system('nctoamr {} {} {}'.format(nc_file, hdf5_file, var_string))
    
class BisiclesData:
    """
    Store data from a BISICLES output and provide
    derived data (e.g speed from velocity)
    """
    
    def __init__(self,file_name, level=0, origin=(0,0), iord=1,
                 croplo = (0,0), crophi = (1,1)):
        """

        """
    
        amrid = amrio.load(file_name)
        self.time = amrio.queryTime(amrid)
        lo,hi = amrio.queryDomainCorners(amrid, level)

        lo_0 = lo
    
        for dir in [0,1]:
            L = hi[dir] - lo[dir]
            hi[dir] = int( lo_0[dir] + crophi[dir]*L )
            lo[dir] = int( lo_0[dir] + croplo[dir]*L )
            
        def read(name):
            return amrio.readBox2D(amrid, level, lo, hi, name, iord)
        
        self.x,self.y,self.thk = read('thickness')


        self.dx = self.x[1] - self.x[0]
        #adjust x,y origin
        self.x += origin[0]
        self.y += origin[1]
        
        x,y,self.usrf = read('Z_surface')
        x,y,self.topg = read('Z_base')
        x,y,self.xvel = read('xVel')
        x,y,self.yvel = read('yVel')
        x,y,self.beta = read('dragCoef')
        x,y,self.frac = read('iceFrac') 
        x,y,self.acab = read('surfaceThicknessSource')
        x,y,self.ocean = read('activeBasalThicknessSource')
 
        amrio.free(amrid)
           
        hf = np.where(self.topg < 0.0, -self.topg*1028.0/917.0, 0.0)
        self.hab = self.thk - hf
        self.speed = np.sqrt(self.xvel**2 + self.yvel**2)
        self.Tb = self.beta * self.speed
    


