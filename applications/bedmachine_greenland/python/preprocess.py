

EXTERNAL_DATA_PATH='../external_data/'
INTERMEDIATE_DATA_PATH='../intermediate_data/'
BEDMACHINE_NC='{}/{}'.format(EXTERNAL_DATA_PATH,'BedMachineGreenland-v5.nc')
MEASURES_NC='{}/{}'.format(EXTERNAL_DATA_PATH,'Greenland_ice_speed_v2017.nc')

def name(s):
    return '{}/greenland_bedmachine_{}.nc'.format(INTERMEDIATE_DATA_PATH,s)

OUTPUT_NC = name('150m')
N_LAYER=24

import os, time

from greenland_projection import add_projection_attr_greenland
from coarsen import coarsen_2D, coarsen_1D
from netCDF4 import Dataset

def create_new_only(file_name, fun):
    if os.path.isfile(file_name):
        mb = int(os.path.getsize(file_name)/1024/1024)
        mt = time.ctime(os.path.getsize(file_name))
        print ('{} ({} Mb {}) exists - delete if you want to recreate it'.format(file_name, mb, mt ))
    else:
        print ('creating {} ...'.format(file_name))
        fun(file_name)
        print ('...done')

def preprocess_150m(output_nc):        
    from preprocess_thk_bed_btrc import preprocess
    preprocess(output_nc, BEDMACHINE_NC, MEASURES_NC)
    return None

def preprocess_therm(output_nc):
    from preprocess_therm_bc import preprocess    
    preprocess(output_nc)
    return None
    
def coarsenc(name, fine_name, v_names = ['thk','topg','umod','btrc','umodc']):
    """
    read uniform cell-centred mesh nc data
    """
    fine_nc = Dataset(fine_name,'r')
    coarse_nc = Dataset(name, 'w')
    
    x_fine = fine_nc.variables['x'][:]
    y_fine = fine_nc.variables['y'][:]
    x_coarse, y_coarse = coarsen_1D(x_fine), coarsen_1D(y_fine)

    xdim = coarse_nc.createDimension('x',size=len(x_coarse))
    ydim = coarse_nc.createDimension('y',size=len(y_coarse))
    
    xv = coarse_nc.createVariable('x','f8',('x'))
    xv[:] = x_coarse

    yv = coarse_nc.createVariable('y','f8',('y'))
    yv[:] = y_coarse
   
    
    for v in v_names:
        vv = coarse_nc.createVariable(v  ,'f8',('y','x'))
        vv[:,:] = coarsen_2D(fine_nc.variables[v][:,:])
    

    add_projection_attr_greenland(coarse_nc, xv, yv)


    
                        
create_new_only(OUTPUT_NC, preprocess_150m)

#create coarse grid versions for later convenience
suffix = ['150m','300m','600m','1200m','2400m','4800m']

for i in range(1,len(suffix)):
    fine = name(suffix[i-1])
    coarse = name(suffix[i])
    def f(file_name):
        from coarsen import coarsenc
        coarsenc(file_name , fine)
        return(None)
    create_new_only(coarse, f)

#create_new_only('greenland_bedmachine_therm_bc_4800m.nc', preprocess_therm)
