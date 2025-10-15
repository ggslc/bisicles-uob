"""

produce coarser resolution versions of netcdf files
for convenience

"""

from netCDF4 import Dataset
import numpy as np

def coarsen_1D(u):
    """
    coarsen u(1D) by a factor of two
    """

    n = len(u)
    uc = 0.5 * ( u[0:n-1:2] + u[1:n:2] )

    #print ('u[1] - u[0] = {} {}'.format(type(u[0]), u[1]-u[0]))
    return uc

def coarsen_2D(u):
    """
    coarsen u(2D) by a factor of two
    """

    n,m = u.shape

    uc = 0.25 * ( u[0:n-1:2, 0:m-1:2]
                + u[1:n:2,   0:m-1:2]
                + u[0:n-1:2, 1:m:2]
                + u[1:n:2,   1:m:2])

    return uc

def add_projection_attr_greenland(ncout, xv, yv):
    """

    Add projections data to a greenland netcdf (assuming EPGS 3143)

    Parameters
    ----------
    ncout : Dataset object 
    xv : x  variable in ncout
    yv : y variable in ncout

    Returns
    -------
    None.

    """
    crsv =  ncout.createVariable('crs','int')
    EPSG = 3143

    ncout.setncattr('Conventions','CF-1.7') 
    crsv.setncattr('EPSG',int(EPSG))
    crsv.setncattr('grid_mapping_name','polar_stereographic')
    crsv.setncattr('latitude_of_projection_origin', 90.0)
    crsv.setncattr('straight_vertical_longitude_from_pole', -45.0)
    crsv.setncattr('scale_factor',1.0)
    crsv.setncattr('standard_parallel',70.0)
    crsv.setncattr('false_easting',0.0)
    crsv.setncattr('false_northing',0.0)
    xv.setncattr('standard_name','projection_x_coordinate')
    xv.setncattr('units','meter')
    yv.setncattr('standard_name','projection_y_coordinate')
    yv.setncattr('units','meter')

    
def coarsenc(name, fine_name, add_projection_f):


    v_names = ['thk','topg','uo','uc','btrc']
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
    

    add_projection_f(coarse_nc, xv, yv)
    coarse_nc.close()
