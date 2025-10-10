import numpy as np
from netCDF4 import Dataset
#provide inverse problem data given 'observations'
#surface elevation (s)
#bedrock elevation (b)
#thickness (h)
#velocity (ux, uy)

def write_nc(file_name, x , y, var_dict):
    
    #ouput netcdf
    ncout = Dataset(file_name,'w')

    #dimensions
    xdim = ncout.createDimension('x',size=len(x))
    ydim = ncout.createDimension('y',size=len(y))
    xv = ncout.createVariable('x','f8',('x'))    
    yv = ncout.createVariable('y','f8',('y'))

    def create2D(name):
        v = ncout.createVariable(name,'f8',('y','x'))
        return v

    for k in var_dict:
        
        vv = create2D(k)
        vv[:,:] = var_dict[k]

    xv[:] = x
    yv[:] = y

    ncout.close()



nc = Dataset('plot.ideal_shelf_synthetic_data.000000.nc','r')

x, y = nc.variables['x'][:], nc.variables['y'][:]
s = nc.variables['Z_surface'][:,:]
ux = nc.variables['xVel'][:,:]
uy = nc.variables['yVel'][:,:]
h = nc.variables['thickness'][:,:]
b = nc.variables['Z_base'][:,:]

#floating cells
flot = s - h > b

umod = np.sqrt((ux*ux + uy*uy))

# assume uniform mesh of squares - should be true
dx = x[1] - x[0]
dsx, dsy = np.gradient(s,  dx)
ds = np.sqrt(dsx*dsx + dsy*dsy)
rho = 918.0
grav = 9.81
beta_init = np.where(flot, 100.0, rho*grav*ds*h / (umod + 1.0e-1)**(1/3))
beta_init = np.where(beta_init > 2000.0, 2000.0, beta_init)

#noisy data in shelf only
DIRT = 10.0
W = x[-1] - x[0] + dx
yy = np.pi * (x + y) * 32.0/W
umod_dirty = umod + np.where(flot, DIRT * np.sin(yy), 0.0)

#speed masks (grounded, all ice)
H_EPS = 1.0
uc_grounded = np.where(flot, 0.0, 1.0)
uc_all = np.where(h < H_EPS, 0.0, 1.0)

write_nc(f'synthetic_observations_{int(dx)}m.nc', x, y, 
         {"umod_clean":umod, 
          "umod_dirty":umod_dirty,
          "beta_init":beta_init,
          "uc_grounded":uc_grounded,
          "uc_all":uc_all,
          "dsmod":ds,
          } 
         )



