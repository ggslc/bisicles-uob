

EXTERNAL_DATA_PATH='../external_data'
BEDMACHINE_NC='{}/{}'.format(EXTERNAL_DATA_PATH,'BedMachineAntarctica-v3.nc')
MEASURES_NC='{}/{}'.format(EXTERNAL_DATA_PATH,'antarctica_ice_velocity_450m_v2.nc')
IMBIE2_BASINS_NC='{}/{}'.format(EXTERNAL_DATA_PATH,'imbie2_basin_numbers_8km.nc')

INTERMEDIATE_DATA_PATH='.'

def add_projection_attr_antarctica(ncout, xv, yv):
    """

    Add projections data to an antarctic netcdf (assuming EPGS 3031)

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
    EPSG = 3031

    ncout.setncattr('Conventions','CF-1.7') 
    crsv.setncattr('EPSG',int(EPSG))
    crsv.setncattr('grid_mapping_name','polar_stereographic')
    crsv.setncattr('latitude_of_projection_origin', -90.0)
    crsv.setncattr('straight_vertical_longitude_from_pole', 0.0)
    crsv.setncattr('scale_factor',1.0)
    crsv.setncattr('standard_parallel',-71.0)
    crsv.setncattr('false_easting',0.0)
    crsv.setncattr('false_northing',0.0)
    xv.setncattr('standard_name','projection_x_coordinate')
    xv.setncattr('units','meter')
    yv.setncattr('standard_name','projection_y_coordinate')
    yv.setncattr('units','meter')




    

def name(s):
    return '{}/antarctica_bedmachine_{}.nc'.format(INTERMEDIATE_DATA_PATH,s)


def imbie_name(s):
    return '{}/antarctica_bedmachine_imbie2_basins_{}.nc'.format(INTERMEDIATE_DATA_PATH,s)

OUTPUT_NC = name('500m')
N_LAYER=24
OUTPUT_TEMP_DX='4km'
OUTPUT_IMBIE2_DX='4km'
OUTPUT_IMBIE2_NC = imbie_name(OUTPUT_IMBIE2_DX)

import os, time

def create_new_only(file_name, fun):
    if os.path.isfile(file_name):
        mb = int(os.path.getsize(file_name)/1024/1024)
        mt = time.ctime(os.path.getsize(file_name))
        print ('{} ({} Mb {}) exists - delete if you want to recreate it'.format(file_name, mb, mt ))
    else:
        print ('creating {} ...'.format(file_name))
        fun(file_name)
        print ('...done')

def preprocess_500m(output_nc):        
    from preprocess_thk_bed_btrc import preprocess
    preprocess(output_nc, BEDMACHINE_NC, MEASURES_NC,
               add_projection_attr_antarctica)
    return None
                        
create_new_only(OUTPUT_NC, preprocess_500m)

#create coarse grid versions for later convenience
suffix = ['500m','1km','2km','4km','8km']

for i in range(1,len(suffix)):
    fine = name(suffix[i-1])
    coarse = name(suffix[i])
    def f(file_name):
        from coarsen import coarsenc
        coarsenc(file_name , fine, add_projection_attr_antarctica)
        return(None)
    create_new_only(coarse, f)


def preprocess_imbie(file_n4ame):
    from preprocess_imbie2 import imbie2_mask_nc
    imbie2_mask_nc(file_name,  IMBIE2_BASINS_NC, name('4km'))

#create_new_only(OUTPUT_IMBIE2_NC, preprocess_imbie)
