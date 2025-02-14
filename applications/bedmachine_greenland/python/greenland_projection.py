"""
Various metadata for Greenland projection, EPSG3143_WKT
WGS 84 / NSIDC Sea Ice Polar Stereographic North

from QGIS and so GDAL

"""

from netCDF4 import Dataset
import numpy as np

EPSG = 3143
EPSG3143_WKT = '''
WKT
PROJCRS["WGS 84 / NSIDC Sea Ice Polar Stereographic North",
    BASEGEOGCRS["WGS 84",
        DATUM["World Geodetic System 1984",
            ELLIPSOID["WGS 84",6378137,298.257223563,
                LENGTHUNIT["metre",1]]],
        PRIMEM["Greenwich",0,
            ANGLEUNIT["degree",0.0174532925199433]],
        ID["EPSG",4326]],
    CONVERSION["US NSIDC Sea Ice polar stereographic north",
        METHOD["Polar Stereographic (variant B)",
            ID["EPSG",9829]],
        PARAMETER["Latitude of standard parallel",70,
            ANGLEUNIT["degree",0.0174532925199433],
            ID["EPSG",8832]],
        PARAMETER["Longitude of origin",-45,
            ANGLEUNIT["degree",0.0174532925199433],
            ID["EPSG",8833]],
        PARAMETER["False easting",0,
            LENGTHUNIT["metre",1],
            ID["EPSG",8806]],
        PARAMETER["False northing",0,
            LENGTHUNIT["metre",1],
            ID["EPSG",8807]]],
    CS[Cartesian,2],
        AXIS["easting (X)",south,
            MERIDIAN[45,
                ANGLEUNIT["degree",0.0174532925199433]],
            ORDER[1],
            LENGTHUNIT["metre",1]],
        AXIS["northing (Y)",south,
            MERIDIAN[135,
                ANGLEUNIT["degree",0.0174532925199433]],
            ORDER[2],
            LENGTHUNIT["metre",1]],
    USAGE[
        SCOPE["unknown"],
        AREA["World - N hemisphere - north of 60°N"],
        BBOX[60,-180,90,180]],
    ID["EPSG",3413]]
'''


def add_projection_attr_greenland(ncout, xv, yv):
    """

    Add projections data to a greenland netcdf (assuming EPGS 3143)

    Parameters
    ----------
    ncout : netcdf Dataset object 
    xv : x  variable in ncout
    yv : y variable in ncout

    Returns
    -------
    None.

    """
    crsv =  ncout.createVariable('crs','int')
    ncout.setncattr('spatial_ref',EPSG3143_WKT)
    ncout.setncattr('Conventions','CF-1.7') 
    crsv.setncattr('EPSG',int(EPSG))
    crsv.setncattr('crs_wkt',EPSG3143_WKT)
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

