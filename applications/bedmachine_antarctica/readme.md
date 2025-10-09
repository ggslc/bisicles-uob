
# BISICLES Antarctic setup with Measures Bedmachine V3

## Steps

1. Data
	
	Obtain the external data, place in the external_data subdir
 
	- antarctica_ice_velocity_450m_v2.nc (NASA measures velocities, https://nsidc.org/data/nsidc-0484/versions/2)
	- BedMachineAntarctica_v3.nc (Bed machine v3, https://nsidc.org/data/nsidc-0756/versions/3)

2. Preprocess

	Compute the fields needed for a basic inverse problem 
	- ice thickness (`thk`) 
        - bedrock elevation (`topg`), 
        - ice speed observations (`uo`)
        - ice speed observations mask (`uc`)
        - initial guess for basal traction coefficient (`btrc`). Assumes m = 1 (linear viscous sliding)

	`> cd intermediate_data`
	`> python3 preprocess.py`

	Produces antarctica_bedmachine_500m.nc, and coarsened verions antarctica_bedmchine_1km.nc (2km,4km,8km)

3. Thermal spinup

	Run a 100,000 year, low (8km) resolution simulation with fixed velocity and climate. 
	Produces a 3D internal energy (enthalpy) field E = Lw + CT  

	`> cd thermal_spinup`
	`> python3 mk_hdf5.py`
	`> bisicles inputs.bm_ant_tspin`

	






