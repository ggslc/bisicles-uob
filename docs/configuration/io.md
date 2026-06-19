# Run-time Configuration, Input, Output, and Checkpoints

## Run-time Configuration (inputs.* files)

The majority of BISICLES options are set in a configuration file which is specified at run time as the first command line argument of the [driver](../executables/driver.md) and solverBenchmark programs, or in the initialization step of the programmable [cdriver interface](../executables/cdriver.md).

For the driver, a typical invocation looks like:

```bash
nohup mpirun -np 4 $BISICLES_HOME/BISICLES/code/exec2D/driver2d.Linux.64.mpic++.gfortran.DEBUG.OPT.MPI.ex inputs.mysim
```

This specifies `inputs.mysim` as the configuration file.

### Option Format

Within the configuration file, options (with optional comments) are written like:

```
amr.plot_interval = 4 #  write a plot file every 4 steps
JFNKSolver.maxIter = 10 # cap the JFNK solver to 10 outer iterations
amr.ref_ratio = 2 2 2 2 2 2 2 2 2 2 2 # dx on level n is twice that of level n-1
```

Later entries supersede earlier entries. For example:

```
amr.plot_interval = 4 #  write a plot file every 4 steps
amr.plot_interval = 2 #  write a plot file every 2 steps
```

is equivalent to just:

```
amr.plot_interval = 2 #  write a plot file every 2 steps
```

### External Options

Some options may be set outside of the configuration file, such as PETSc solver options (which can be set in a file `.petsrc` in the working directory).

### Time Units

BISICLES can run in user-defined time units. The default is the tropical year. To use seconds instead, set:

```
constants.seconds_per_unit_time = 1.0
```

When using seconds, make sure that surface mass balance is specified in m/s (not m/a), and ensure that the rate factor (A in Glen's law) and basal traction coefficient are given in the correct units.

## Log Output (pout.* files)

Non-MPI executables write logs to stdout, while MPI executables write to files called `pout.0`, `pout.1`, etc. (one for each processor). These files are fairly verbose by default; `grep` is your friend, and they are the most effective way to ensure that a simulation is behaving correctly. The progress of the velocity solver (the hardest and most time-consuming part of the code) is reported in some detail.

Optionally specify the output file name (for MPI binaries) with:

```
main.poutBaseName = pout.something
```

This produces `pout.something.0`, `pout.something.1`, etc. Either absolute or relative paths can be specified.

## Plot Files (plot.*.2d.hdf5 files)

Ice sheet snapshots containing spatially-varying fields such as ice thickness and velocity can be written either at regular time steps or times.

### Interval-based Output

To write plot files every N timesteps:

```
amr.plot_interval = 32
amr.plot_prefix = plot.mysim.
```

This produces files like `plot.mysim.000000.2d.hdf5`, `plot.mysim.000032.2d.hdf5`, `plot.mysim.000064.2d.hdf5`, etc.

Note that in most cases the time step is variable, so the file names will not provide much information about the times between them.

### Fixed Time Step (not recommended)

You can fix the time step (though this is not recommended):

```
amr.fixed_dt = 0.125
```

### Time-interval Output

To have plot files written at particular time intervals:

```
amr.plot_time_interval = 1.0 # time between plot files in years
```

Do not set `amr.plot_interval` when using time-interval output.

If the time interval is {... 1/4, 1/2, 1, 2, 4, ...}, you might also want to restrict the time step to match these values:

```
amr.time_step_ticks = 1
```

### Plot File Contents

Several options determine which fields are included in plot files. Default values are:

```
amr.reduced_plot = false; # if true, write a reduced set of fields
amr.write_solver_rhs = false; # if true, include RHS of stress equation (gravity)
amr.write_dHDt = true; # if true, write ice thickening rate
amr.write_fluxVel = true; # if true, write vertically averaged mass flux velocity
amr.write_viscousTensor = false; # if true, write effective drag, viscosity, stress tensor
amr.write_baseVel = true; # if true, write velocity at ice base
amr.write_internal_energy = false; # if true, include internal energy density per layer
amr.write_thickness_sources = false; # if true, write SMB and melt-rate
amr.write_layer_velocities = false; # if true, write velocity at each vertical layer
amr.write_mask = false; # if true, write grounded/floating/land/sea mask
```

### Output File Numbering

By default, files are numbered by timestep (e.g., `plot.mysim.000032.2d.hdf5` for timestep 32). Other options are:

```
amr.output_file_numbering = time_yyyymmdd_360 # e.g., 20010401 for 1 April 2001
amr.output_file_numbering = time_years # time as long integer in years
amr.output_file_numbering = time_seconds # time as long integer in seconds
```

## CF Plot Files (plot.*.CF.2d.hdf5 files)

CF (Climate and Forecast) plot files are an alternative to standard plot files, designed to produce [Climate and Forecast (CF)](http://cfconventions.org/) compatible data. They differ in several respects:

1. They include a single mesh level; data from coarser or finer levels is interpolated/coarsened as needed
2. Fields are time-averaged rather than snapshots; each file contains time means since the previous file was written
3. Time sequences of domain-wide data (such as ice volume) are included, spanning the interval since the previous file was written
4. Output fields are labeled with CF standard names

The files produced are not themselves CF-compliant (because they are hdf5 rather than netcdf), but the [flatten file tool](../executables/filetools.md#flatten) can create CF-compliant netcdf files. To convert `plot.X.CF.2d.hdf5`:

```bash
$BISICLES_HOME/BISICLES/code/filetools/flatten2d.Linux.64.g++.gfortran.DEBUG.OPT.ex plot.X.CF.2d.hdf5 plot.X.nc 0
```

To enable CF plot files:

```
amr.plot_style_cf = true
amr.plot_style_amr = false # omit to have both sorts of files
CFIO.level = 1 # default is 0 (coarsest level)
CFIO.lithk = true # for land_ice_thickness
CFIO.velbase = true # for land_ice_basal_velocity
CFIO.orog = true # for surface_altitude
CFIO.topg = true # for bedrock_altitude

# Specify coordinate system (optional for BISICLES, needed for CF compliance)
CRS.EPSG = 3031 # EPSG is the only system supported for now
CRS.origin_x = 1.234 # coordinates of point (0,0)
CRS.origin_y = 5.678 # in metres
```

## Checkpoints and Restarts

### Writing Checkpoints

BISICLES periodically writes checkpoint files. Although these are hdf5 files, they cannot be viewed in VisIt and are intended only to allow simulations to resume from a previous point. To write a checkpoint every 32 time steps:

```
amr.check_interval = 32
amr.check_prefix = chk.mysim.
amr.check_overwrite = 0
```

This produces files like `chk.mysim.000000.2d.hdf5`, `chk.mysim.000032.2d.hdf5`, `chk.mysim.000064.2d.hdf5`, etc.

These files can be large and I/O can be expensive, so avoid writing too many.

### Overwriting Checkpoints (not recommended)

You can overwrite the checkpoint file each time:

```
amr.check_overwrite = 1
```

This produces files called `chk.mysim.2d.hdf5` that overwrite the previous file. **We advise against this** because:

1. If the program aborts midway through I/O (running out of time on a cluster, or disk space), the file will be corrupt
2. Keeping a series of checkpoints allows you to construct branching simulations not necessarily envisaged to begin with

### Restarting from a Checkpoint

To restart a simulation from a checkpoint, add to the configuration file:

```
amr.restart_file = <filename>
```

Checkpoints are useful for creating simulation branches after relaxation. In such cases, you might want to reset the simulation time:

```
amr.restart_time = 0.0
amr.restart_set_time = true
```

### Cluster Job Submission Strategy

An effective strategy for cluster simulations with time limits is:

1. Write checkpoint files frequently
2. Use job submission scripts that look for the latest checkpoint and append options to the input file
3. Ensure you also append:

```
amr.restart_set_time = false
```

#### ARCHER Cluster Script Template

The following template is for ARCHER-style clusters. The script is submitted from a given directory (setting `$PBS_O_WORKDIR`). When it runs, each job creates its own numbered working subdirectory for separated log output, while plot and checkpoint files are written to the parent directory.

**Important:** Replace all variables in `<angle brackets>`

```bash
#!/bin/bash
#PBS -l walltime=08:00:00
#PBS -l select=<number of nodes>
#PBS -A <allocation code>
#PBS -j oe

export JOBNO=$PBS_JOBID
export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)

module swap PrgEnv-cray PrgEnv-intel
module load cray-hdf5-parallel
module swap anaconda python-compute

export DRIVER=<path to>/driver2d.Linux.64.CC.ftn.OPT.MPI.INTEL.ex

BASEDIR=$PBS_O_WORKDIR
RUNDIR=$BASEDIR/$PBS_JOBID
mkdir -p $RUNDIR
cd $RUNDIR

export INFILEBASE=$BASEDIR/<configuration file>
export INFILE=$INFILE.$JOBNO
cp $INFILEBASE $INFILE

# Work out what the latest checkpoint file is (if it exists)
if test -n "$(find ../ -maxdepth 1 -name 'chk.<simulation name>.??????.2d.hdf5' -print -quit)"
then
    LCHK=`ls -th ../chk.<simulation name>.??????.2d.hdf5 | head -n 1`
    echo "" >> $INFILE # ensure line break
    echo "amr.restart_file=$LCHK" >> $INFILE
    echo "amr.restart_set_time=false" >> $INFILE
    echo "" >> $INFILE # ensure line break
fi

aprun -n <number of cores> $DRIVER $INFILE
```
