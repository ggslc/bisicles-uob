# BISICLES Run-time configuration, input, output, and checkpoints

## [Run-time configuration (inputs. * files)](#config)

The majority of BISICLES options are set in a configuration file which
is specified at run time as the first command line argument of the
[driver](driver.md) and [solverBenchmark](benchmark.html) programs, or
in the initialization step of the programmable [cdriver
interface](cdriver.md). In the case of driver, something like

    nohup mpirun -np 4 $BISICLES_HOME/BISICLES/code/exec2D/driver2d.Linux.64.mpic++.gfortran.DEBUG.OPT.MPI.ex inputs.mysim

would specify inputs.mysim to be the configuration file. Some options
may be set outside of this file, such as petsc solver options (which can
be set in a file .petsrc in the working directory). Within the
configuration file, options (with optional comments) are written like

    amr.plot_interval = 4 #  write a plot file every 4 steps
    JFNKSolver.maxIter = 10 # cap the JFNK solver to 10 outer iterations
    amr.ref_ratio = 2 2 2 2 2 2 2 2 2 2 2 # dx on level n is twice that of level n-1, for all n

Later entries supersede earlier entries, e.g

    amr.plot_interval = 4 #  write a plot file every 4 steps
    amr.plot_interval = 2 #  write a plot file every 2 steps

is the same as just

    amr.plot_interval = 2 #  write a plot file every 2 steps

Note that BISICLES can run in user-defined time units: the default is
the tropical year. To chose seconds instead, set

    constants.seconds_per_unit_time = 1.0

In that case, make sure that e.g surface mass balance is specified in
m/s (as opposed to m/a). Likewise, ensure that the rate factor (A in
Glen 's law) and the basal traction coefficient are given in the correct
units.

## [Log output (pout. * files)](#pout)

Non-mpi executables write a log is to stdout, while mpi executables
write to files called (by default) pout.0, pout.1, one for each
processor. These are fairly verbose by default: grep is your friend, and
are the most effective way to be sure that a simulation is behaving (e.g
the progress of the velocity solver, which is the hardest and most time
consuming part of the code is reported in some detail). Optionally
specify the file name (for mpi binaries) with

    main.poutBaseName = pout.something

to get pout.something.0, pout.something.1, etc. Either absolute or
relative paths can be specified.

## [Plot files (plot. *.2d.hdf5 files)](#plot)

Ice sheet snapshots, containing spatially varying fields such as ice
thickness and velocity can be written either at regular time steps or
times. For example, set

    amr.plot_interval = 32
    amr.plot_prefix = plot.mysim.

to produce a file every 32 timesteps, with a names like
plot.mysim.000000.2d.hdf5, plot.mysim.000032.2d.hdf5,
plot.mysim.000064.2d.hdf5 , etc Note that, in most cases the time step
is variable, so the file names will not given much information about the
times between them. It is possible, but not recommended, to fix the time
step, with, e.g

    amr.fixed_dt = 0.125

To have plot files written at particular time interval, set

    amr.plot_time_interval = 1.0 # time between plot file in years

and do not set (e.g)

    amr.plot_interval = 32

Provided the time interval is { ... 1/4 , 1/2, 1, 2, 4,  ...} it might
also be worth setting

    amr.time_step_ticks = 1

which will restrict the time step to a value in { ... 1/4 , 1/2, 1, 2,
4,  ...} ,

There are several options that determine which fields are included in
the plot files. With default values they are

    amr.reduced_plot = false; # if true, write a reduced set of fields
    amr.write_solver_rhs = false; # if true, include the right hand side of the stress equation (gravity)
    amr.write_dHDt = true; # if true, write the ice thickening rate
    amr.write_fluxVel = true; # if true, write the vertically averaged mass flux velocity
    amr.write_viscousTensor = false; # if true, write the effective drag and viscosity, plus components of the stress tensor
    amr.write_baseVel = true; # if true, write the velocity at the ice base
    amr.write_internal_energy = false; # if true, include the internal energy density of each vertical layer, plus surface and base 
    amr.write_thickness_sources = false; # if true, write the surface and basal thickness fluxes (SMB and melt-rate)
    amr.write_layer_velocities = false; # if true, write the horizontal velocity at each vertical layer
    amr.write_mask = false;# if true, write the grounded/floating/land/sea mask

By default, files are numbered with the timestep (e.g
plot.mysim.000032.2d.hdf5 corresponds timestep 32), but it is possible
to have them numbered otherwise: the options are

      amr.output_file_numbering = time_yyyymmdd_360 # divides the year into 12 months of 20 days e.g 20010401 for 1 April 2001
      amr.output_file_numbering = time_years # a long integer time in years
      amr.output_file_numbering = time_seconds # a long integer time in seconds

## [CF Plot files (plot. *.CF.2d.hdf5 files)](#plotcf)

CF plot files are alternative (or complement) to the [standard plot file
format](#plot) that are designed to produce [Climate and Forecast
(CF)](http://cfconventions.org/) compatible data. They differ from the
standard files in several respects:

1.  They include a single mesh level, and data from coarser or finer
    levels is interpolated/coarsened as needed)
2.  The fields are by default time-averaged rather than snapshots: each
    file contains the time means since the previous file was written.
    Alternatively, you can specify snapshots (final values)
3.  Time sequences of domain wide data (such as ice volume) are
    included: these span the time interval since the previous file was
    written
4.  The output fields (ice thickness etc) are labeled with CF short
    names

The files produced are not themselves CF-compliant, because they are
hdf5 rather than netcdf files, but the flatten [file
tool](filetools.md) can create CF compliant netcdf files from the a
plot. *.CF.2d.hdf5 file. To convert plot.X.CF.2d.hdf5, run (for example)

      $BISICLES_HOME/BISICLES/code/filetools/flatten2d.Linux.64.g++.gfortran.DEBUG.OPT.ex plot.X.CF.2d.hdf5 plot.X.nc 0

which will also construct some metadata, e.g specify the projection. To
enable CF plot files, disable the standard files and include some
typical data:

    amr.plot_style_cf = true
    amr.plot_style_amr = false # leave out to have both sorts of file
    CFIO.level = 1 # default is 0 (coarsest level); set to desired level

    CFIO.lithk = true # ouput land ice thickness Note the form: CFIO.short_name = true
    CFIO.lithk.time_integration = 1 # output time step final value for lithk. Default is 0 (means)

    CFIO.acabf = true # surface mass balance
    CFIO.acabf.time_integration = 0 # default is 0 (time mean)

    CFIO.velbase = true # for land_ice_basal_velocity
    CFIO.orog = true # for surface_altitude
    CFIO.topg = true # for bedrock_altitude

    #specify the coordinate system - optional for BISICLES but needed for CF compliance
    CRS.EPSG = 3031 # EPSG is the only system supported for now
    CRS.origin_x = 1.234 # the coordinates of the point (0,0) on (all) BISICLES levels 
    CRS.origin_y = 5.678 # (measured in metres)

    #CFIO requires additional data on restarts (see below for restarts in general).
    #when plot_style_cf = true, a *.CF.2d.hdf5 file is written in addition to the usual *.2d.hdf5 checkpoint 
    CFIO.imply_restart = true # default false. If true, imply a filename  (chk.*.CF.2d.hdf5) from amr.check_prefix (chk.*.2d.hdf5)  
    #amr.restart_file = choose a specific CF filename
    CFIO.checkpoint_strict = true # default true. rause error on restarts without a cf checkpoint
    #if false, ignore missing checkpointt data but time integration lower limit will be the restart time   

## [Checkpoints (chk. *.2d.hdf5 files) and restarts](#restart)

BISICLES will periodically write a checkpoint file. Although these are
hdf5 files, they cannot be viewed in visit etc and are intended only
allow simulations to pick up from a point during a previous simulation.
To write a checkpoint every 32 time steps, for example set

    #check points
    amr.check_interval = 32
    amr.check_prefix = chk.mysim.
    amr.check_overwrite = 0

which will produce files chk.mysim.000000.2d.hdf5,
chk.mysim.000032.2d.hdf5, chk.mysim.000064.2d.hdf5, etc. These files can
be large, and I/O can be expensive, so avoid writing too many. By
setting

    amr.check_overwrite = 1

each checkpoint will be called chk.mysim.2d.hdf5 and overwrite the
previous file. We advise avoiding this option. Firstly, if the program
aborts midway through I/O (through running out of time on a cluster, or
disk space), the file will be corrupt, and at least if there are a
series you can go back to the previous file. Secondly, keeping a
reasonable number of checkpoints allows you to construct s set of
branching simulations not necessarily envisaged to begin with.

To restart a simulation from a checkpoint add the option

    amr.restart_file = <filename>

to the configuration file. An obvious use for checkpoints is to create a
number of simulation branches, having carried out some sort of
relaxation. In that case, it might be desirable to set

    amr.restart_time = 0.0
    amr_restart_set_time = true

An effective strategy when running simulations on clusters that impose a
time limit is to write checkpoint files frequently, and then write job
submission scripts that look for the latest checkpoint file and append
options to the input file. Make sure you also append

    amr_restart_set_time = false

The following script template was designed for ARCHER. The script is
submitted from a given directory (which will set  $PBS_O_WORKDIR). When
it runs, each job creates its own numbered working sub-directory, so
that log output is separated. plot and checkpoint files are written to
../ - which will be  $PBS_O_WORKDIR. Subsequent runs of the same script
will be restarted from the last checkpoint. If you copy this script,
make sure you replace all the  <variables >

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
    export INFILE=@INFILE.$JOBNO
    cp $INFILEBASE $INFILE

    #work out what the latest checkpoint file is (if it exists)
    if test -n "$(find ../ -maxdepth 1 -name 'chk.<simulation name>.??????.2d.hdf5' -print -quit)"
        then
        LCHK=`ls -th ../chk.<simulation name>.??????.2d.hdf5 | head -n 1`
        echo "" >> $INFILE #ensure line break
        echo "amr.restart_file=$LCHK" >> $INFILE
        echo "amr.restart_set_time=false" >> $INFILE
        echo "" >> $INFILE #ensure line break
    fi

    aprun -n <number of cores> $DRIVER $INFILE
