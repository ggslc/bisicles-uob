# BISICLES solver benchmark executable

Since the vast majority of the CPU time in a typical BISICLES run is
spent in the velocity solver, it can be useful to run the solver as a
standalone code to examine performance and convergence issues. The code
in BISICLES/code/benchmark is designed to read in a BISICLES plotfile
written out as a part of a standard BISICLES run and perform the
momentum balance velocity solve based on the configuration in the
plotfile. In order for this to work, the option

    amr.write_solver_rhs = 1

must be set in the input file used by the driver run. Then, the option

    main.filename =  hdf5 filename 

must be set in the inputs file used to run the benchmark code, which
tells the code where to look for the problem setup. To see how this all
works (assuming that Chombo, etc are all set up already), try the
following test:

1.  cd BISICLES/code/benchmark
2.  type  "make regression ". This causes the following to happen:
    1.  If necessary, the driver code in BISICLES/exec2D is compiled
    2.  the driver executable is run with
        BISICLES/benchmark/inputs.regression as the inputs file. In the
        process, a plotfile named
         "stream.regression.preSolve.000000.2d.hdf5 " is written out.
         "preSolve " means that it 's written immediately before
        conducting a velocity solve
    3.  the solverBenchmark executable is then run using the same inputs
        file
    4.  the code reads in the inputs file and the hdf5 file and performs
        what should be an identical velocity solve to the one which was
        performed in the driver.

