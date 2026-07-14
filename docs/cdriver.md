# [BISICLES cdriver interface](#top)

## [Overview](#over)

The BISICLES cdriver interface allows a C or FORTRAN (or indeed, any
language that can call simple C functions) to control BISICLES without
any direct interaction with C++. The
[UniCiCles](https://puma.nerc.ac.uk/svn/UniCiCles_svn/UniCiCles)
repository packages BISICLES and Chombo with a version of the
Glimmer-CISM ice sheet model that is able to interface with BISICLES
through cdriver. In effect, Glimmer-CISM hands over ice sheet dynamics
to BISICLES but continues to compute GIA, interact with atmosphere
models through glint, etc. FORTRAN models able to interact with
Glimmer-CISM automatically gain the ability to interact with BISICLES.
Source code documentation can be found at
[cwrapper_8H.md](../code/doc/doxygen/html/cwrapper_8H.html) if you
gave built the doxygen documentation.

In terms of the interface itself, it is possible to:

-   create multiple BISICLES instances (for example, to simulate several
    ice sheets in a global model);
-   instruct BISICLES to read e.g surface mass balance data from a
    rectangular array, potentially distributed across MPI ranks in
    parallel programs;
-   initialize BISICLES, including mesh generation;
-   advance the ice sheet to a certain time;
-   read data from BISICLES into rectangular arrays, again, potentially
    distributed across MPI ranks.

A calling program will typically:

1.  Set up MPI, if applicable
2.  Call bisicles_new_instance (C) or f_bisicles_new_instance (Fortran),
    obtaining an integer key that will be used in all subsequent calls.
3.  Allocate memory for data to BISICLES to read, such as surface mass
    balance
4.  Make multiple calls to (f\_)bisicles_set_2d_data to instruct
    BISICLES to read e.g surface mass balance from a particular array.
5.  Call (f\_)bisicles_init_instance: after this point it is possible to
    change the (e.g) surface mass balance data but **not its memory
    address**.
6.  Call f_bisicles_get_2d_data to e.g the ice sheet upper surface
    elevation into rectangular arrays.
7.  Modify (e.g) the surface mass balance data.
8.  Call (f\_)bisicles_advance to instruct BISICLES to advance in time
9.  Repeat steps 6-8 till complete.
10. Shut down the BISICLES instance by calling
    (f\_)bisicles_free_instance.
11. Shut down MPI , if applicable.

## [FORTRAN 90 example](#fortran)

An example FORTRAN 90 program which controls BISICLES via the cdriver
interface is provided with the source code
[testwrapper.F90](../code/cdriver/testwrapper.F90). The usual
combination of binary types (debug,optimized,mpi,petsc) can be built and
linked with the appropriate compiled objects. For example, a
non-optimized serial version can be compiled (assuming CXX=g++ and
FC=gfortran) and run with:

      > cd $BISICLES_HOME/code/cdriver
      > make ftestwrapper DEBUG=TRUE OPT=FALSE USE_PETSC=FALSE MPI=FALSE
      > ./ftestwrapper.2d.Linux.64.g++.gfortran.DEBUG.ex 

For a optimized parallel (non-petsc) version, assuming MPICXX=mpiCC and
FC=gfortan:

      > cd $BISICLES_HOME/code/cdriver
      > make ftestwrapper DEBUG=TRUE OPT=TRUE MPI=TRUE USE_PETSC=FALSE
      > mpirun -np 4 ftestwrapper.2d.Linux.64.mpiCC.gfortran.DEBUG.OPT.MPI.ex
