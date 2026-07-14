# BISICLES driver executable

The BISICLES _driver_ executable is used to run standalone time-dependent
ice sheet simulations. At the start of a simulation, it carries out the
following operations

1.  Read a [run-time configuration](io.md#config) file, specified as
    the first command line argument
2.  Compute initial thickness, topography, internal energy, stiffness
    and basal friction coefficient data over the whole domain on the
    base level (l = 0) mesh. This data might be loaded from a file
    through the [LevelData interface](leveldatainterface.md), computed
    through the [Python interface](pythoninterface.md) or through a
    number of special purpose modules.
3.  Set L = 0 and solve the stress-balance equations to find the
    velocity field on the l = 0 mesh
4.  Generate a new mesh with L + 1 levels on which some regions are more
    finely resolved according to the refinement criteria
5.  Compute the thickness, topography, internal energy, damage and basal
    friction coefficient data on the L + 1 mesh
6.  Solve the stress-balance equations to find the velocity field on the
    L + 1 mesh
7.  Write a plot file, including the thickness, topography, surface
    elevation, velocity etc
8.  if L has not reached its maximum value, set L = L + 1 and repeat 3-7

It then runs the following operations until the simulation has completed

1.  Write a check-point file, if required
2.  If time t \> simulation time, exit, else continue
3.  Compute surface fluxes (accumulation and melting) given the geometry
    and velocity field
4.  Compute horizontal thickness and internal energy fluxes given the
    geometry, velocity field and surface fluxes
5.  Evolve the thickness and internal energy fields for a time
    determined by the CFL condition
6.  Solve the stress-balance equations to find a new velocity field
7.  Write a plot file, if required
8.  Update the mesh according to the refinement criteria
9.  If needed, solve the stress-balance equations to find a new velocity
    field
10. repeat from step 1

It is also possible to [restart](io.md#restart) a simulation from a
checkpoint file, for example to create a number of simulation branches,
having carried out some sort of relaxation
