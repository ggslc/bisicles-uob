# BISICLES Python Interface

BISICLES ' python interface is intended to make some of less performance
intensive parts of BISICLES easier to program. It is fairly crude : at
certain points during a run BISICLES will use an embedded python
interpreter to evaluate various fields. At the moment, it is possible to
specify [surface fluxes](#flux) (that is, accumulation and melting), the
[initial geometry](#ibc) (topography and thickness), the [initial
temperature](#ibctemp) and the [basal traction coefficient](#btrc). The
most likely application of the Python interface is in defining idealized
problems, where the use of the [LevelData
interface](leveldatainterface.md) would lead to undesirable numerical
error

You should have included the python interface when [building
BISICLES](readme.md), otherwise you will need to do so. You will need
to run the appropriate [make clean](readme.md#makeclean)

To run BISICLES with the python interface, you need to ensure that the
environment variable PYTHONPATH includes the directory where your python
modules (.py files) are stored. If this is the working directory, then
in serial, run BISICLES like so

    PYTHONPATH=`pwd` driver2d.Linux.64.mpic++.gfortran.DEBUG.OPT.ex inputs.file sout.0 2>err.0 &

and in parallel, assuming 4 processors,

    nohup mpirun -np 4 -x PYTHONPATH=`pwd` driver2d.Linux.64.mpic++.gfortran.DEBUG.OPT.MPI.ex inputs.file &

## [Surface Fluxes](#flux)

To specify a [surface flux](surfaceflux.md) through a python function
write a python function that takes five scalar arguments
(x,y,t,thickness,topography). Let 's say that you have created a file
foo.py

    #file foo.py
    def acab(x,y,t,thck,topg):
        return 1.0e-3 * (thck + topg)

then you would include lines like

    surfaceFlux.type = pythonFlux
    surfaceFlux.module = foo
    surfaceFlux.function = acab

in the input file. It is possible to add a few keyword args (kwargs) to
surfaceFlux.function, e.g, if you want to compute a melt rate that
depends on distance from the grounding line, you could have

    #file foo.py
    import math
    def melt(x,y,t,thck,topg,gl_proximity=0.0,gl_proximity_scale=1.0):
        d = - math.log(gl_proximity + 1.0e-10) * gl_proximity_scale
        return -  (d/1.0e+3)**2 * gl_proximity 

with input file lines

    basalFlux.type = maskedFlux
    basalFlux.floating.type = pythonFlux
    baslaFlux.floating.module = foo
    basalFlux.floating.function = acab
    basalFlux.floating.n_kwargs = 2
    basalFlux.floating.kwargs = gl_proximity gl_proximity_scale

So far the supported list is short

1.  gl_proximity : solution to p + scale ^2  * grad ^2 = 1 (grounded) or
    0 (floating) with natural BCs
2.  gl_proximity_scale : scale in the above

## [Initial Geometry and boundary conditions](#ibc)

To set the initial thickness and topography using python functions,
create a python module contain two scalar functions of (x,y), e.g

    #file foo.py
    def thck(x,y):
         return 1.5e+2

    def topg(x,y):
        return -1.0e+2 - x/1.0e+3

and specify it in the configuration file like so

    geometry.problem_type = Python
    PythonIBC.module = foo
    PythonIBC.thicknessFunction = thck
    PythonIBC.topographyFunction = topg

x and y are double precision floating point numbers, they are given in
meters, and will correspond to the centers of cells. Note that they will
be on the BISICLES co-ordinate system, ie the bottom left corner of the
mesh is (0,0)

The choice of

    geometry.problem_type = Python

implies **boundary conditions as well as initial conditions**. By
default, reflection boundary conditions (ice divides) are imposed on all
four domain edges. This can be changed to periodic boundary conditions
in the usual way, e.g

    amr.is_periodic = 0 1 0

selects reflection boundaries at the x-faces of the domain, but periodic
boundaries at the y-faces. Alternatively , a mixture of Reflection (Ice
Divide) and No Slip conditions can be imposed through the
PythonIBC.bc_lo and PythonIBC.bc_hi options. Choose 0 for reflection ,
and 1 for no-slip. For example

    PythonIBC.bc_lo = 0 1 # ice divide at x = 0 , no slip at y = 0
    PythonIBC.bc_hi = 1 1 # no-slip at x = L , no slip at y = W

To set a marine boundary condition, use DomainEdgeCalvingModel, e.g

    CalvingModel.type = DomainEdgeCalvingModel
    CalvingModel.front_hi = 1 0    #impose a marine boundary at the high x-face 
    CalvingModel.front_lo = 0 0

Alternatively, it is possible to maintain a calving front inside the
domain by setting a large negative accumulation: Adding a function

    #file foo.py
    def melt(x,y,t,thck,topg):
        melt = 0.0
        if (x^2 + y^2 > 1.0e+12):
          melt = -1.0e+3
        return melt

to a python module and the configuration lines

    basalFlux.type = pythonFlux
    basalFlux.module = foo
    basalFlux.function = melt

would lead to a quarter-circular calving front centered on the bottom
left corner with radius 1000 km

## [Initial Geometry and temperature boundary conditions](#ibctemp)

To set the initial temperature using python functions, create a python
module containing a scalar function of (x,y,thickness,topography,sigma),
e.g

    #file foo.py
    def temperature(x,y,thck,topg,sigma):
        T = 263.0
        if (abs(y-64.0e+3) < 8.0e+3):
            T = T + sigma*10.0
        return T

and specify it in the configuration file like so

    temperature.type = Python
    PythonIceTemperatureIBC.module = foo
    PythonIceTemperatureIBC.function = temperature

## [Basal Friction Coefficient](#btrc)

You can specify a basal friction coefficient in much the same way as a
surface flux, e.g add

    #file foo.py
    import math

    def friction(x,y,t,thck,topg):
        friction = 1.01e+3 + 1.0e+3 * math.sin(x * 1.0e+3)*math.sin(y * 1.0e+3)
        return friction   

in a python module and

    geometry.beta_type = Python
    PythonBasalFriction.module = foo
    PythonBasalFriction.function = friction

in the configuration file.
