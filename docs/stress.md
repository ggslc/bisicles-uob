# BISICLES stresses

Every BISICLES simulation requires both [basal](#basal) and
[englacial](#basal) stress models to be defined. These relate basal and
englacial stresses to the velocity field. There are also options to
modify the [driving](#driving) stress

## [Driving stress](#driving)

The driving stress is determined by the geometry, however there are some
options that can be used to modify it. Input geometry may lead to
isolated ares of implausible driving stress, for example at tall ice
cliffs. The option

      velocity_rhs.max_rhs_dx = 1.0e+10 #(Pa m) :any floating point number > 1

is disabled by default but will ensure that rho  * g  * h  *  |grad(s) |
 * dx  < velocity_rhs.max_rhs_dx 1.0e+10 Pa m: roughly speaking a value
of1.0e+10 Pa m: corresponds to a 1 km high ice grounded cliff.

Driving stresses are modified immediately upstream and downstream of the
grounding line to avoid the combination of steep slopes and zero basal
friction that would otherwise arise. This is the one-sided scheme
described in Cornford, J. Comput. Phys. 2013, and is enabled by default.
It can be disabled with.

      velocity_rhs.gl_correction = false # true is the default

## [Basal stresses](#basal)

The basal stress model is broken down into a spatially varying basal
friction [coefficient](#btrc) beta ^2(x,y) and temperature dependent
[rate factor](#bratefac) A(T), which do not depend on velocity and a
basal friction [relation](#btrr) f ( C, u), where C =   beta ^2/A(T),
which does. A (T) = 1.0 by default. The simplest model is the linear
viscous sliding law tau_b = C u, more complex rules include the usual
power law tau_b =   beta ^2  |u | ^(-2/3) u , and rules that depend on
the effective pressure. The code also include an additional linear basal
traction term, so that the final rule applied is tau_b = f( C, u) + C_0
u with C_0 used to impose, for example, drag from rocky fjord walls.
Note that whatever friction coefficient is chosen, it will only be
applied to grounded ice or floating ice immediately adjacent to an ice
free region whose upper surface lies above the lower surface of the ice
(a fjord wall).

### [Basal friction coefficient](#btrc)

If you have built the doxygen code documentation, the C++ class
hierarchy underlying the basal friction coefficient is described at
[classBasalFriction.md](../code/doc/doxygen/html/classBasalFriction.html).

#### Floating and grounded ice

By default, BISICLES sets beta ^2 = 0 for every cell whose center is
grounded. Setting

      basal_friction.grounding_line_subdivision = 4 # (any integer > 0, but 4 is recommended)

selects a sub-grid interpolation scheme, which proves useful in some
circumstances (see e.g Cornford et at, Ann, Glac. 2016). The thickness
above flotation is interpolated between cell centers to compute a
grounded-ice fraction w, with  < w  < 1, which is used to weight
beta ^2. At some point we may make this the default.

The subgrid scheme typically means that a factor of two coarser mesh can
be used. Note that BISICLES does not suffer from the much more severe
truncation errors noted in Seroussi, The Cryopshere, 2014 when the
subgrid scheme is not used: we attribute this to the one-sided
difference scheme for the driving stress described in Cornford, J.
Comput. Phys. 2013. This is enabled by default, but can be disabled with

      velocity_rhs.gl_correction = false # true is the default

#### Drag from Rocky Walls

Floating Ice flowing between rocky walls (e.g Petermann Glacier)
experiences some drag by default, the coefficient C_0 is computed from
the basal friction coefficient that would apply if the ice was grounded,
and from the contact areas between ice and rock. It is possible to
increase this drag with a scalar parameter.

    wall_drag.basic = true # default true
    wall_drag.extra = 0.0 # default 0.0

#### Drag in thin ice regions

Realistic problems are sometimes made more difficult by the appearance
of fast flowing thin ice, often far from the regions of interest, which
tends to reduce the stable time step. This can sometimes be mitigated by
imposing some additional drag in thin ice regions. Be careful with this
parameter - it can easily slow down ice shelves.

    thin_ice_drag.extra = 10.0 # default 0.0
    thin_ice_drag.thickness = 10.0 # default 0.0

#### Constant Friction

The simplest meaningful beta ^2 is constant in space and time. E.g to
set beta ^2 = 1000 on all grounded ice:

    geometry.beta_type = constantBeta
    geometry.betaValue = 1.0e+3

#### Python Basal Friction

See also the description of the [python
interface](pythoninterface.md).

If beta ^2 can be expressed as a simple function of **local** thickness
and topography, the python interface can compute it. The major advantage
of this method is that the python expression can be readily evaluated on
whatever meshes the AMR scheme throws up without the need for
interpolation. For example, create a function

    #file foo.py
    import math

    def friction(x,y,t,thck,topg):
        friction = 1.01e+3 + 1.0e+3 * math.sin(x * 1.0e+3)*math.sin(y * 1.0e+3)
        return friction   

in a python module and set

    geometry.beta_type = Python
    PythonBasalFriction.module = foo
    PythonBasalFriction.function = friction

in the configuration file.

#### LevelData Basal Friction

See also the description of the [LevelData
interface](leveldatainterface.md).

If beta ^2 cannot be expressed as a simple function, it might be
represented as data on a uniform mesh. Typically, this will be the case
if beta ^2 was computed from observations in some way. BISICLES will
average or interpolate the data as needed as meshes are generated, with
one caveat: the data mesh spacing must be compatible with the AMR
scheme. Typically, the **data level grid** will usually have a
resolution that coincides with one of the AMR levels. For example, to
read the field btrc from a file basal.2d.hdf5:

    geometry.beta_type = LevelData
    inputLevelData.frictionFile = basal.2d.hdf5
    inputLevelData.frictionName = btrc

It is also possible to use a number of files to make a time-dependent
beta ^2. E.g,

    geometry.beta_type = LevelData
    inputLevelData.frictionFileFormat = basal%4d.2d.hdf5
    inputLevelData.frictionFileStartTime = 0.0
    inputLevelData.frictionFileTimeStep = 1.0
    inputLevelData.frictionName = btrc

will interpolate data in time between the files
basal0000.2d.hdf5,basal0001.2d.hdf5, ...

#### MultiLevelData Basal Friction

Essentially the same as having a LevelData basal friction coefficient,
but with a non-uniform mesh. Tricky in practice, but, e.g

    geometry.beta_type = MultiLevelData
    inputLevelData.frictionFile = basal.2d.hdf5
    inputLevelData.frictionName = btrc

#### Others

Other beta ^2 options include sinusoidalBeta, sinusoidalBetay,
twistyStreamx, singularStream, gaussianBump. These exist for testing to
be carried out without the python interface: don 't worry about them.

### [Basal rate factor](#bratefac)

If you have built the doxygen code documentation, the C++ class
hierarchy underlying the rate factor is described at
[classRateFactor.md](../code/doc/doxygen/html/classRateFactor.html)

### [Basal friction relation](#btrr)

If you have built the doxygen code documentation, the C++ class
hierarchy underlying the basal friction relation is described at
[classBasalFrictionRelation.md](../code/doc/doxygen/html/classBasalFrictionRelation.html)

#### [Power law basal friction relation](#btrrpwr)

The power law basal friction relation covers rules of the form tau_b = C
 |u | ^(m-1) u, which includes linear viscous sliding (m=1), and hard
bed sliding (m=1/3), and yield stress sliding (m=0.9999, though this
will not work well). For example, specify the common third power law by
setting tau_b = C  |u | ^(-2/3) u

    main.basalFrictionRelation = powerLaw
    BasalFrictionPowerLaw.m = 0.3333 # for m = 1/3
    BasalFrictionPowerLaw.includeEffectivePressure = false #optional false is the default

If

    BasalFrictionPowerLaw.includeEffectivePressure = true

Then a factor hab ^m (where hab is the thickness above flotation) is
introduced.

#### [Pressure limited basal friction relation](#btrrlim)

The pressure limited law modifies another basal friction relation f(C,u)
to ensure that basal traction cannot exceed the Coulomb friction,
proportional to the effective pressure N. Two forms are implemented, the
version from Tsai 2015, where  |tau_b | = min( a  * N, f(C,u) ), and the
version from Leguy 2014 (also Schoof 2005 and Gagliardini 2007), where
 |tau_b | =

To choose Tsai 2015, set

    main.basalFrictionRelation = pressureLimitedLaw
    BasalFrictionPressureLimitedLaw.coefficient = 0.5 #for example. Coulomb friction coefficient. should be 0 < a < 1 
    BasalFrictionPressureLimitedLaw.model = Tsai
    BasalFrictionPressureLimitedLaw.basalFrictionRelation = powerLaw
    BasalFrictionPowerLaw.m = 0.3333

To choose Leguy 2014, set

    main.basalFrictionRelation = pressureLimitedLaw
    BasalFrictionPressureLimitedLaw.coefficient = 8.0e12  #for example
    BasalFrictionPressureLimitedLaw.model = Leguy
    BasalFrictionPressureLimitedLaw.basalFrictionRelation = powerLaw
    BasalFrictionPowerLaw.m = 0.3333

## [Englacial stresses](#englacial)

The englacial stress model is broken down into a spatially varying
[stiffness factor](#mucoef) phi(x,y) and temperature dependent [rate
factor](#ratefac) A(T), which do not depend on velocity and a
[relation](#mu) f ( A, grad u), which does, to provide an effective
viscosity phi f (A, grad u). The simplest model is a linear rheology
(constant f). Ice sheet models will usually need either Glen 's flow law
or a rule based on it such as the L1L2 rule (Schoof and Hindmarsh 2010).
There is also an experimental rule that combines any other relation with
a continuum damage model.

### [Englacial stiffness (mu) coefficient](#mucoef)

If you have built the doxygen code documentation, the C++ class
hierarchy underlying the stiffness factor is described at
[classMuCoefficient.md](../code/doc/doxygen/html/classMuCoefficient.html)

### [Englacial rate factor](#ratefac)

If you have built the doxygen code documentation, the C++ class
hierarchy underlying the rate factor is described at
[classRateFactor.md](../code/doc/doxygen/html/classRateFactor.html)

### [Englacial Constitutive relation](#mu)

The englacial constitutive relation has three roles. First, it computes
the effective viscosity (mu) given the strain rate (e_ij) and rate
factor (A), so that the stress components are given by t_ij = mu  *
e_ij. Secondly, it computes the rate of strain heating (usually mu  *
e_ij  * e_ji). Thirdly, it is currently responsible for reconstructing
the velocity field u(x,y,z) from the basal velocity u_b(x,y), though
this is slated for change when non-vertically integrated stresses are
supported.

If you have built the doxygen code documentation, the C++ class
hierarchy underlying the englacial constitutive relation is described at
[classConstitutiveRelation.md](../code/doc/doxygen/html/classConstitutiveRelation.html)

#### Glen 's flow law

To choose Glen 's flow law, set

    main.ConstitutiveRelation = GlensLaw 
    GlensLaw.n = 3.0 # 3.0 is the default
    GlensLaw.epsSqr0 = 1.0e-12
    GlensLaw.delta = 0.0 #

#### The LlL2 flow law

The L1L2 law is a variant of Glen 's flow law used to approximate
vertical shear strains in vertically integrated models.

    main.ConstitutiveRelation = L1L2
    L1L2.n = 3.0 # 3.0 is the default
    L1L2.epsSqr0 = 1.0e-12
    L1L2.delta = 0.0 #
    L1L2.solverTolerance =1.0e-6
    L1L2.effectiveViscositySIAGradSLimit = 1.0e-2 # limit grad(s) in the mass flux calculation
    L1L2.additionalVelocitySIAGradSLimit = 1.0e+10 # limit grad(s) in the effective viscosity
    L1L2.additionalVelocitySIAOnly = false # default is true
    L1L2.startFromAnalyticMu = false       # can be true for n = 3
