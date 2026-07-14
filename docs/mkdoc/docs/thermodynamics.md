::: (#head)
-   [Index page](index.md)
-   [Top of page](#top)

# Contents

1.  [Temperature and water fraction](#Tw)
2.  [Till water](#till)
3.  [Coupling with the stress model](#coupling)
    1.  [Englacial Stress](#stress)
    2.  [Basal Stress (Weertman)](#bstressW)
    3.  [Basal Stress (Coulomb)](#bstressC)
:::

::: (#main)
# BISICLES Thermodynamics

Temperature and water fraction

BISICLES has an optional thermodynamics component along the lines of
Aschwanden et al, 2012 (doi: 10.3189/2012JoG11J088). It can be enabled
by setting

    amr.isothermal = false

Following the Aschwanden et al 2012 model, the state variable, an energy
density E (measured in J/kg), is composed from the temperature T and the
water fraction w, such that E = CT + Lw, where C is the specific heat
capacity, and L is the specific latent heat of fusion. E is a 3D field:
the ice sheet is subdivided into n layers and n values for E stored at
the centre of each mesh cell. The uppermost layer is labelled layer 0,
and the vertical distribution of layer interfaces is specified through
the amr.sigma option. E.g, for 10 evenly spaced layers, set

      amr.sigma = 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 

There are relatively few options to set for the thermodynamics
component: the initial temperature, and the boundary conditions at the
upper and lower surface. The initial temperature is set in the same
manner as a constant-in-time temperature. In the Pine Island Glacier
example, temperature is read from an hdf5 file with temperatures
specified at the centers of layers which subdivide the ice thickness
evenly, with the options

      temperature.type = LevelData
      inputLevelData.temperatureFile = pig-bisicles-1km.2d.hdf5
      inputLevelData.temperatureName = temp000000

For a uniform initial temperature, set

      temperature.type = constant
      temperature.value = 268

The file contains 10 variables, temp000000 to temp000009, specifying the
temperature at sigma = 0.5,1.5, etc

By default, the boundary condition at the upper surface prescribes zero
heat flux. At the lower surface, the default geothermal flux is also
zero, but the heat flux due to sliding friction is always imposed. The
same means used to specify [surface flux](surfaceflux.md)

for mass can be employed to set these fluxes, with

      surfaceHeatBoundaryData.type = ...
      basalHeatBundaryData.type = ...

The fluxes are measured in J/a/m\^2. It is also possible to set an upper
surface temperature, rather than a heat flux, by setting

      surfaceHeatBoundaryData.Dirichlett = true  # false by default
      surfaceHeatBoundaryData.Temperature = true # true by default

In that case the field specified by SurfaceHeatBoundaryData is taken to
be as a temperature rather than a heat flux.

Drainage and till

Englacial water does not increase in volume indefinitely: once the water
fraction w grows beyond a certain value it is drained to a till layer.
The current model is crude: it is similar (but cruder than) the drainage
model of Aschwanden 2010. The relevant parameters, and their default
values, are:

      ColumnThermodynamics.water_fraction_drain = 0.01 
      ColumnThermodynamics.water_drain_factor = 0.02
      CoulmnThermodynamics.water_fraction_max = 0.05

The first parameter (water_fraction_drain) specifies a water fraction
below which there is no drainage. The second (water_drain_factor)
governs that rate of drainage - which is proportional to the water
fraction, for water fractions up to to third parameter
(water_fraction_max). water above this limit is immediately transferred
to the till.

Water in the till is itself transported elsewhere by the basal hydrology
model. In the simplest, default case (see e.g van Pelt and Oerlemans,
2012, Bueler and van Pelt, 2015) it is simply lost (to a putative ground
water system) at a rate proportional to the till water depth, and
limited to a maximum value. These can be set with

      ColumnThermodynamics.till_water_drain_factor = 0.001 #(default 0.001 1/a)
      ColumnThermodynamics.till_water_max = 4.0 # #(default 4.0 m)

It is possible to set a spatially variable till water drain factor using
the [surface flux classes](surfaceflux.md). The following example
imposes vary rapid drainage in a disc around the origin, preventing till
water from accumulating there.

      #in the inputs.* file
      tillWaterDrainFactor.type = pythonFlux
      tillWaterDrainFactor.module = twc
      tillWaterDrainFactor.function =  till_water_drain_factor

      #in twc.py
      def till_water_drain_factor(x,y,*etc):
        R2 = (240e+3)**2
        factor = 0.005
        fast_factor = 1.0e+3
        r2 = x**2 + y**2
        if (r2 < R2):
            factor = fast_factor

        return factor

A more sophisticated basal hydrology model is in development.

## [](#coupling)Coupling with the stress model

The thermodynamics model is coupled with the stress balance model
through an (optional) temperature dependence of the rate factor that
appear in Glen\'s Law etc ( A(T) ), through an (optional) similar factor
in the basal traction that affects Weertman friction rules, and through
a relationship between till water depth and effective pressure that
affects Coulomb friction rules.

Englacial stress

See also [](stress.md#englacial)

### [Basal stress (Weertman)](#bstressW)

See also [](stress.md#basal)

### [Basal stress (Coulomb)](#bstressC)

See also [](stress.md#basal)
:::
