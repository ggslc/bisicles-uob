# BISICLES LevelData interface

BISICLES LevelData interface can load various kinds of gridded data from
files in order to specify some of 2D and 3D fields. Currently, these
are:

[TOC]

All of these require the same kind of input file, namely a Chombo amr
file stored as hdf5, with a single level of data covering the whole
computational domain. This **data level grid** will usually have a
resolution that coincides with one of the AMR levels specified by the
input file parameters amr.num_cells,amr.maxLevel and amr.ref_ratio. More
generally, it can have a resolution which is 2 ^n (where n is an
integer) times coarser or finer than any of the AMR levels. The simplest
way to produce the correct kind of file is to create a netcdf file and
convert it with the [nctoamr file tool](filetools.md#nctoamr), as in
the [Pine Island Glacier example](pineisland.md).

## Initial geometry and boundary conditions

One of the most common data sources in ice sheet problems is the digital
elevation map (DEM). To load a DEM into BISICLES, determine the bedrock
topography and ice thickness, and create an appropriate hdf5 file from
those fields. In the [Pine Island Glacier example](pineisland.md), the
following lines are used to select an initial geometry loaded from a
file named  'pig-bisicles-1km.2d.hdf5 ', in which the relevant fields
are called  "thk " and  "topg "

    geometry.problem_type = LevelData
    inputLevelData.geometryFile = pig-bisicles-1km.2d.hdf5
    inputLevelData.thicknessName = thk
    inputLevelData.topographyName = topg

It is important to note the choice geometry.problem_type = LevelData
also implies the **lateral boundary conditions as well as initial
conditions**. By default, reflection boundary conditions (ice divides)
are imposed on all four domain edges. This can be changed to periodic
boundary conditions in the usual way, e.g

    amr.is_periodic = 0 1 0

selects reflection boundaries at the x-faces of the domain, but periodic
boundaries at the y-faces. To set a marine boundary condition at the
domain edge, use DomainEdgeCalvingModel, e.g

    CalvingModel.type = DomainEdgeCalvingModel
    CalvingModel.front_hi = 1 0    #impose a marine boundary at the high x-face 
    CalvingModel.front_lo = 0 0

On the other hand, if there is a calving front in the DEM, it can be
fixed, e.g it, use

    CalvingModel.type = FixedFrontCalvingModel
    CalvingModel.min_thickness = 1.0

as in the [Pine Island Glacier example](pineisland.md). This would
force any regions of the domain that begin ice-free to remain ice free,
and prevent the ice thickness from dropping below 1m in the rest of the
domain.

## Initial temperature

See also the description of [thermodynamics](thermodynamics.md)

Temperature is a three dimensional field which is discretized over a
fixed number N of layers in the current version of BISICLES. Each layer
comprises a multi-level 2D field of cell-centered values which is
located at the layer midpoint. The data hdf5 file needs to have
temperature data stored as N consecutive components, and the input file
should give the name of the first. BISICLES will then derive the
temperature of the top layer from that first component, and the
temperature of the remaining N-1 layers from the following N-1
components.

In the [Pine Island Glacier example](pineisland.md), the following
lines are relevant.

    amr.num_cells = 64 96 10
    amr.sigma = 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
    temperature.type = LevelData
    inputLevelData.temperatureFile = pig-bisicles-1km.2d.hdf5
    inputLevelData.temperatureName = temp000000

Note the 11 values of sigma; these are the values at the layer faces.
The temperature values are stored at the layer midpoint, so 10 are
needed, with the first called  "temp000000 ".

## Basal traction coefficient

See also the the description of [stresses](stress.md)

## Englacial stiffness coefficient,

See also the the description of [stresses](stress.md)

## Surface fluxes

See also the the description of [surface fluxes](surfaceflux.md)

## Inputs to the [inverse problem](##inverse)

See also the the description of [inverse
problem](velocity.md#inversevi)
