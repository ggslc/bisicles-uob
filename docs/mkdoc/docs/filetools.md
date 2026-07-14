::: (#head)
-   [Index page](index.md)
-   [Top of page](#top)

# Contents

1.  [extract](#extract)
2.  [merge](#merge)
3.  [pythonf](#pythonf)
4.  [stats](#stats)
5.  [flatten](#flatten)
6.  [nctoamr](#nctoamr)
7.  [diagnostics](#diagnostics)
:::

::: (#main)
# BISICLES file tools

BISICLES includes some tools used to perform miscellaneous actions on
BISICLES/Chombo hdf5 files (amr files). They are located in the
code/filetools directory. The most commonly used are described below.

When built (e.g by typing \'make all\' in the code/filetools directory)
a number of executables, with names like e.g
extract2d.Linux.64.g++.gfortran.DEBUG.ex are created. Here, we will use
the base name (e.g) \'extract\' in examples; in practice, you would need
e.g extract2d.Linux.64.g++.gfortran.DEBUG.ex. All of the filetools will
provide a usage string when executed with no (or the incorrect number)
of command line arguments.

## [extract](#extract)

Extract one or more fields from one amr file and write to a new amr
file.

    usage: extract <input_file> <output_file> <var 1> [<var 2> [...]]

## [merge](#merge)

combine fields from two amr files but on the same mesh into one amr
file.

    usage: merge <input_file a> <input_file b> <output_file>

## [pythonf](#pythonf)

Evaluate a python cell by cell over an amr file and save the results to
a new amr file.

    usage: pythonf <input_file> <output_file> <python script> <python function> <input tuple> <output tuple>

pythonf can be used to carry out a wide range of post- and
pre-processing tasks. For example, imagine a file
\'linear_traction.2d.hdf5\' that contains fields C (a drag coefficient),
u (x-velocity), and v (y-velocity) that allowed the basal traction to be
computed for a linear viscous law (T_x = C \* u, T_y = C \* v). This
would be a typical result from the inverse problem. It is possible to
create an hdf5 designed to work with a third-power law ( T_i = D \*
\|u\|\^{-2/3} \* u_i) that should produce the same initial velocity by
setting D = C \* \|u\|\^2/3. Create a python file \'ctransform.py\'

    #ctransform.py
    def c_third(c_one,u,v,*etc):
        uu = (u*u + v*v)**(0.5)
        c_third = c_one * uu**(2.0/3.0)
        return c_third

then run

    $BISICLES_HOME/BISICLES/code/filetools/pythonf2d.Linux.64.g++.gfortran.DEBUG.ex  linear_traction.2d.hdf5 third_traction.2d.hdf5 ctransform c_third C,u,v D,u,v

## [stats](#stats)

Compute summary statistics (e.g volume) from a plot file. This tool has
been around for a long time, and we will not remove it, but you might
prefer its more powerful alternative: [diagnostics](#diagnostics)

    usage: stats <plot file> <ice_density> <water_density> <gravity> [mask_file] [mask_no_start = 0] [mask_no_end = mask_no_start] 

The stats tool is primarily used to compute quantities like ice volume
and volume above flotation. At minimum, supply a file name and
\<ice_density\> \<water_density\> \<gravity\> arguments. e.g, for the
[Pine Island Glacier](pineisland.md) example:

    ../../code/filetools/stats2d.Linux.64.g++.gfortran.DEBUG.ex plot.pigv5.1km.l1l2.2lev.000000.2d.hdf5 918 1028 9.81 | grep time
     time = 0.000000000000e+00 iceVolumeAll = 1.327650177413e+14  iceVolumeAbove = 6.305522125119e+13  groundedArea = 8.783100000000e+10  floatingArea = 5.154000000000e+09  totalArea = 9.298500000000e+10  groundedPlusOpenLandArea = 8.783100000000e+10  Total Melt = -1.995254865428e+11.

It is also possible to calculate the same statistics for a number of
sub-regions. For that, you need an hdf5 (e.g mask.2d.hdf5) file (created
with [nctoamr](#nctoamr), for example) including a single field (named
e.g \'mask\') that identifies each area by a unique double precision
representation of an integer. If the Pine Island Glacier domain was
split into sub-regions 0-3 inclusive, stats for regions 1,2 could be
computed by running

    ../../code/filetools/stats2d.Linux.64.g++.gfortran.DEBUG.ex plot.pigv5.1km.l1l2.2lev.000000.2d.hdf5 918 1028 9.81 mask.2d.hdf5 1 2

Output would include one line per sub-region.

## [flatten](#flatten)

Produce a uniform mesh hdf5 or netcdf file by conservative
averaging/interpolation data of from an amr file.

    usage: flatten <input_file> <output_file> level [x0 [y0 [z0]]]

The main benefit of BISICLES AMR is in the solution of PDEs to compute
the ice sheet state. The state itself is usually amenable to analysis at
more modest resolution than is required around the grounding line, and
on top of that the hdf5 amr files are a specialized format that many
find awkward. The flatten tool can be used to make a uniform mesh
version of the amr file, either as an hdf5 file that BISICLES can read
for other purposes, or a netcdf file for wider dissemination. All of the
interpolation is locally conservative.

When writing a NetCDF file, the location of the problem domain origin
can be specified as a command line argument. For example, the [Pine
Island Glacier](pineisland.md) example is set on a domain whose origin
corresponds to the point (-1707 km, -384 km ) on the usual Antarctic
Polar Stereographic Projection (EPSG 3031). To make a 4km resolution
file that can be directly compared with other data on that projection,
run (e.g)

    $BISICLES_HOME/BISICLES/code/filetools/flatten2d.Linux.64.g++.gfortran.DEBUG.ex plot.pigv5.1km.l1l2.2lev.000000.2d.hdf5 pig-epsg3031.nc 0 -1707000 -384000

For data at 1 km resolution, run

    $BISICLES_HOME/BISICLES/code/filetools/flatten2d.Linux.64.g++.gfortran.DEBUG.ex plot.pigv5.1km.l1l2.2lev.000000.2d.hdf5 pig-epsg3031.nc 2 -1707000 -384000

## [nctoamr](#nctoamr)

Convert a uniform mesh netcdf file into a uniform mesh hdf5 file.

    usage: nctoamr <input_file> <output_file> <var 1> [<var 2>, ...] 

If you have a NetCDF file called data.nc, which contains 2D fields named
thk and topg arranged on a regular grid (with equal spacing in x and y),
and must also contain a 1D field called x, you can run

    $BISICLES_HOME/BISICLES/code/filetools/nctoamr2d.Linux.64.g++.gfortran.DEBUG.ex data.nc data.2d.hdf5 thk topg

to obtain a file called data.2d.hdf5 that BISICLES can read as (say)
thickness and topography data.

## [diagnostics](#diagnostics)

Compute summary statistics including discharge (e.g volume) from a plot
file. Each statistic is a quantity integrated over one of several
regions the entire (sub)domain, the grounded ice, the floating ice, all
ice, and ice free regions. It is possible to supply a mask file, that
will allow the intergations to be carried out over each of several
subdomains (e.g drainage basins).

This tool is a more powerful alternative to [stats](#stats). It uses
keyword (name=value) arguments rather than positional arguments. It has
one mandatory argument (plot_file=???), the name of the plot file, e.g

      /path/to/diagnostics2d.Linux.64.g++.gfortran.DEBUG.OPT.ex plot_file=plot.camip_fullthule_expt34_tc2.063331.2d.hdf5

  name          default             purpose                                                          
  ------------- ------------------- ------------------------------------------------- -- -- -- -- -- --
  out_file      out_file=\"\"       specify output to a file rather than pout/stdio                  
  ice_density   ice_density=918.0   ice density (kg / m\^3)                                          

Output is to stdio / pout.X by defauly, and comprises multiple lines of
csv data. Each line includes one quantity, integrated over one of the
subdomains, at a given time. The data in each line are

-   the string csvheader or csvdata, to allow filtering
-   time : time in the plot file time unit (usually years)
-   maskNo : number of the masked sub-domain. 0 imples the whole domain
-   region: entire (sub)domain, grounded portion, floating portion etc
-   quantity: name of the quantity, e.g SMB (surface mass balance 1/2
    timestep previous), discharge (computed from ice thickness and
    velocity, flxDivReconstr (flux divergence computed from ice
    thickness and velocity), flxDivFile (flux divergence in the file, if
    present - will differ from fluxDivReconstr because it is computed
    1/2 timestep previous).
-   unit: unit of measurement
-   value: e.g SMB inegrated over the region

```{=md}
<!-- -->
```
      csvheader,time,maskNo,region,quantity,unit,value
      ...
      csvdata,1.451000000000e+03,0,grounded,SMB,m3/a,2.661530348348e+11
      csvdata,1.451000000000e+03,0,grounded,BMB,m3/a,0.000000000000e+00
      csvdata,1.451000000000e+03,0,grounded,dhdt,m3/a,-4.574041539644e+10
      csvdata,1.451000000000e+03,0,grounded,flxDivFile,m3/a,3.118848162762e+11
      csvdata,1.451000000000e+03,0,grounded,calving,m3/a,8.633955017224e+06
      csvdata,1.451000000000e+03,0,grounded,flxDivReconstr,m3/a,3.121113211384e+11
      csvdata,1.451000000000e+03,0,grounded,discharge,m3/a,3.118966365235e+11
      ...
:::
