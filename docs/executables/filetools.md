# BISICLES File Tools

BISICLES includes several tools used to perform miscellaneous actions on BISICLES/Chombo hdf5 files (amr files). These tools are located in the `code/filetools` directory.

When built (e.g., by typing `make all` in the `code/filetools` directory), executables with names like `extract2d.Linux.64.g++.gfortran.DEBUG.ex` are created. In the examples below, we use the base name (e.g., `extract`) for simplicity; in practice, you would need the full executable name.

All filetools provide a usage string when executed with no or incorrect command line arguments.

## extract

Extract one or more fields from one amr file and write to a new amr file.

```
usage: extract <input_file> <output_file> <var 1> [<var 2> [...]]
```

## merge

Combine fields from two amr files (on the same mesh) into one amr file.

```
usage: merge <input_file a> <input_file b> <output_file>
```

## pythonf

Evaluate a Python function cell-by-cell over an amr file and save the results to a new amr file.

```
usage: pythonf <input_file> <output_file> <python script> <python function> <input tuple> <output tuple>
```

### Example: Converting Drag Coefficient Laws

The `pythonf` tool can carry out a wide range of post- and pre-processing tasks. For example, imagine a file `linear_traction.2d.hdf5` that contains fields:
- `C` (a drag coefficient)
- `u` (x-velocity)
- `v` (y-velocity)

For a linear viscous law (T_x = C * u, T_y = C * v), you can convert this to a third-power law where (T_i = D * |u|^{-2/3} * u_i) by setting D = C * |u|^{2/3}.

Create a Python file `ctransform.py`:

```python
def c_third(c_one, u, v, *etc):
    uu = (u*u + v*v)**(0.5)
    c_third = c_one * uu**(2.0/3.0)
    return c_third
```

Then run:

```bash
$BISICLES_HOME/BISICLES/code/filetools/pythonf2d.Linux.64.g++.gfortran.DEBUG.ex \
  linear_traction.2d.hdf5 third_traction.2d.hdf5 ctransform c_third C,u,v D,u,v
```

## stats

Compute summary statistics (e.g., volume) from a plot file.

```
usage: stats <plot file> <ice_density> <water_density> <gravity> [mask_file] [mask_no_start = 0] [mask_no_end = mask_no_start]
```

**Note:** This tool has been around for a long time and will not be removed, but you might prefer its more powerful alternative: [diagnostics](#diagnostics).

The stats tool is primarily used to compute quantities like ice volume and volume above flotation. At minimum, supply a file name and the three density/gravity arguments. For example, for the Pine Island Glacier example:

```bash
../../code/filetools/stats2d.Linux.64.g++.gfortran.DEBUG.ex \
  plot.pigv5.1km.l1l2.2lev.000000.2d.hdf5 918 1028 9.81 | grep time
```

Output:
```
time = 0.000000000000e+00 iceVolumeAll = 1.327650177413e+14  iceVolumeAbove = 6.305522125119e+13  groundedArea = 8.783100000000e+10  floatingArea = 5.154000000000e+09  totalArea = 9.298500000000e+10
```

### Using Masks for Sub-regions

You can calculate statistics for sub-regions using a mask file. The mask file (created with [nctoamr](#nctoamr), for example) should contain a single field that identifies each area by a unique double-precision integer. 

For example, if the Pine Island Glacier domain is split into sub-regions 0-3, compute stats for regions 1 and 2:

```bash
../../code/filetools/stats2d.Linux.64.g++.gfortran.DEBUG.ex \
  plot.pigv5.1km.l1l2.2lev.000000.2d.hdf5 918 1028 9.81 mask.2d.hdf5 1 2
```

Output will include one line per sub-region.

## flatten

Produce a uniform mesh hdf5 or netcdf file by conservative averaging/interpolation of data from an amr file.

```
usage: flatten <input_file> <output_file> level [x0 [y0 [z0]]]
```

The main benefit of BISICLES AMR is in solving PDEs to compute the ice sheet state. The state itself is usually amenable to analysis at more modest resolution than required around the grounding line. The hdf5 amr files are also a specialized format that many find awkward. The `flatten` tool can convert an amr file to either:
- An hdf5 file that BISICLES can read for other purposes
- A netcdf file for wider dissemination

All interpolation is locally conservative.

### NetCDF Output with Projections

When writing a NetCDF file, you can specify the location of the problem domain origin as command-line arguments. For example, the Pine Island Glacier example is set on a domain whose origin corresponds to (-1707 km, -384 km) on the Antarctic Polar Stereographic Projection (EPSG 3031).

To make a 4km resolution file for direct comparison with other EPSG 3031 data:

```bash
$BISICLES_HOME/BISICLES/code/filetools/flatten2d.Linux.64.g++.gfortran.DEBUG.ex \
  plot.pigv5.1km.l1l2.2lev.000000.2d.hdf5 pig-epsg3031.nc 0 -1707000 -384000
```

For 1 km resolution:

```bash
$BISICLES_HOME/BISICLES/code/filetools/flatten2d.Linux.64.g++.gfortran.DEBUG.ex \
  plot.pigv5.1km.l1l2.2lev.000000.2d.hdf5 pig-epsg3031.nc 2 -1707000 -384000
```

## nctoamr

Convert a uniform mesh netcdf file into a uniform mesh hdf5 file.

```
usage: nctoamr <input_file> <output_file> <var 1> [<var 2> [...]]
```

If you have a NetCDF file called `data.nc` containing 2D fields named `thk` and `topg` arranged on a regular grid (with equal spacing in x and y), and also containing a 1D field called `x`, you can run:

```bash
$BISICLES_HOME/BISICLES/code/filetools/nctoamr2d.Linux.64.g++.gfortran.DEBUG.ex data.nc data.2d.hdf5 thk topg
```

This produces a file called `data.2d.hdf5` that BISICLES can read as thickness and topography data.

## diagnostics

Compute summary statistics including discharge from a plot file. Each statistic is a quantity integrated over one of several regions: the entire (sub)domain, grounded ice, floating ice, all ice, and ice-free regions. You can optionally supply a mask file to allow integrations over subdomains (e.g., drainage basins).

This tool is a more powerful alternative to [stats](#stats). It uses keyword (name=value) arguments rather than positional arguments.

Mandatory argument:
- `plot_file=???` - the name of the plot file

Example:

```bash
/path/to/diagnostics2d.Linux.64.g++.gfortran.DEBUG.OPT.ex \
  plot_file=plot.camip_fullthule_expt34_tc2.063331.2d.hdf5
```

### Output Format

Output is to stdio/pout.X by default and comprises multiple lines of CSV data. Each line includes one quantity integrated over one of the subdomains at a given time. The data in each line are:

- `csvheader` or `csvdata` - allows filtering
- `time` - time in the plot file time unit (usually years)
- `maskNo` - number of the masked sub-domain (0 implies the whole domain)
- `region` - entire (sub)domain, grounded portion, floating portion, etc.
- `quantity` - name of the quantity (e.g., SMB for surface mass balance, discharge, flxDivReconstr)
- `unit` - unit of measurement
- `value` - e.g., SMB integrated over the region

Example output:

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
```
