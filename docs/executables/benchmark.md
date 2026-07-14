# BISICLES solver benchmark executable

The `solverBenchmark` executable is used to run the ice velocity solver as a standalone tool to examine performance and convergence.

## Background

The vast majority of CPU time in a typical BISICLES simulation is spent in the velocity solver. Running the solver as a standalone tool allows you to:
- Examine solver performance characteristics
- Debug convergence issues
- Analyze solver behavior independently

## Usage

### Step 1: Enable solver output in the driver run

In the configuration file for the driver run, add:

```
amr.write_solver_rhs = 1
```

This writes the solver problem setup to an HDF5 file.

### Step 2: Configure the benchmark executable

In the configuration file for the benchmark run, specify the HDF5 file from Step 1:

```
main.filename = <hdf5 filename>
```

This tells the benchmark code where to find the problem setup.

## Example: Running the regression test

To see how this works (assuming Chombo and dependencies are already set up):

```bash
cd BISICLES/code/benchmark
make regression
```

This performs the following steps:
1. Compiles the driver code in `BISICLES/exec2D` if necessary
2. Runs the driver with `BISICLES/benchmark/inputs.regression` as the configuration file, producing a plot file named `stream.regression.preSolve.000000.2d.hdf5`
3. Runs the `solverBenchmark` executable with the same configuration file
4. The benchmark reads the configuration file and HDF5 file, then performs a velocity solve that should be identical to the one performed by the driver

This provides a useful way to isolate and analyze the solver's behavior.
