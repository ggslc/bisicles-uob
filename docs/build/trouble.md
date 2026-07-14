# BISICLES common issues

This document describes common issues encountered when building and running BISICLES, with troubleshooting guidance.

## The code does not compile

See the [build instructions](../readme.html) for system requirements and dependency setup.

## The code compiles but does not run

### ImportError: No module named xyz

If the Python interpreter cannot find a module:

```
ImportError: No module named rift
failed to import Python module rift
failed to import Python module  !!!
```

This means either:
- The Python file does not exist, or
- It is in a directory outside `$PYTHONPATH`

By default, `$PYTHONPATH` does not include the current working directory, which is a common mistake.

**Solution:** Add the current directory to your Python path:

```bash
export PYTHONPATH=$PWD:$PYTHONPATH
```

## The code compiles and runs but crashes

Sometimes the code starts but then terminates with a cryptic error in the log files.

### maxFaceVelocity > 0.5 * HUGE_VEL

This is a very common error indicating that the velocity solver has produced an unreasonably large speed somewhere in the domain.

**Troubleshooting steps:**
1. Review the [velocity calculation](../velocity.html) documentation
2. Check your input data, especially:
   - [Geometry](../geometry.html) (thickness, topography, boundary conditions)
   - [Basal traction](../stress.html#basal) coefficient

The velocity solver can fail for many reasons, often due to problematic inputs rather than the solver itself.

## The code compiles and runs to completion but the results are wrong

If your simulation completes but produces unexpected results, there are several things to check. See the individual configuration documentation pages for guidance on setting up physically realistic problems.
