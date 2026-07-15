# BISICLES common issues

## The code does not compile

## The code compiles but does not run

1.  **ImportError: No module named xyz**. For some reason the python
    interpreter was unable to read xyz.py, e.g

        ImportError: No module named rift
        failed to import Python module rift
        failed to import Python module  !!!
          

    In this case, either rift.py does not exists, or, it is in a
    directory outside  $PYTHONPATH. By default  $PYTHONPATH does not
    include the current working directory, so this is an easy mistake to
    make. Try e.g

            export PYTHONPATH=$PWD:$PYTHONPATH
          

## The code compiles and runs but crashes or fails in some spectacular way

Sometimes the code starts, but then terminates, leaving a more or less
cryptic error in one other other log.

1.  **maxFaceVelocity  > 0.5  * HUGE_VEL**. A very common error. It
    implies that the velocity solver - which can fail for many reasons,
    has produced an enormous speed somewhere. Start off by looking at
    the [velocity calculation](velocity.md), but be aware that the
    velocity solver might be struggling because of its inputs -
    especially the [geometry](geometry.md) and the [basal
    traction](stress.md#basal).

## The code compiles and runs to completions but the results are drivel

Arguably, all results are like this :). But there are some things to
watch out for
