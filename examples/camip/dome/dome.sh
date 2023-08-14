for TAGCAP in -1 0 1 2 3
do
    sed -e s/@TAGCAP/$TAGCAP/ inputs.dome_spin.template > inputs.dome_spin.tc$TAGCAP
done

export PYTHONPATH=$PWD/../python/:$PYTHONPATH
for TAGCAP in -1 0 1 2 3
do
    mpirun -np 40 ./driver2d.Linux.64.mpiCC.gfortran.DEBUG.OPT.MPI.ex inputs.dome_spin.tc$TAGCAP
done
