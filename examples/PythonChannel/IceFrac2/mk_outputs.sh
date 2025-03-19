for EXPT in rp ri rn ap ai an
do
    for TC in -1 0 1 2 3;
    do
	./driver2d.Linux.64.g++.gfortran.DEBUG.OPT.ex inputs.channel_frac.$EXPT.tc$TC > sout.channel_frac.$EXPT.tc$TC &
    done
done
