for EXPT in rp ri rn 
do
    for TC in -1 0 1 2 3;
    do
	sed -e s/@EXPT/$EXPT/ -e s/@TC/$TC/ inputs.channel_frac.template > inputs.channel_frac.$EXPT.tc$TC
    done
done
