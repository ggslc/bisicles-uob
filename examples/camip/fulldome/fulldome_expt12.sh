INPUTS_BASE=inputs.fulldome_expt12
TEMPLATE=$INPUTS_BASE.template
JOB_TEMPLATE=job_fde12_bp

for TAGCAP in -1 0 1 2 3
do
    
    INPUTS= $INPUTS_BASE"_tc"$TAGCAP
    JOB=$JOB_TEMPLATE"_tc"$TAGCAP".sh"
    sed -e s/@TAGCAP/$TAGCAP/ $TEMPLATE > $INPUTS
    SPIN=`ls -th chk.camip_fulldome_spin_tc$TAGCAP.??????.2d.hdf5 | head -n 1`
    echo "\n" >> $INPUTS
    echo "amr.restart_file = $LAST" >> $INPUTS
    echo "amr.restart_time = 0.0" >> $INPUTS
    echo "amr.restart_set_time = true"  >> $INPUTS
    echo "\n" >> $INPUTS

    sed -e s/@TAGCAP/$TAGCAP/ $JOB_TEMPLATE > $JOB
done

