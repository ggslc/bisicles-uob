for TAGCAP in -1 0 1 2 3
do
    sed -e s/@TAGCAP/$TAGCAP/ inputs.fulldome_spin.template > inputs.fulldome_spin_tc$TAGCAP
    sed -e s/@TAGCAP/$TAGCAP/ job_fds_bp > job_fds_tc$TAGCAP.sh
done

