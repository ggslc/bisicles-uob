for TAGCAP in -1 0 1 2 3
do
    sed -e s/@TAGCAP/$TAGCAP/ inputs.fullthule_spin.template > inputs.fullthule_spin_tc$TAGCAP
    sed -e s/@TAGCAP/$TAGCAP/ job_fts_bp > job_fts_tc$TAGCAP.sh
done

for TAGCAP in -1 0 1 2 3
do
    sed -e s/@TAGCAP/$TAGCAP/ inputs.fullthule_expt34.template > inputs.fullthule_expt34_tc$TAGCAP
    LAST=`ls -th chk.camip_fullthule_spin_tc$TAGCAP.??????.2d.hdf5 | head -n 1`
    if [ "$LAST" != "" ]; then
	echo $LAST
	echo "\n" >> inputs.fullthule_expt34_tc$TAGCAP
	echo "amr.restart_file = $LAST" >> inputs.fullthule_expt34_tc$TAGCAP
	echo "\n" >> inputs.fullthule_expt34_tc$TAGCAP
    fi
    sed -e s/@TAGCAP/$TAGCAP/ job_fte34_bp > job_fte34_tc$TAGCAP.sh


done
