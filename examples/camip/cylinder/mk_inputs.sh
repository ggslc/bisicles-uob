
for CFL in 0.25 0.125 0.0625 0.03125
do
    for CR in 0 250 500 750
    do
	for LEV in 0 1 2 3 4
	do
	    NAME="cyl_"$LEV"lev_cr"$CR"_cfl"$CFL
	    sed -e s/@CR/$CR/ -e s/@CFL/$CFL/ -e s/@LEV/$LEV/ inputs.cyl_template > inputs.$NAME
	done
    done
done
