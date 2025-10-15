getcre()
{
  case $smod in
    l1l2) cre=L1L2;;
    ssa)  cre=GlensLaw;;
    *) echo "unknown stress model"
  esac
}


#forward problems
for smod in l1l2 ssa;
  do
  for lev in 0 1 2 3 4 5; 
    do
    tagcap=$(( lev - 1 )) 
    getcre
    
    sed  -e s/@MAXLEVEL/$lev/ -e s/@TAGCAP/$tagcap/ -e s/@SMOD/$smod/ -e s/@CRE/$cre/ inputs.pigv5.1km.template > inputs.pigv5.1km.$smod.l$lev    
  done
done

#inverse problems
for BASALBETA in m1 m3 m3_u100;
do
    case $BASALBETA in
	m1) BASALM=1.0; BASALU0=-1.0;;
	m3) BASALM=0.3333; BASALU0=-1.0;;
	m3_u100) BASALM=0.3333; BASALU0=100.0;;
	*) echo "unknown basal friction model"
    esac
    sed  -e s/@BASALBETA/$BASALBETA/ -e s/@BASALM/$BASALM/ -e s/@BASALU0/$BASALU0/ inputs.pig.ctrl > inputs.pig.ctrl$BASALBETA 
done
