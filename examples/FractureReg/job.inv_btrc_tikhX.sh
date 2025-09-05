#!/bin/bash 
#
# job script for isambard 3
# job array for single parameter (e.g L-curve) exploration
#
#SBATCH -J ideal_shelf_inv_btrc_tikh

#SBATCH --time=16:00:00

#wholde
#SBATCH --nodes=1
#SBATCH --tasks-per-node=8
#18 jobs in the array, 8 cores each
#SBATCH --array=0-17

module load craype-network-ofi
module load PrgEnv-gnu
module load cray-python/3.11.7

NCORE=$((18 * $SLURM_JOB_NUM_NODES))
echo "n_node: $SLURM_JOB_NUM_NODES , n_core: $NCORE"
DRIVER=$PWD/driver2d.Linux.64.CC.ftn.DEBUG.OPT.MPI.PETSC.GNU.ex
NAME=ideal_shelf_inv_btrc_tikh
RUNDIR=$SLURM_SUBMIT_DIR/$SLURM_ARRAY_JOB_ID"-"$NAME/$SLURM_ARRAY_TASK_ID
mkdir -p $RUNDIR
INFILEBASE=inputs.$NAME.X
INFILE=inputs.$NAME.$SLURM_ARRAY_JOB_ID"-"$SLURM_ARRAY_TASK_ID

#work out parameter from $SLURM_ARRAY_TASK_ID
PARM=TIKHC
EXPO=$(( SLURM_ARRAY_TASK_ID - 9))
MANT=1
VAL=$MANT".0e"$EXPO
echo $SLURM_ARRAY_JOB_ID" "$SLURM_ARRAY_TASK_ID" "$PARM"="$VAL
sed -e s/@"$PARM"/"$VAL"/ $SLURM_SUBMIT_DIR/$INFILEBASE > $RUNDIR/$INFILE
cp $SLURM_SUBMIT_DIR/.petscrc $RUNDIR/
cd $RUNDIR
ln -s $SLURM_SUBMIT_DIR/synthetic_observations_500m.2d.hdf5

export CH_TIMER=1
export CH_OUTPUT_INTERVAL=999
export PYTHONPATH=$PWD:$SLURM_SUBMIT_DIR:$PYTHONPATH

echo "srun $DRIVER $INFILE"
srun $DRIVER $INFILE
