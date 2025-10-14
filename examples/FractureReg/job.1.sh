#!/bin/bash 
#
# job script for isambard 3
#
#SBATCH -J ideal_shelf_1
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=18
module load craype-network-ofi
module load PrgEnv-gnu
module load cray-python/3.11.7

NCORE=$((18 * $SLURM_JOB_NUM_NODES))
echo "n_node: $SLURM_JOB_NUM_NODES , n_core: $NCORE"
DRIVER=$PWD/driver2d.Linux.64.CC.ftn.DEBUG.OPT.MPI.PETSC.GNU.ex
NAME=1_generate_data
RUNDIR=$SLURM_SUBMIT_DIR/$SLURM_JOB_ID-$NAME
mkdir -p $RUNDIR
INFILEBASE=inputs.$NAME
INFILE=$INFILEBASE"."$SLURM_JOB_ID
cp $SLURM_SUBMIT_DIR/$INFILEBASE $RUNDIR/$INFILE
cp $SLURM_SUBMIT_DIR/.petscrc $RUNDIR/
cd $RUNDIR

export CH_TIMER=1
export CH_OUTPUT_INTERVAL=999
export PYTHONPATH=$PWD:$SLURM_SUBMIT_DIR:$PYTHONPATH

echo "srun $DRIVER $INFILE"
srun $DRIVER $INFILE
