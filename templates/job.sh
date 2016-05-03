#!/bin/bash
#PBS -l walltime=<<WALLTIME>>
#PBS -l nodes=1:ppn=12
#PBS -l vmem=60gb
#PBS -m ae
#PBS -j oe 

NTHREADS=12
export OMP_NUM_THREADS=$NTHREADS

cd $PBS_O_WORKDIR

mpirun -np 12 NCSD
