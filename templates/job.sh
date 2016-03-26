#!/bin/bash
#PBS -l walltime=<<WALLTIME>>
#PBS -l nodes=1:ppn=${NTHREADS}
#PBS -l vmem=60gb
#PBS -m ae
#PBS -M dilyn.fullerton@alumni.ubc.ca
#PBS -j oe 

NTHREADS=12
export OMP_NUM_THREADS=$NUM_THREADS

cd $PBS_O_WORKDIR

mpirun -np 12 NCSD
