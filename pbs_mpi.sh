#!/bin/bash
### Job Name
#PBS -N mpi_test

### Project code
#PBS -q intel
#PBS -j oe
#PBS -l select=6:ncpus=2
#PBS -o logs

module load openmpi3/gcc/64/3.1.4
module load pbs
module load julia/1.3.1

cd $PBS_O_WORKDIR

export OMP_NUM_THREADS=2
export UCX_WARN_UNUSED_ENV_VARS=n


# MPI will automatically use the PBS_NODEFILE to spread the jobs

mpiexec -np 6 julia bcast_load_genomes.jl
