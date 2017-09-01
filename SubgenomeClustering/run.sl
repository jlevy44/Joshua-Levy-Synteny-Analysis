#!/bin/bash -l
#SBATCH --nodes=2
#SBATCH --time=03:30:00
#SBATCH --partition=regular
#SBATCH --license=SCRATCH
#SBATCH --constraint=haswell
#SBATCH --workdir=.
srun ./runCluster.sh