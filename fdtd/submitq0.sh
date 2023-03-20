#!/bin/bash
#SBATCH -A p31865               # Allocation
#SBATCH -p normal                # Queue
#SBATCH -t 24:00:00             # Walltime/duration of the job
#SBATCH -N 8                    # Number of Nodes
#SBATCH --ntasks-per-node=28
#SBATCH --mem=0               # Memory per node in GB needed for a job. Also see --mem-per-cpu
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ysu2024@u.northwestern.edu

# unload any modules that carried over from your command line session
module purge

# add a project directory to your PATH (if needed)

# load modules you need to use
module load mpi

# A command you actually want to execute:

cfgfile=chromatin_on_glass_NAi55_a20_dx10.cfg

mpirun /projects/p31209/angora/source/angora-0.22.5/angora $cfgfile
