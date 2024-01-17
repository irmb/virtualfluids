#!/bin/bash -l

# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder

#SBATCH --partition=standard
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --job-name=VirtualFluids
#SBATCH --ntasks-per-node=4
#SBATCH --output=virtualfluids-mpich-bind/virtualfluids.out

date

echo "SLURM_JOBID="$SLURM_JOBID
echo "SLURM_JOB_NODELIST"=$SLURM_JOB_NODELIST
echo "SLURM_NNODES"=$SLURM_NNODES
echo "SLURM_TASKS_PER_NODE"=$SLURM_TASKS_PER_NODE
echo "working directory = "$SLURM_SUBMIT_DIR
echo ""

module purge
module load singularity/3.9.9
module load mpi/mpich/mpich_3.2
module list

cd virtualfluids-mpich-bind/
mkdir -p results

export MPI_DIR="/cluster/mpi/mpich"
srun --mpi=pmi2 singularity exec --bind "$MPI_DIR" rockylinux9-mpich-bind.sif /build/bin/FlowAroundCylinder cylinder.cfg 
