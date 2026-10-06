#!/bin/bash
#SBATCH --account=andrea.perin
#SBATCH --job-name=1024_cleanup
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=andrea.perin@irz.uni-hannover.de
#SBATCH --mail-type=ALL

julia --project -e 'using Pkg; Pkg.instantiate(); include("./pce_internal_scheduler.jl")'
