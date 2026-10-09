#!/bin/bash
#SBATCH --account=andrea.perin
#SBATCH --job-name=1024_cleanup
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=andrea.perin@irz.uni-hannover.de
#SBATCH --mail-type=ALL

# Submit from the repository root: sbatch scripts/surrogates/pce.sh
export THERMOOPTIPLAN_DATA=/work/andrea.perin/ThermoOptiPlan/data

julia --project=. -e 'using Pkg; Pkg.instantiate(); include("scripts/surrogates/pce.jl")'
