#!/bin/bash
#SBATCH --ntasks=1 
#SBATCH -p compsci
#SBATCH --exclude=linux[1-30]
#SBATCH -o slurm-%A_%a.out -e slurm-%A_%a.err
#SBATCH --cpus-per-task=40
#SBATCH --mem=150000
#SBATCH --time="7-0"

# NOTE: Currently the array job listing is taken care of in test_wrapper.py
# 	alternative sbatch array line: --array=1-5
# 	alternative sbatch array line: --array=1-50%5 

script_loc=/usr/project/dlab/Users/gth/projects/osprey_test_suite/scripts/run_design.py
#script_loc=/usr/project/dlab/Users/gth/projects/osprey_test_suite/scripts/run_design_single_pfunc.py

# try to set JAVA_HOME correctly
if [ ! -x "$JAVA_HOME/bin/java" ]; then
	export JAVA_HOME="/usr/lib/jvm/java-8-openjdk-amd64"
fi

#Add parameters
params=()
params+=(-e 0.68)
params+=(--slurm-out="slurm-"$SLURM_ARRAY_JOB_ID"_"$SLURM_ARRAY_TASK_ID".out")
params+=(--slurm-err="slurm-"$SLURM_ARRAY_JOB_ID"_"$SLURM_ARRAY_TASK_ID".err")
# Determine which algorithm to use
params+=(--algo="$2")

# activate venv
source /usr/project/dlab/Users/gth/code/osprey/virtual_envs/sharkstar_parallelism/bin/activate

declare -a lines
readarray -t lines < "$1"
python "$script_loc" -c 40 "${lines[$((SLURM_ARRAY_TASK_ID-1))]}" "output_$SLURM_ARRAY_TASK_ID.json" "${params[@]}"
