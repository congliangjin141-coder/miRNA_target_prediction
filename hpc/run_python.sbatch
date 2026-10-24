#!/bin/bash
# Run one script of the pipeline on a compute node.
# Usage, from the repository root:
#   mkdir -p logs
#   sbatch hpc/run_python.sbatch analysis/02_build_ref_score_table.py
# Adjust the conda environment name and the resources to your cluster.
#SBATCH --job-name=mirna_bench
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

SCRIPT="$1"
if [ -z "$SCRIPT" ]; then echo "Usage: sbatch hpc/run_python.sbatch <path/to/script>.py"; exit 1; fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mirna
export MPLBACKEND=Agg

python -u "$SCRIPT"
