#!/bin/bash
#SBATCH --output=logs/16S_Q2_%j.out
#SBATCH --error=logs/16S_Q2_%j.err
#SBATCH --job-name=q2_16S
#SBATCH --partition=Bytesnet
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16


set -euo pipefail

# go to the directory
cd /export/microlab/users/ENAT/AEST

## Load conda env
echo "Activating conda environment for QIIME2...."
source /export/microlab/miniconda3/etc/profile.d/conda.sh
conda activate /export/microlab/miniconda3/envs/gut-to-soil

# Run Snakemake
echo "Starting snakemake for QIIME2 analysis..."
snakemake --cores 16 --use-conda --rerun-incomplete --printshellcmds --configfile config.yaml