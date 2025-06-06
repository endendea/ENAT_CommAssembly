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

## merging rep-seqs.qza files from YXIN and AEST dataset
echo "Merging rep-seqs.qza files from YXIN and AEST dataset....."
qiime feature-table merge-seqs \
  --i-data /export/microlab/users/ENAT/YXIN_proj1_Q16712/QZA_input_files/XYIN_16S_515F_926R_07092023_Q16712_XYIN_representative_sequences.qza \
  --i-data /export/microlab/users/ENAT/AEST/results/AEST_rep-seqs.qza \
  --o-merged-data results/YXIN_AEST_merged-rep-seqs.qza

echo "Summarizing merged representative sequences....."
qiime feature-table tabulate-seqs \
  --i-data results/YXIN_AEST_merged-rep-seqs.qza \
  --o-visualization results/YXIN_AEST_merged-rep-seqs.qzv

echo "Generating new phylogeny tree....."
qiime phylogeny align-to-tree-mafft-fasttree \
          --i-sequences results/YXIN_AEST_merged-rep-seqs.qza \
          --o-alignment results/YXIN_AEST_merged-aligned-rep-seqs.qza \
          --o-masked-alignment results/YXIN_AEST_merged-masked-aligned-rep-seqs.qza \
          --o-tree results/YXIN_AEST_merged-unrooted-tree.qza \
          --o-rooted-tree results/YXIN_AEST_merged-rooted-tree.qza