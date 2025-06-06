#!/bin/bash
#SBATCH --output=16S_Q2_%j.out
#SBATCH --job-name=q2_16S
#SBATCH --partition=Bytesnet
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16


# go to the directory
cd /export/microlab/users/ENAT/AEST
project="AEST"

## On Linux
echo "Activating QIIME2 amplicon deployment version 2024.10 amplicon"
source /export/microlab/miniconda3/etc/profile.d/conda.sh
conda activate /export/microlab/miniconda3/envs/qiime2-amplicon-2024.5

# importing data
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path raw_data/ \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza

# Visualize sample-metadata
qiime metadata tabulate \
    --m-input-file QZA_input_files/${project}@metadata_formatted.txt \
    --o-visualization ${project}@metadata_formatted.qzv

# Summarize the feature-table
qiime feature-table summarize-plus \
    --i-table QZA_input_files/${project}_table.qza \
    --m-metadata-file QZA_input_files/${project}@metadata_formatted.txt \
    --o-summary ${project}_table.qzv \
    --o-sample-frequencies ${project}_sample_frequencies.qza \
    --o-feature-frequencies ${project}_feature_frequencies.qza

# Summarize the feature-data
qiime feature-table tabulate-seqs \
    --i-data QZA_input_files/${project}_representative_seqs.qza \
    --o-visualization ${project}_representative_seqs.qzv

# Filter out seq/feature that are only found in  1 samples
# first, filter the feature table
qiime feature-table filter-features \
    --i-table QZA_input_files/${project}_table.qza \
    --p-min-samples 2 \
    --o-filtered-table ${project}_table_ms2.qza

# then, filter the feature data
qiime feature-table filter-seqs \
    --i-data QZA_input_files/${project}_representative_sequences.qza \
    --i-table ${project}_table_ms2.qza \
    --o-filtered-data ${project}_representative_sequences_ms2.qza

# summarize the new feature table
qiime feature-table summarize-plus \
    --i-table ${project}_table_ms2.qza \
    --m-metadata-file QZA_input_files/${project}@metadata_formatted.txt \
    --o-summary ${project}_table_ms2.qzv \
    --o-sample-frequencies ${project}_sample_frequencies_ms2.qza \
    --o-feature-frequencies ${project}_feature_frequencies_ms2.qza

# apply classifier to the dataset
qiime feature-classifier classify-sklearn \
    --i-classifier QZA_input_files/${project}_silva-138-99-nb-classifier.qza \
    --i-reads ${project}_representative_sequences_ms2.qza \
    --o-classification ${project}_taxonomy.qza

# calculating alpha and beta diversity metrics
# generate kmer feature table
qiime kmerizer seqs-to-kmers \
    --i-sequences ${project}_representative_sequences_ms2.qza \
    --i-table ${project}_table_ms2.qza \
    --o-kmer-table ${project}_kmer_table.qza

# summarize kmer table
qiime feature-table summarize-plus \
    --i-table ${project}_kmer-table.qza \
    --m-metadata-file QZA_input_files/${project}@metadata_formatted.txt \
    --o-summary ${project}_kmer-table.qzv \
    --o-sample-frequencies ${project}_kmer-sample-frequencies.qza \
    --o-feature-frequencies ${project}_kmer-feature-frequencies.qza

# Diversity metrics are sensitive to differences uneven sequencing effort (seq depth).
# Here, we perform rarefaction and bootstrap before calculate alpha and beta diversity
# Rarefaction to determine the depth;
# I choose around median of feature per sample (21x10^6) as seq depth
qiime diversity alpha-rarefaction \
    --i-table ${project}_kmer_table.qza \
    --p-max-depth 21000000 \
    --m-metadata-file QZA_input_files/${project}@metadata_formatted.txt \
    --o-visualization ${project}_alpha_rarefaction.qzv

# Computing bootstrapped alpha and beta diversity metrics
# I choose 14000000 to include all the samples
qiime boots core-metrics --i-table ${project}_kmer_table.qza \
    --m-metadata-file QZA_input_files/${project}@metadata_formatted.txt \
    --p-sampling-depth 14000000 \
    --p-n 10 \
    --p-replacement \
    --p-alpha-average-method median \
    --p-beta-average-method medoid \
    --output-dir ${project}_boots-core-metrics

# creating PCoAs
# unweighted-PCoAs based on Jaccard metrics
qiime vizard scatterplot-2d \
    --m-metadata-file QZA_input_files/${project}@metadata_formatted.txt ${project}_boots-core-metrics/pcoas/jaccard.qza ${project}_boots-core-metrics/alpha_diversities/shannon.qza \
    --o-visualization jaccard-diversity-scatterplot.qzv

# Weighted-PCoAs based on Braycurtis metrics
qiime vizard scatterplot-2d \
    --m-metadata-file QZA_input_files/${project}@metadata_formatted.txt ${project}_boots-core-metrics/pcoas/braycurtis.qza ${project}_boots-core-metrics/alpha_diversities/shannon.qza \
    --o-visualization braycurtis-diversity-scatterplot.qzv

# Differential abundance testing with ANCOM-BC
# filter sample to exclude control
qiime feature-table filter-samples \
    --i-table ${project}_table_ms2.qza \
    --m-metadata-file QZA_input_files/${project}@metadata_formatted.txt \
    --p-where '[sludge_source] IN ("Almere", "Bath")' \
    --o-filtered-table ${project}_table_ms2_sludge_source.qza

# Collapse the representative seqs (ASVs) into genera
qiime taxa collapse \
    --i-table ${project}_table_ms2_sludge_source.qza \
    --i-taxonomy ${project}_taxonomy.qza \
    --p-level 6 \
    --o-collapsed-table ${project}_genus_table_ms2_sludge_source.qza

# Apply differential abundance testing
qiime composition ancombc \
    --i-table ${project}_genus_table_ms2_sludge_source.qza \
    --m-metadata-file QZA_input_files/${project}@metadata_formatted.txt \
    --p-formula sludge_source \
    --p-reference-levels 'sludge_source::Almere' \
    --o-differentials ${project}_genus_ancombc.qza

# Visualize the abundance testing result
qiime composition da-barplot \
    --i-data ${project}_genus_ancombc.qza \
    --p-significance-threshold 0.001 \
    --p-level-delimiter ';' \
    --o-visualization ${project}_genus_ancombc.qzv

# Next, we try to do abundance testing based on substrate
qiime composition ancombc \
    --i-table ${project}_genus_table_ms2_sludge_source.qza \
    --m-metadata-file QZA_input_files/${project}@metadata_formatted.txt \
    --p-formula substrate \
    --p-reference-levels 'substrate::acetate' \
    --o-differentials ${project}_genus_ancombc_substrate.qza

qiime composition da-barplot \
    --i-data ${project}_genus_ancombc_substrate.qza \
    --p-significance-threshold 0.001 \
    --p-level-delimiter ';' \
    --o-visualization ${project}_genus_ancombc_substrate.qzv