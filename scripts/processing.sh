#!/bin/bash

###########################################################################
# METABARCODING PIPELINE BASED ON QIIME2's DADA2 workflow
#
# This file runs the commands to process the raw fastq files.
#
# It needs the conda environments :
#   - bioinf_util_env with cutadapt (4.8), fastqc (0.12.1) and multiqc (1.21)
#   - qiime2_env with the workflow (2024.2, downloaded from QIIME2 website)
#
# qiime2_env is from https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.2-py38-linux-conda.yml
# Actual conda environment specifications are in conda_environments/
# 
# It assumes:
#   - the raw data files are in raw_data/
#   - the SILVA database is in databases/qiime2-silva138_ssuref_nr99_full/
#   - this command file is run in the working directory with these two folders as subfolders
#
# The database is 138 SSURef NR99 full-length (sequences and taxonomy qza files)
# and downloaded from the QIIM2 website
# 
# It will create the output directories:
#   - clean_data/ with processed fastq files and fastqc statistics
#   - qiime2_pip/ with the intermediary files from QIIME2
# 
# author: Michaël Pierrelée
# date: April 2024
###########################################################################

###########################################################################
# CUSTOM PARAMETERS
###########################################################################

THREADS=4
CUT_MIN_LENGTH=180
CUT_ERROR_RATE=0.2
CUT_OVERLAP=15
CUT_QUALITY_CUTOFF=22
CUT_ADAPTER_g=ACTCCTACGGGAGGCAGCAG
CUT_ADAPTER_G=GGACTACHVGGGTWTCTAAT
DADA_TRUNC_LEN_F=0
DADA_TRUNC_LEN_R=0
DADA_TRUNC_LEFT_F=0
DADA_TRUNC_LEFT_R=0
DADA_POOLING=pseudo
DADA_CHIMERA=pooled
DADA_MIN_FOLD=1
SAMPLING_DEPTH=93819
BLAST_PERC_IDENTITY=0.8

###########################################################################
# CREATE DIRECTORIES
###########################################################################

mkdir clean_data
mkdir clean_data/1-fastqc_raw
mkdir clean_data/2-trimmed
mkdir clean_data/2-removed
mkdir clean_data/3-fastqc_trimmed
mkdir clean_data/3-fastqc_removed

DIR_IMPORT=qiime2_pip/1-import
mkdir $DIR_IMPORT
DIR_DADA=qiime2_pip/2-dada
mkdir $DIR_DADA
DIR_PHYLO=qiime2_pip/3-phylogeny
mkdir $DIR_PHYLO
DIR_DIV=qiime2_pip/4-diversity
mkdir $DIR_DIV
DIR_TAX=qiime2_pip/5-taxonomy
mkdir $DIR_TAX

DIR_DB=databases/qiime2-silva138_ssuref_nr99_full

###########################################################################
# CLEAN DATA
###########################################################################

conda activate bioinf_util_env

# Run fastqc on raw data

fastqc raw_data/*fq.gz \
  --threads ${THREADS} \
  --outdir clean_data/1-fastqc_raw

multiqc clean_data/1-fastqc_raw/ --outdir clean_data/1-fastqc_raw

# Trim reads

function cut_fct () { \
cutadapt --cores=${THREADS} \
        --minimum-length ${CUT_MIN_LENGTH} \
        --error-rate ${CUT_ERROR_RATE} \
        --overlap ${CUT_OVERLAP} \
        --quality-cutoff ${CUT_QUALITY_CUTOFF} \
        -g ${CUT_ADAPTER_g} \
        -G ${CUT_ADAPTER_G} \
        -o clean_data/2-trimmed/$1_1.fq \
        -p clean_data/2-trimmed/$1_2.fq \
        --untrimmed-output clean_data/2-removed/$1_1.fq \
        --untrimmed-paired-output clean_data/2-removed/$1_2.fq \
        raw_data/$1_1.fq.gz \
        raw_data/$1_2.fq.gz; \
}

cut_fct 20.rawdata > clean_data/2-trimmed/cutadapt_20.log
cut_fct 60.rawdata > clean_data/2-trimmed/cutadapt_60.log
cut_fct 100.rawdata > clean_data/2-trimmed/cutadapt_100.log

# Run fastqc on trimmed and removed data

fastqc clean_data/2-trimmed/*.fq \
  --threads ${THREADS} \
  --outdir clean_data/3-fastqc_trimmed

multiqc clean_data/3-fastqc_trimmed/ --outdir clean_data/3-fastqc_trimmed

fastqc clean_data/2-removed/*.fq \
  --threads ${THREADS} \
  --outdir clean_data/3-fastqc_removed

multiqc clean_data/3-fastqc_removed/ --outdir clean_data/3-fastqc_removed

# Save to gz files

gzip clean_data/2-trimmed/*.fq
gzip clean_data/2-removed/*.fq

###########################################################################
# PREPARE DATA FOR QIIME2 WORKFLOW
###########################################################################

conda activate qiime2_env

# generate qiime file

qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path qiime_fastq_manifest.tsv \
    --output-path ${DIR_IMPORT}/kapil.qza \
    --input-format PairedEndFastqManifestPhred33V2

# Summarize demultiplexed and trimmed reads

qiime demux summarize \
  --i-data ${DIR_IMPORT}/kapil.qza \
  --o-visualization ${DIR_IMPORT}/kapil.qzv

qiime tools extract \
  --input-path ${DIR_IMPORT}/kapil.qzv \
  --output-path ${DIR_IMPORT}/summary

###########################################################################
# RUN DADA2
###########################################################################

# DADA2 denoising

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs ${DIR_IMPORT}/kapil.qza \
  --p-n-threads ${THREADS} \
  --p-trunc-len-f ${DADA_TRUNC_LEN_F} \
  --p-trunc-len-r ${DADA_TRUNC_LEN_R} \
  --p-trim-left-f ${DADA_TRUNC_LEFT_F} \
  --p-trim-left-r ${DADA_TRUNC_LEFT_R} \
  --p-pooling-method ${DADA_POOLING} \
  --p-chimera-method ${DADA_CHIMERA} \
  --p-min-fold-parent-over-abundance ${DADA_MIN_FOLD} \
  --o-table ${DIR_DADA}/table.qza \
  --o-representative-sequences ${DIR_DADA}/rep-seqs.qza \
  --o-denoising-stats ${DIR_DADA}/stats.qza

# extract ASV frequencies over all samples

qiime feature-table summarize \
  --i-table ${DIR_DADA}/table.qza \
  --o-visualization ${DIR_DADA}/table.qzv \
  --m-sample-metadata-file qiime_metadata.tsv

qiime tools extract \
  --input-path ${DIR_DADA}/table.qzv \
  --output-path ${DIR_DADA}/table

# extract ASV sequences

qiime feature-table tabulate-seqs \
  --i-data ${DIR_DADA}/rep-seqs.qza \
  --o-visualization ${DIR_DADA}/rep-seqs.qzv

qiime tools extract \
  --input-path ${DIR_DADA}/rep-seqs.qzv \
  --output-path ${DIR_DADA}/rep-seqs

# stats chimeras and used reads

qiime metadata tabulate \
  --m-input-file ${DIR_DADA}/stats.qza \
  --o-visualization ${DIR_DADA}/stats.qzv

qiime tools extract \
  --input-path ${DIR_DADA}/stats.qzv \
  --output-path ${DIR_DADA}/stats

# extract ASVs for each sample

qiime tools export \
  --input-path ${DIR_DADA}/table.qza \
  --output-path ${DIR_DADA}/table-biom

biom convert \
  -i ${DIR_DADA}/table-biom/feature-table.biom \
  -o ${DIR_DADA}/table-biom/asv-table.tsv \
  --to-tsv
  
###########################################################################
# MAKE PHYLOGENETIC ANALYSIS
###########################################################################

# phylogenetic tree

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences ${DIR_DADA}/rep-seqs.qza \
  --o-alignment ${DIR_PHYLO}/aligned-rep-seqs.qza \
  --o-masked-alignment ${DIR_PHYLO}/masked-aligned-rep-seqs.qza \
  --o-tree ${DIR_PHYLO}/unrooted-tree.qza \
  --o-rooted-tree ${DIR_PHYLO}/rooted-tree.qza

qiime tools extract \
  --input-path ${DIR_PHYLO}/aligned-rep-seqs.qza \
  --output-path ${DIR_PHYLO}/aligned-rep-seqs

qiime tools extract \
  --input-path ${DIR_PHYLO}/masked-aligned-rep-seqs.qza \
  --output-path ${DIR_PHYLO}/masked-aligned-rep-seqs

qiime tools extract \
  --input-path ${DIR_PHYLO}/unrooted-tree.qza \
  --output-path ${DIR_PHYLO}/unrooted-tree

qiime tools extract \
  --input-path ${DIR_PHYLO}/rooted-tree.qza \
  --output-path ${DIR_PHYLO}/rooted-tree

###########################################################################
# MAKE DIVERSITY ANALYSIS
###########################################################################

# diversity analysis with rarefaction

qiime diversity core-metrics-phylogenetic \
  --i-table ${DIR_DADA}/table.qza \
  --i-phylogeny ${DIR_PHYLO}/rooted-tree.qza \
  --p-sampling-depth ${SAMPLING_DEPTH} \
  --m-metadata-file qiime_metadata.tsv \
  --output-dir ${DIR_DIV}/core-metrics-results

# Richness

qiime tools extract \
  --input-path ${DIR_DIV}/core-metrics-results/observed_features_vector.qza \
  --output-path ${DIR_DIV}/observed_features

# Shannon

qiime diversity alpha-group-significance \
  --i-alpha-diversity ${DIR_DIV}/core-metrics-results/shannon_vector.qza \
  --m-metadata-file qiime_metadata2.tsv \
  --o-visualization ${DIR_DIV}/shannon_vector.qzv

qiime tools extract \
  --input-path ${DIR_DIV}/core-metrics-results/shannon_vector.qza \
  --output-path ${DIR_DIV}/shannon_vector_notest

qiime tools extract \
  --input-path ${DIR_DIV}/shannon_vector.qzv \
  --output-path ${DIR_DIV}/shannon_vector

# Clustering

qiime emperor plot \
  --i-pcoa ${DIR_DIV}/core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file qiime_metadata.tsv \
  --o-visualization ${DIR_DIV}/unweighted-unifrac-emperor.qzv

qiime tools extract \
  --input-path ${DIR_DIV}/core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --output-path ${DIR_DIV}/unweighted_unifrac_pcoa_results

qiime tools extract \
  --input-path ${DIR_DIV}/core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --output-path ${DIR_DIV}/unweighted_unifrac_distance_matrix

qiime tools extract \
  --input-path ${DIR_DIV}/unweighted-unifrac-emperor.qzv \
  --output-path ${DIR_DIV}/unweighted-unifrac-emperor

###########################################################################
# MAKE TAXONOMIC ANALYSIS
###########################################################################

# Apply Consensus BLAST

qiime feature-classifier classify-consensus-blast \
  --p-num-threads ${THREADS} \
  --i-query ${DIR_DADA}/rep-seqs.qza \
  --i-reference-reads ${DIR_DB}/silva-138-99-seqs.qza \
  --i-reference-taxonomy ${DIR_DB}/silva-138-99-tax.qza \
  --p-perc-identity ${BLAST_PERC_IDENTITY} \
  --o-classification ${DIR_TAX}/blast-taxa.qza \
  --o-search-results ${DIR_TAX}/blast-tophits.qza \
  --verbose

# Extract taxa of ASVs

qiime metadata tabulate \
  --m-input-file ${DIR_TAX}/blast-taxa.qza \
  --o-visualization ${DIR_TAX}/blast-taxa.qzv

qiime tools extract \
  --input-path ${DIR_TAX}/blast-taxa.qzv \
  --output-path ${DIR_TAX}/blast-taxa-tab

# Aggregate ASV frequencies per taxon (barplots)

qiime taxa barplot \
  --i-table ${DIR_DADA}/table.qza \
  --i-taxonomy ${DIR_TAX}/blast-taxa.qza \
  --m-metadata-file qiime_metadata.tsv \
  --o-visualization ${DIR_TAX}/blast-taxa-barplots.qzv

qiime tools extract \
  --input-path ${DIR_TAX}/blast-taxa-barplots.qzv \
  --output-path ${DIR_TAX}/blast-taxa-barplots
