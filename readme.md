# Taxonomic profiling of 16S rDNA metabarcoding data

* _author: Michaël Pierrelée_
* _date: April 2024_
* _Technical University of Denmark, DTU Sustain - Lyngby, DENMARK_

The pipeline is based on QIIME2's workflow using DADA2 (https://docs.qiime2.org/2024.2/)

Main versions:
- preprocessing: cutadapt (4.8), fastqc (0.12.1) and multiqc (1.21)
- processing: QIIME2 (2024.2)
- reporting: python (3.12.2), scipy (1.12.0), pandas (2.1.4), numpy (1.26.4), matplotlib (3.8.0), seaborn (0.12.2)

Folders:
- **raw_data/**: Raw fastq files downloaded from SRA (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1111705): `20.rawdata\_1.fq.gz`, `20.rawdata\_2.fq.gz`, `60.rawdata\_1.fq.gz`, `60.rawdata\_1.fq.gz`, `100.rawdata\_1.fq.gz`, `100.rawdata\_2.fq.gz`
- **clean_data/**: Trimmed fastq files by cutadapt with the statistics from fastqc and multiqc
- **databases/**: SILVA database in the subfolder `qiime2-silva138_ssuref_nr99_full/`
- **qiime2_pip/**: intermediary files from QIIME2
- **scripts/**: bash script and python notebook
- **results/**: figures and summary table
- **conda_environments/**: yaml specifications files of the actual conda environments (the qiime2 env was downloaded from https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.2-py38-linux-conda.yml)

Replication:
1. Clone the repository and cd there, as your working directory
2. Install the environments (and your own jupyter environment to run the notebook) from `conda_environments/`
3. Create the folders `raw_data/`, `databases/qiime2-silva138_ssuref_nr99_full/` and `qiime2_pip/`
4. Download the raw data and the SILVA database
5. Update the file `qiime_fastq_manifest.tsv` with the paths to the expected location of the trimmed fastq files (eg `/pathToTheWorkingDirectory/clean_data/2-trimmed/20.rawdata_2.fq.gz`)
6. Run the bash file `scripts/processing.sh`
7. Run the python notebook `scripts/reporting.ipynb`
