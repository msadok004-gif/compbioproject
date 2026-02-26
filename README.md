# COMP 383 Python Pipeline Project

Author: Manel Sadok



## Project Overview

This project implements a fully reproducible bioinformatics workflow using Snakemake to analyze Human Cytomegalovirus (HCMV) sequencing data.

For each sample, the pipeline performs:

1. Download of the HCMV reference genome from NCBI
2. Construction of a Bowtie2 index
3. Alignment of paired-end reads to the HCMV reference
4. Filtering to retain only mapped read pairs
5. Genome assembly using SPAdes
6. Contig statistics calculation for contigs >1000 bp
7. Extraction of the longest contig
8. BLAST of the longest contig against a Betaherpesvirinae database
9. Automated report generation

All steps are automated and reproducible through Snakemake.



## Data Source (NCBI)

All reference genomes and BLAST database sequences are obtained from:

- NCBI Genome database
- NCBI SRA
- Betaherpesvirinae genomes downloaded using the NCBI datasets CLI

Reference genome accession used:

GCF_000845245.1 (HCMV)



## Repository Structure

compbioproject/
│
├── Snakefile
├── config.yaml
├── README.md
├── data/
│   └── test/                 Small FASTQ test files (<50MB each)
├── results/                  Generated automatically during pipeline runs
└── Sadok_PipelineReport.txt  Full-run report



## Dependencies

The following software must be installed and available in your PATH:

- snakemake
- bowtie2
- spades.py
- makeblastdb
- blastn
- NCBI datasets CLI
- unzip
- awk
- bash



## Installation (Recommended: Conda)

Create and activate a conda environment with required tools:

conda create -n cmv_pipeline snakemake bowtie2 spades blast ncbi-datasets-cli unzip -c bioconda -c conda-forge
conda activate cmv_pipeline

Verify installation:

snakemake --version
bowtie2 --version
spades.py --version
blastn -version
datasets version



## How to Run

Clone the repository:

git clone https://github.com/msadok004-gif/compbioproject.git
cd compbioproject



Run Test Mode:

snakemake --cores 4 PipelineReport.txt

This runs the pipeline using the included small FASTQ test files.
Each test FASTQ file is under 50MB.
The test run completes in under 2 minutes.

Output:
PipelineReport.txt



Run Full Mode (All Input Reads):

snakemake --cores 4 Sadok_PipelineReport.txt

Output:
Sadok_PipelineReport.txt



## Output Description

Each report includes:

- Read pairs before and after filtering
- Assembly output directory
- Number of contigs >1000 bp
- Total base pairs in contigs >1000 bp
- Top 5 BLAST hits in tabular format

