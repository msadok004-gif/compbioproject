# COMP 383 Python Pipeline Project  
**Author:** Manel Sadok  


## Project Overview

This project implements a fully reproducible bioinformatics workflow using **Snakemake** to analyze Human Cytomegalovirus (HCMV) sequencing data.

For each sample, the pipeline performs:

- Download of the HCMV reference genome from NCBI  
- Construction of a Bowtie2 index  
- Alignment of paired-end reads to the HCMV reference  
- Filtering to retain only mapped read pairs  
- Genome assembly using SPAdes  
- Contig statistics calculation for contigs >1000 bp  
- Extraction of the longest contig  
- BLAST of the longest contig against a Betaherpesvirinae database  
- Automated report generation  

All steps are automated and reproducible through Snakemake.


## Data Source (NCBI)

All reference genomes and BLAST database sequences are obtained from:

- NCBI Genome database  
- NCBI SRA  
- Betaherpesvirinae genomes downloaded using the NCBI datasets CLI  

Reference genome accession used:

```
GCF_000845245.1
```


## Dependencies

The following software must be installed and available in your `PATH`:

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

```bash
conda create -n cmv_pipeline snakemake bowtie2 spades blast ncbi-datasets-cli unzip -c bioconda -c conda-forge
conda activate cmv_pipeline
```

Verify installation:

```bash
snakemake --version
bowtie2 --version
spades.py --version
blastn -version
datasets version
```


## How to Run

###1. Clone the repository

```bash
git clone https://github.com/msadok004-gif/compbioproject.git
cd compbioproject
```


### 2. Run Test Mode

```bash
snakemake --cores 4 PipelineReport.txt
```

This runs the pipeline using the included small FASTQ test files.

- Each test FASTQ file is under 50MB  
- The test run completes in under 2 minutes  

**Output:**

```
PipelineReport.txt
```


### 3. Run Full Mode (All Input Reads)

```bash
snakemake --cores 4 Sadok_PipelineReport.txt
```

**Output:**

```
Sadok_PipelineReport.txt
```


## Output Description

Each report includes:

- Read pairs before and after filtering  
- Assembly output directory  
- Number of contigs >1000 bp  
- Total base pairs in contigs >1000 bp  
- Top 5 BLAST hits in tabular format  


## Reproducibility

After cloning this repository and installing dependencies, the entire workflow (Steps 2–5 of the project instructions) can be executed with a **single Snakemake command** as shown above.

No file movement or path editing is required.
