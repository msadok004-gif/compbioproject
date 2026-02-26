# COMP 383 Python Pipeline Project

## Overview
This repo contains a Snakemake pipeline that automates analysis steps 2–5 and produces `PipelineReport.txt`.


- snakemake
- python3
- bowtie2
- spades.py (SPAdes)
- blast+ (makeblastdb, blastn)
- NCBI datasets CLI OR Biopython (for retrieving HCMV reference)

## Quick start (test data should take around two minutes)
1) Clone repo
2) Install dependencies
3) Run:
   snakemake -j 1 --use-conda

##  Run
Place full FASTQ files in data/raw/ (see Step 1 instructions below), then run the same snakemake command.

## Step 1 (documented, not automated)
Download these SRA runs and convert to paired FASTQ:
- SRR5660030
- SRR5660033
- SRR5660044
- SRR5660045


