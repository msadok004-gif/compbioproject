# CMV Genome Assembly & Betaherpesvirinae BLAST Pipeline  
Author: Manel  Sadok  


## Project Overview

This project implements a fully reproducible bioinformatics workflow using **Snakemake** to analyze Human Cytomegalovirus (HCMV) sequencing data.

The pipeline performs:

1. Host read filtering using Bowtie2  
2. De novo genome assembly using SPAdes  
3. Contig statistics calculation (>1000 bp)  
4. Longest contig extraction  
5. BLAST analysis against the Betaherpesvirinae database  
6. Automated report generation  

All steps are automated and reproducible.



## Data Source: NCBI

All reference genomes, sequencing data, and BLAST databases were obtained from:

**National Center for Biotechnology Information (NCBI)**  
https://www.ncbi.nlm.nih.gov/

### Reference Genome
Human Cytomegalovirus reference genome:  
GCF_000845245.1  

### Sequence Read Archive (SRA) Accessions
The following SRA accessions were analyzed:

- SRR5660030  
- SRR5660033  
- SRR5660044  
- SRR5660045  

### BLAST Database
Betaherpesvirinae nucleotide sequences downloaded from NCBI.


## Installation

It is recommended to use **conda** to install dependencies.

### Option 1:  Using Conda

```
conda create -n cmv_pipeline snakemake bowtie2 spades blast -c bioconda -c conda-forge
conda activate cmv_pipeline
```

### Option 2:  Manual Installation

Ensure the following tools are installed and available in your PATH:

- snakemake  
- bowtie2  
- spades.py  
- makeblastdb  
- blastn  

You can verify installation with:

```
snakemake --version
bowtie2 --version
spades.py --version
blastn -version
```



## Directory Structure

```
data/
    ref/
    test/
    full/

results/
    bowtie2_filtered/
    assembly/
    stats/
    blast/
```


## How to Run

After cloning the repository:

```
git clone https://github.com/msadok004-gif/compbioproject.git
cd compbioproject
```

### Run Test Mode

```
snakemake --cores 4 PipelineReport.txt
```

This will run the pipeline using the included sample test data.

### Run Full Mode

```
snakemake --cores 4 Sadok_PipelineReport.txt
```


## Output

The final report file:

```
Sadok_PipelineReport.txt
```

This report includes:

- Read counts before and after host filtering  
- Number of contigs >1000 bp  
- Total base pairs in contigs >1000 bp  
- Top 5 BLAST hits against Betaherpesvirinae  


## Reproducibility

This workflow is fully automated using Snakemake.  

All results can be regenerated from raw sequencing data using:

- Snakefile  
- config.yaml  
- Included sample test data  

Parallel execution is supported using multiple cores.

No absolute file paths are hardcoded; all paths are relative to the project directory.
