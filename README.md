# CMV Genome Assembly & Betaherpesvirinae BLAST Pipeline  
Author: Manel  Sadok  

---

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

---

## Data Source (NCBI)

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

---

## Software & Tools Used

- Snakemake  
- Bowtie2  
- SPAdes  
- NCBI BLAST+  
- awk  
- Bash  

---

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

---

## How to Run

### Test Mode

```
snakemake --cores 4 PipelineReport.txt
```

### Full Mode

```
snakemake --cores 4 Sadok_PipelineReport.txt
```

---

## Output

Final report file:

```
Sadok_PipelineReport.txt
```

This report includes:

- Read counts before and after host filtering  
- Number of contigs >1000 bp  
- Total base pairs in contigs >1000 bp  
- Top 5 BLAST hits against Betaherpesvirinae  

---

## Reproducibility

This workflow is fully automated using Snakemake.  
All results can be regenerated from raw sequencing data using the provided:

- Snakefile  
- config.yaml  

Parallel execution is supported using multiple cores.
