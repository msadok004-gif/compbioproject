import os

configfile: "config.yaml"

SAMPLES = config["samples"]
THREADS = int(config["params"]["threads"])
K = int(config["params"]["spades_k"])

RAW = config["paths"]["raw_dir"]
TEST = config["paths"]["test_dir"]
REFDIR = config["paths"]["ref_dir"]
RESULTS = config["paths"]["results_dir"]

HCMV_FNA = config["reference"]["hcmv_fna"]
HCMV_INDEX = config["reference"]["bowtie2_index_prefix"]

BETA_FNA = config["reference"]["betaherpes_fna"]
BLAST_DB = config["reference"]["blast_db_prefix"]

LAST = config.get("last_name", "LASTNAME")

BLAST_COLS = "sacc pident length qstart qend sstart send bitscore evalue stitle"


def in_fastq(sample, mode, mate):
    if mode == "test":
        return f"{TEST}/{sample}_test_{mate}.fastq"
    return f"{RAW}/{sample}_{mate}.fastq"


rule all:
    input:
        "PipelineReport.txt",
        f"{LAST}_PipelineReport.txt"


########################################
# Step 2 prereq: HCMV reference + Bowtie2 index
########################################

rule fetch_hcmv_ref:
    output:
        HCMV_FNA
    params:
        acc="GCF_000845245.1",
        refdir=REFDIR
    shell:
        r"""
        mkdir -p {params.refdir}
        rm -rf {params.refdir}/hcmv_dataset
        rm -f {params.refdir}/hcmv.zip

        datasets download genome accession {params.acc} --include genome -o {params.refdir}/hcmv.zip
        unzip -o {params.refdir}/hcmv.zip -d {params.refdir}/hcmv_dataset >/dev/null

        fna=$(find {params.refdir}/hcmv_dataset -name "*genomic.fna" | head -n 1)
        if [ -z "$fna" ]; then
          echo "ERROR: Could not find genomic.fna for HCMV." >&2
          exit 1
        fi
        cp "$fna" {output}
        """

rule build_hcmv_bowtie2_index:
    input:
        HCMV_FNA
    output:
        HCMV_INDEX + ".1.bt2"
    params:
        index_prefix=HCMV_INDEX
    shell:
        r"""
        mkdir -p $(dirname {params.index_prefix})
        bowtie2-build {input} {params.index_prefix}
        """


########################################
# Step 5 prereq: Betaherpesvirinae FASTA + BLAST DB
########################################

rule fetch_betaherpes_fna:
    output:
        BETA_FNA
    params:
        refdir=REFDIR
    shell:
        r"""
        mkdir -p $(dirname {output})
        rm -rf {params.refdir}/beta_dataset
        rm -f {params.refdir}/beta.zip

        datasets download genome taxon Betaherpesvirinae --include genome -o {params.refdir}/beta.zip
        unzip -o {params.refdir}/beta.zip -d {params.refdir}/beta_dataset >/dev/null

        find {params.refdir}/beta_dataset -name "*genomic.fna" -print0 | xargs -0 cat > {output}

        if [ ! -s {output} ]; then
          echo "ERROR: Betaherpesvirinae FASTA is empty." >&2
          exit 1
        fi
        """

rule make_betaherpes_blastdb:
    input:
        BETA_FNA
    output:
        BLAST_DB + ".nin"
    params:
        blast_prefix=BLAST_DB
    shell:
        r"""
        mkdir -p $(dirname {params.blast_prefix})
        makeblastdb -in {input} -dbtype nucl -out {params.blast_prefix}
        """


########################################
# Step 2: Bowtie2 filter (keep mapped pairs only)
########################################

rule bowtie2_filter:
    input:
        idx=HCMV_INDEX + ".1.bt2",
        r1=lambda wc: in_fastq(wc.sample, wc.mode, 1),
        r2=lambda wc: in_fastq(wc.sample, wc.mode, 2)
    output:
        r1=f"{RESULTS}/bowtie2_filtered/{{mode}}/{{sample}}.filtered.1.fastq",
        r2=f"{RESULTS}/bowtie2_filtered/{{mode}}/{{sample}}.filtered.2.fastq"
    threads: THREADS
    params:
        results=RESULTS,
        index_prefix=HCMV_INDEX
    shell:
        r"""
        mkdir -p {params.results}/bowtie2_filtered/{wildcards.mode}

        bowtie2 -x {params.index_prefix} \
          -1 {input.r1} -2 {input.r2} \
          --threads {threads} --very-sensitive \
          --al-conc {params.results}/bowtie2_filtered/{wildcards.mode}/{wildcards.sample}.tmp_aligned.fastq \
          -S /dev/null

        mv {params.results}/bowtie2_filtered/{wildcards.mode}/{wildcards.sample}.tmp_aligned.1.fastq {output.r1}
        mv {params.results}/bowtie2_filtered/{wildcards.mode}/{wildcards.sample}.tmp_aligned.2.fastq {output.r2}
        """


########################################
# Step 2: Read-count stats (before/after filtering)
########################################

rule read_counts:
    input:
        before_r1=lambda wc: in_fastq(wc.sample, wc.mode, 1),
        after_r1=f"{RESULTS}/bowtie2_filtered/{{mode}}/{{sample}}.filtered.1.fastq"
    output:
        reads_tsv=f"{RESULTS}/stats/{{mode}}/{{sample}}.reads.tsv"
    params:
        results=RESULTS
    shell:
        r"""
        mkdir -p {params.results}/stats/{wildcards.mode}

        before_lines=$(wc -l < {input.before_r1})
        after_lines=$(wc -l < {input.after_r1})

        before_pairs=$((before_lines / 4))
        after_pairs=$((after_lines / 4))

        echo -e "{wildcards.sample}\t$before_pairs\t$after_pairs" > {output.reads_tsv}
        """


########################################
# Step 3: SPAdes assembly (output contigs.fasta)
########################################

rule spades:
    input:
        r1=f"{RESULTS}/bowtie2_filtered/{{mode}}/{{sample}}.filtered.1.fastq",
        r2=f"{RESULTS}/bowtie2_filtered/{{mode}}/{{sample}}.filtered.2.fastq"
    output:
        contigs=f"{RESULTS}/assembly/{{mode}}/{{sample}}_spades/contigs.fasta"
    threads: THREADS
    params:
        results=RESULTS
    shell:
        r"""
        mkdir -p {params.results}/assembly/{wildcards.mode}
        outdir={params.results}/assembly/{wildcards.mode}/{wildcards.sample}_spades

        spades.py --only-assembler -k {K} -t {threads} \
          -1 {input.r1} -2 {input.r2} \
          -o "$outdir"

        test -s {output.contigs}
        """


########################################
# Step 4: Contig stats (>1000 bp)
########################################

rule contig_stats_gt1000:
    input:
        contigs=f"{RESULTS}/assembly/{{mode}}/{{sample}}_spades/contigs.fasta"
    output:
        stats=f"{RESULTS}/stats/{{mode}}/{{sample}}.contigs_gt1000.tsv"
    params:
        results=RESULTS
    shell:
        r"""
        mkdir -p {params.results}/stats/{wildcards.mode}

        count=$(awk '/^>/ {{if (seqlen > 1000) c++; seqlen=0; next}} {{seqlen += length($0)}} END {{if (seqlen > 1000) c++; print c+0}}' {input.contigs})
        total=$(awk '/^>/ {{if (seqlen > 1000) t+=seqlen; seqlen=0; next}} {{seqlen += length($0)}} END {{if (seqlen > 1000) t+=seqlen; print t+0}}' {input.contigs})

        echo -e "{wildcards.sample}\t$count\t$total" > {output.stats}
        """


########################################
# Step 5: Longest contig extraction
########################################

rule longest_contig:
    input:
        contigs=f"{RESULTS}/assembly/{{mode}}/{{sample}}_spades/contigs.fasta"
    output:
        longest=f"{RESULTS}/assembly/{{mode}}/{{sample}}_spades/longest_contig.fasta"
    shell:
        r"""
        awk '
        /^>/ {{
            if (seqlen > maxlen) {{maxlen=seqlen; maxseq=seq; maxheader=header}}
            header=$0; seq=""; seqlen=0; next
        }}
        {{
            seq=seq $0; seqlen+=length($0)
        }}
        END {{
            if (seqlen > maxlen) {{maxlen=seqlen; maxseq=seq; maxheader=header}}
            print maxheader
            print maxseq
        }}' {input.contigs} > {output.longest}
        """


########################################
# Step 5: BLAST longest contig (top 5 hits)
########################################

rule blast_longest:
    input:
        db=BLAST_DB + ".nin",
        longest=f"{RESULTS}/assembly/{{mode}}/{{sample}}_spades/longest_contig.fasta"
    output:
        out=f"{RESULTS}/blast/{{mode}}/{{sample}}.blast.tsv"
    params:
        results=RESULTS,
        blast_prefix=BLAST_DB
    shell:
        r"""
        mkdir -p {params.results}/blast/{wildcards.mode}

        blastn \
          -query {input.longest} \
          -db {params.blast_prefix} \
          -out {output.out} \
          -outfmt "6 {BLAST_COLS}" \
          -max_target_seqs 5 \
          -max_hsps 1
        """


########################################
# Reports (PYTHON - no brace conflicts)
########################################

def read_tsv(path):
    with open(path, "r") as f:
        line = f.readline().strip()
    return line.split("\t")

rule report_test:
    input:
        reads=expand(f"{RESULTS}/stats/test/{{s}}.reads.tsv", s=SAMPLES),
        contigs=expand(f"{RESULTS}/stats/test/{{s}}.contigs_gt1000.tsv", s=SAMPLES),
        blasts=expand(f"{RESULTS}/blast/test/{{s}}.blast.tsv", s=SAMPLES)
    output:
        "PipelineReport.txt"
    run:
        with open(output[0], "w") as out:
            out.write("Pipeline Report (TEST RUN)\n\n")
            for s in SAMPLES:
                out.write(f"Sample: {s}_test\n\n")

                out.write("Step 2: Bowtie2 filtering (HCMV GCF_000845245.1)\n")
                _, before, after = read_tsv(f"{RESULTS}/stats/test/{s}.reads.tsv")
                out.write(f"Read pairs BEFORE filtering: {before}\n")
                out.write(f"Read pairs AFTER  filtering: {after}\n\n")

                out.write(f"Step 3: SPAdes assembly (k={K})\n")
                out.write(f"Assembly directory: {RESULTS}/assembly/test/{s}_spades/\n\n")

                out.write("Step 4: Contig stats (contigs > 1000 bp)\n")
                _, c, t = read_tsv(f"{RESULTS}/stats/test/{s}.contigs_gt1000.tsv")
                out.write(f"Number of contigs >1000 bp: {c}\n")
                out.write(f"Total bp in contigs >1000 bp: {t}\n\n")

                out.write("Step 5: BLAST (top 5 hits; Betaherpesvirinae DB)\n")
                out.write(f"Columns: {BLAST_COLS}\n")
                blast_path = f"{RESULTS}/blast/test/{s}.blast.tsv"
                if os.path.exists(blast_path) and os.path.getsize(blast_path) > 0:
                    with open(blast_path, "r") as bf:
                        out.write(bf.read())
                else:
                    out.write("(no BLAST hits)\n")
                out.write("\n----------------------------------------\n\n")

rule report_full:
    input:
        reads=expand(f"{RESULTS}/stats/full/{{s}}.reads.tsv", s=SAMPLES),
        contigs=expand(f"{RESULTS}/stats/full/{{s}}.contigs_gt1000.tsv", s=SAMPLES),
        blasts=expand(f"{RESULTS}/blast/full/{{s}}.blast.tsv", s=SAMPLES)
    output:
        f"{LAST}_PipelineReport.txt"
    run:
        with open(output[0], "w") as out:
            out.write("Pipeline Report (FULL RUN)\n\n")
            for s in SAMPLES:
                out.write(f"Sample: {s}\n\n")

                out.write("Step 2: Bowtie2 filtering (HCMV GCF_000845245.1)\n")
                _, before, after = read_tsv(f"{RESULTS}/stats/full/{s}.reads.tsv")
                out.write(f"Read pairs BEFORE filtering: {before}\n")
                out.write(f"Read pairs AFTER  filtering: {after}\n\n")

                out.write(f"Step 3: SPAdes assembly (k={K})\n")
                out.write(f"Assembly directory: {RESULTS}/assembly/full/{s}_spades/\n\n")

                out.write("Step 4: Contig stats (contigs > 1000 bp)\n")
                _, c, t = read_tsv(f"{RESULTS}/stats/full/{s}.contigs_gt1000.tsv")
                out.write(f"Number of contigs >1000 bp: {c}\n")
                out.write(f"Total bp in contigs >1000 bp: {t}\n\n")

                out.write("Step 5: BLAST (top 5 hits; Betaherpesvirinae DB)\n")
                out.write(f"Columns: {BLAST_COLS}\n")
                blast_path = f"{RESULTS}/blast/full/{s}.blast.tsv"
                if os.path.exists(blast_path) and os.path.getsize(blast_path) > 0:
                    with open(blast_path, "r") as bf:
                        out.write(bf.read())
                else:
                    out.write("(no BLAST hits)\n")
                out.write("\n----------------------------------------\n\n")
