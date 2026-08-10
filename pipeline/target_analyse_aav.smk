# ============================================================
# AAV integration detection pipeline
#
# Realigns reads to AAV vector references and quantifies AAV
# integrations at the CRISPR target sites:
#
#   1. prepare AAV reference fastas (one per sgRNA) + BWA index
#   2. BWA-MEM alignment of the target-sequencing reads to each
#      AAV reference (input: fastp/BWA output of target_analyse_2.smk)
#   3. coverage of the AAV vector (samtools depth)
#   4. genome positions of AAV-mapped reads
#      (filter the genome BAM by AAV-mapped read names)
#   5. genome coverage of AAV-mapped reads
#   6. insertion-ratio tables (count_aav_ratio, count_aav_ratio2)
#   7. QC plots (plot_aav_insertion)
#
# Config keys:
#   samples_table : TSV with columns read1, read2, sample_name
#   primers_table : CSV with columns name, chr, pos, length, strand
#   conda_env     : conda environment used for script rules
#   min_mapq      : minimum mapping quality for AAV-mapped reads
#
# The AAV reference sequences are expected under
# data/raw/aav/aav_{sgRNA}.fa (included in this repository).
# ============================================================

import sys
from pathlib import Path

import pandas as pd

base_dir = workflow.workdir_init
if base_dir not in sys.path:
    sys.path.insert(0, base_dir)

from src.config import Location

# ============================================================
# Parameters
# ============================================================

DEFAULT_PARAMS = {
    "samples_table": "configs/target_sequence/samples.tsv",
    "primers_table": "configs/target_sequence/target_primers.csv",
    "min_mapq"     : 10,
}
params = {**DEFAULT_PARAMS, **config}

COND_ENV = config.get("conda_env", "long-term-crispr-in-liver")

sgRNA_list = ["angptl3_g4", "con", "pcsk9_g1", "pcsk9_g3"]
sample_df  = pd.read_table(params["samples_table"], sep="\t")  # three columns: read1, read2, sample_name
primers_df = pd.read_csv(params["primers_table"], dtype={"chr": str})

# ============================================================
# Preprocessing inputs (outputs of target_analyse_2.smk)
# ============================================================

last_processed_dir = Location.processed_data / "targetsequence"
processed_dir      = Location.processed_data / "targetsequence_aav"

processed_read1 = str(last_processed_dir / "{sample_name}" / "fastp" / "trimmed_R1.fastq.gz")
processed_read2 = str(last_processed_dir / "{sample_name}" / "fastp" / "trimmed_R2.fastq.gz")

# ============================================================
# Rules
# ============================================================

rule all:
    input:
        aav_reference = [
            processed_dir / "aav_reference" / f"aav_{sgRNA}" / f"aav_{sgRNA}.fa.amb"
            for sgRNA in sgRNA_list
        ],
        mapped_res = [
            processed_dir / "mapping" / sample_name / f"{sgRNA}.coverage.txt"
            for sample_name in sample_df.sample_name
            for sgRNA in ["con"]
        ],
        genome_filtered_coverage = [
            processed_dir / "mapping" / sample_name / f"{sgRNA}.genome.filtered.coverage.txt"
            for sample_name in sample_df.sample_name
            for sgRNA in ["con", "angptl3_g4"]
        ],
        genome_coverage = [
            processed_dir / "mapping" / sample_name / f"{sgRNA}.genome.coverage.txt"
            for sample_name in sample_df.sample_name
            for sgRNA in ["con", "angptl3_g4"]
        ],
        aav_insertion_plots = [
            Path("reports/figures/aav") / sample_name / f"{sgRNA}.aav_insertion.svg"
            for sample_name in sample_df.sample_name
            for sgRNA in ["con", "angptl3_g4"]
        ],
        aav_insertion_ratios = [
            Path(processed_dir / "mapping" / sample_name / f"{sgRNA}.aav_insertion_ratio.txt")
            for sample_name in sample_df.sample_name
            for sgRNA in ["con", "angptl3_g4", "pcsk9_g1", "pcsk9_g3"]
        ],
        aav_insertion_ratios2 = [
            Path(processed_dir / "mapping" / sample_name / f"{sgRNA}.aav_insertion_ratio.2.txt")
            for sample_name in sample_df.sample_name
            for sgRNA in ["pcsk9_g1", "pcsk9_g3", "angptl3_g4"]
        ]
# ---

rule prepare_aav_fasta:
    input:
        fasta = "data/raw/aav/aav_{sgRNA}.fa"
    output:
        fasta = processed_dir / "aav_reference" / "aav_{sgRNA}" / "aav_{sgRNA}.fa",
        fai   = processed_dir / "aav_reference" / "aav_{sgRNA}" / "aav_{sgRNA}.fa.fai"
    shell:
        '''
        cp {input.fasta} {output.fasta}
        samtools faidx {output.fasta}
        '''
# ---

rule prepare_aav_bwa_index:
    input:
        fasta = processed_dir / "aav_reference" / "aav_{sgRNA}" / "aav_{sgRNA}.fa"
    output:
        index_files = expand(processed_dir / "aav_reference" / "aav_{{sgRNA}}" / "aav_{{sgRNA}}.fa.{ext}", ext=["amb", "ann", "bwt", "pac", "sa"])
    shell:
        '''
        bwa index {input.fasta}
        '''
# ---

rule bwa_mapping:
    input:
        refgenome = processed_dir / "aav_reference" / "aav_{sgRNA}" / "aav_{sgRNA}.fa",
        R1        = processed_read1,
        R2        = processed_read2
    output:
        bam       = str(processed_dir / "mapping" / "{sample_name}" / "{sgRNA}.mapped.bam"),
        bam_index = str(processed_dir / "mapping" / "{sample_name}" / "{sgRNA}.mapped.bam.bai")
    threads: 8
    shell:
        '''
        bwa mem -Y -a -t {threads} {input.refgenome} {input.R1} {input.R2} | samtools sort -@ {threads} -o {output.bam} && samtools index -@ {threads} {output.bam}
        '''
# ---

rule samtools_coverage:
    input:
        bam = str(processed_dir / "mapping" / "{sample_name}" / "{sgRNA}.mapped.bam")
    output:
        coverage = str(processed_dir / "mapping" / "{sample_name}" / "{sgRNA}.coverage.txt")
    params:
        min_mapq = params["min_mapq"]
    shell:
        '''
        samtools depth -q {params.min_mapq} {input.bam} > {output.coverage}
        '''
# ---

# Filter the genome BAM to reads that also map to the AAV vector
rule filter_mapped_bam:
    input:
        aav_bam     = processed_dir / "mapping" / "{sample_name}" / "{sgRNA}.mapped.bam",
        genome_bam  = last_processed_dir / "{sample_name}" / "mapped.bam"
    output:
        filtered_reads = processed_dir / "mapping" / "{sample_name}" / "{sgRNA}.mapped.readnames.txt",
        filtered_bam   = processed_dir / "mapping" / "{sample_name}" / "{sgRNA}.genome.filtered.bam"
    params:
        min_mapq = params["min_mapq"]
    shell:
        '''
        # extract read names from the AAV-mapped BAM
        samtools view -q {params.min_mapq} {input.aav_bam} | cut -f1 > {output.filtered_reads}
        # filter the genome BAM by those read names
        samtools view -b -N {output.filtered_reads} {input.genome_bam} > {output.filtered_bam}
        samtools index {output.filtered_bam}
        '''
# ---

rule samtools_aav_insertion_coverage:
    input:
        bam = processed_dir / "mapping" / "{sample_name}" / "{sgRNA}.genome.filtered.bam"
    output:
        coverage = processed_dir / "mapping" / "{sample_name}" / "{sgRNA}.genome.filtered.coverage.txt"
    shell:
        '''
        samtools depth {input.bam} > {output.coverage}
        '''
# ---

rule samtools_genome_coverage:
    input:
        bam = last_processed_dir / "{sample_name}" / "mapped.bam"
    output:
        coverage = processed_dir / "mapping" / "{sample_name}" / "{sgRNA}.genome.coverage.txt"
    params:
        min_mapq = params["min_mapq"]
    shell:
        '''
        # filter mapping quality >= min_mapq
        samtools depth -q {params.min_mapq} {input.bam} > {output.coverage}
        '''
# ---

rule plot_aav_insertion:
    input:
        full_genome_cov = processed_dir / "mapping" / "{sample_name}" / "{sgRNA}.genome.coverage.txt",
        genome_cov      = processed_dir / "mapping" / "{sample_name}" / "{sgRNA}.genome.filtered.coverage.txt",
        aav_cov         = processed_dir / "mapping" / "{sample_name}" / "{sgRNA}.coverage.txt"
    output:
        full_genome_plot = Path("reports/figures/aav") / "{sample_name}" / "{sgRNA}.full_genome_coverage.svg",
        genome_plot      = Path("reports/figures/aav") / "{sample_name}" / "{sgRNA}.aav_insertion.svg",
        aav_plot         = Path("reports/figures/aav") / "{sample_name}" / "{sgRNA}.aav_coverage.svg"
    conda: COND_ENV
    params:
        sample_name = "{sample_name}",
        sgRNA       = "{sgRNA}"
    script:
        "scripts/aav_plot_insertion.R"
# ---

rule count_aav_ratio:
    input:
        genome_bam   = last_processed_dir / "{sample_name}" / "mapped.bam",
        filtered_bam = processed_dir / "mapping" / "{sample_name}" / "{sgRNA}.genome.filtered.bam"
    output:
        ratio_file = processed_dir / "mapping" / "{sample_name}" / "{sgRNA}.aav_insertion_ratio.txt"
    conda: COND_ENV
    shell:
        '''
        genome_reads_count=$(samtools view -c -F 4 {input.genome_bam})
        filtered_reads_count=$(samtools view -c -F 4 {input.filtered_bam})
        ratio=$(echo "scale=8; $filtered_reads_count / $genome_reads_count" | bc)
        echo -e "genome_reads_count\tfiltered_reads_count\taav_insertion_ratio" > {output.ratio_file}
        echo -e "$genome_reads_count\t$filtered_reads_count\t$ratio" >> {output.ratio_file}
        '''


def get_f_edit_reads_file(wildcards):
    sgRNA = wildcards.sgRNA.replace("_", "-")
    return last_processed_dir / f"{wildcards.sample_name}" / f"{sgRNA}-f" / "intarget_qnames__cover.txt"


def get_r_edit_reads_file(wildcards):
    sgRNA = wildcards.sgRNA.replace("_", "-")
    return last_processed_dir / f"{wildcards.sample_name}" / f"{sgRNA}-r" / "intarget_qnames__cover.txt"


rule count_aav_ratio2:
    input:
        aav_reads    = processed_dir / "mapping" / "{sample_name}" / "con.mapped.readnames.txt",
        f_edit_reads = get_f_edit_reads_file,
        r_edit_reads = get_r_edit_reads_file
    output:
        ratio_file = processed_dir / "mapping" / "{sample_name}" / "{sgRNA}.aav_insertion_ratio.2.txt"
    conda: COND_ENV
    script:
        "scripts/aav_count_insertion_ratio.py"
