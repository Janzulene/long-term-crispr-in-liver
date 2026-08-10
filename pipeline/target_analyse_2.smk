# ============================================================
# Targeted enrichment sequencing pipeline: translocation and
# short-indel detection at CRISPR/Cas9 target sites.
#
# Steps:
#   1. fastp quality control / adapter trimming
#   2. BWA-MEM alignment to GRCm39
#   3. read counting, on-target read filtering (primer match,
#      strand, complexity)
#   4. translocation detection (check_translocation)
#   5. short-indel counting (count_indel_pairs)
#   6. final translocation table (filter_translocation)
#
# Config keys (see configs/target_sequence_240501/snakemake_config.yaml):
#   samples_table     : TSV with columns read1, read2, sample_name
#   primers_table     : CSV with columns name, chr, pos, length, strand
#   project_name      : batch name used under data/processed and data/final
#   genome_index      : path to the GRCm39 BWA index (external reference)
#   conda_env         : conda environment used for script rules
#   fastp_trim / fastp_stitch / primer_filter_strategy / min_mapq /
#   min_complexity / min_map_length : analysis parameters
# ============================================================

import sys
from pathlib import Path

import pandas as pd
import siuba as xb
from siuba import _

base_dir = workflow.workdir_init
if base_dir not in sys.path:
    sys.path.insert(0, base_dir)

from src.config import Location

# ============================================================
# Parameters
# ============================================================

DEFAULT_PARAMS = {
    "fastp_trim"            : True,
    "fastp_stitch"          : True,
    "primer_filter_strategy": "cover",  # strategy used for the final results
    "min_mapq"              : 0,
    "min_complexity"        : 1.5,
    "min_map_length"        : 50,
}
params = {**DEFAULT_PARAMS, **config}

use_fastp = params["fastp_trim"] or params["fastp_stitch"]

# GRCm39 BWA index (external reference data; must be provided by the user)
GENOME_INDEX = config.get("genome_index", None)
if GENOME_INDEX is None:
    raise ValueError(
        "Please provide 'genome_index' (path to the GRCm39 BWA index) in the snakemake config."
    )

# Conda environment used by the script rules
COND_ENV = config.get("conda_env", "long-term-crispr-in-liver")

# Input tables
sample_df     = pd.read_table(params["samples_table"])  # three columns: read1, read2, sample_name
processed_dir = Location.processed_data / params["project_name"]
final_dir     = Location.final_data / params["project_name"]
primers_df    = pd.read_csv(config["primers_table"], dtype={"chr": str})

# ============================================================
# Wildcards
# ============================================================

raw_read1 = lambda wildcards: str(
    sample_df
    >> xb.filter(_.sample_name == wildcards.sample_name)
    >> _.read1.iloc[0]
)
raw_read2 = lambda wildcards: str(
    sample_df
    >> xb.filter(_.sample_name == wildcards.sample_name)
    >> _.read2.iloc[0]
)

fastp_read1 = str(processed_dir / "{ sample_name }" / "fastp" / "trimmed_R1.fastq.gz")
fastp_read2 = str(processed_dir / "{ sample_name }" / "fastp" / "trimmed_R2.fastq.gz")

# Feed the BWA rule with trimmed reads when fastp is enabled
if use_fastp:
    processed_read1 = fastp_read1
    processed_read2 = fastp_read2
else:
    processed_read1 = raw_read1
    processed_read2 = raw_read2

filter_strategy = params["primer_filter_strategy"]

# ============================================================
# Rules
# ============================================================

rule all:
    input:
        [
            str(processed_dir / sample_name / "mapping_stat.tsv")
            for sample_name in sample_df.sample_name
        ],
        [
            str(processed_dir / sample_name / primer_name / f"translocation_stat__{filter_strategy}.tsv")
            for sample_name in sample_df.sample_name
            for primer_name in primers_df.name
        ],
        [
            str(processed_dir / sample_name / primer_name / f"indel_stat__{filter_strategy}.tsv")
            for sample_name in sample_df.sample_name
            for primer_name in primers_df.name
        ],
        [
            final_dir / sample_name / f"translocation_stat__{filter_strategy}.filtered.tsv"
            for sample_name in sample_df.sample_name
        ]
# ---

rule fastp:
    input:
        R1 = raw_read1,
        R2 = raw_read2
    output:
        R1   = fastp_read1,
        R2   = fastp_read2,
        json = str(processed_dir / "{ sample_name }" / "fastp" / "fastp.json"),
        html = str(processed_dir / "{ sample_name }" / "fastp" / "fastp.html")
    threads: 16
    shell:
        '''
        fastp \
        -i {input.R1}\
        -I {input.R2}\
        -o {output.R1}\
        -O {output.R2}\
        --thread {threads}\
        -j {output.json}\
        -h {output.html}
        '''
# ---

rule bwa_mapping:
    input:
        refgenome = GENOME_INDEX,
        R1        = processed_read1,
        R2        = processed_read2
    output:
        bam       = str(processed_dir / "{ sample_name }" / "mapped.bam"),
        bam_index = str(processed_dir / "{ sample_name }" / "mapped.bam.bai")
    threads: 16
    shell:
        '''
        bwa mem -Y -a -t {threads} {input.refgenome} {input.R1} {input.R2} | samtools sort -@ {threads} -o {output.bam} && samtools index -@ {threads} {output.bam}
        '''
# ---

rule count_mapped_reads:
    input:
        bam = str(processed_dir / "{ sample_name }" / "mapped.bam")
    output:
        bam_stat = str(processed_dir / "{ sample_name }" / "mapping_stat.tsv")
    conda: COND_ENV
    params:
        min_mapq = params["min_mapq"]
    script:
        "scripts/count_mapped_reads.py"
# ---

rule find_intarget_reads:
    input:
        bam           = str(processed_dir / "{ sample_name }" / "mapped.bam"),
        primers_table = config["primers_table"]
    output:
        qname_list = str(processed_dir / "{sample_name}" / "{primer_name}" / f"intarget_qnames__{filter_strategy}.txt"),
        bam_stat   = str(processed_dir / "{sample_name}" / "{primer_name}" / f"intarget_stat__{filter_strategy}.tsv")
    conda: COND_ENV
    params:
        strategy    = filter_strategy,
        primer_name = "{primer_name}"
    script:
        "scripts/find_intarget_reads.py"
# ---

rule filter_intarget_reads:
    input:
        bam        = str(processed_dir / "{sample_name}" / "mapped.bam"),
        qname_list = str(processed_dir / "{sample_name}" / "{primer_name}" / f"intarget_qnames__{filter_strategy}.txt")
    output:
        tmp_bam = str(processed_dir / "{sample_name}" / "{primer_name}" / f"intarget__{filter_strategy}.filtered.bam"),
        bam     = str(processed_dir / "{sample_name}" / "{primer_name}" / f"intarget__{filter_strategy}.sorted.bam")
    threads: 16
    shell:
        '''
        samtools view -@ threads -N {input.qname_list} -o {output.tmp_bam} {input.bam} && samtools sort -n -@ threads -o {output.bam} {output.tmp_bam}
        '''
# ---

rule check_translocation:
    input:
        bam           = str(processed_dir / "{sample_name}" / "{primer_name}" / f"intarget__{filter_strategy}.sorted.bam"),
        primers_table = config["primers_table"]
    output:
        translocation_stat = str(processed_dir / "{sample_name}" / "{primer_name}" / f"translocation_stat__{filter_strategy}.tsv")
    conda: COND_ENV
    params:
        primer_name    = "{primer_name}",
        min_mapq       = params["min_mapq"],
        min_complexity = params["min_complexity"],
        min_map_length = params["min_map_length"]
    script:
        "scripts/check_translocations.py"
# ---

# Short-indel counting on the on-target reads
rule count_indel_pairs:
    input:
        bam           = str(processed_dir / "{sample_name}" / "{primer_name}" / f"intarget__{filter_strategy}.sorted.bam"),
        primers_table = config["primers_table"]
    output:
        indel_stat = str(processed_dir / "{sample_name}" / "{primer_name}" / f"indel_stat__{filter_strategy}.tsv")
    conda: COND_ENV
    params:
        primer_name    = "{primer_name}",
        min_mapq       = params["min_mapq"],
        min_complexity = params["min_complexity"],
        min_map_length = params["min_map_length"]
    script:
        "scripts/count_indel_pairs.py"
# ---

# Aggregate the per-primer translocation calls into the final table
rule filter_translocation:
    input:
        lambda wildcards: [
            processed_dir / wildcards.sample_name / primer_name / f"translocation_stat__{filter_strategy}.tsv"
            for primer_name in primers_df.name
        ]
    output:
        final_dir / "{sample_name}" / f"translocation_stat__{filter_strategy}.filtered.tsv"
    conda: COND_ENV
    script:
        "scripts/filter_translocation.R"
