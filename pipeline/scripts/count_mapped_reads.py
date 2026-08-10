from snakemake.script import snakemake

from src.utils import count_total_reads

def main(
        bam_file_path    : str,
        min_mapq         : int,
        output_stats_file: str
    ) -> None:
    read1_count, read2_count = count_total_reads(bam_file_path, min_mapq)
    with open(output_stats_file, "w") as f:
        f.write(f"read1_count\tread2_count\ttotal_count\n")
        f.write(f"{read1_count}\t{read2_count}\t{read1_count + read2_count}\n")

main(snakemake.input.bam, snakemake.params.min_mapq, snakemake.output.bam_stat)
