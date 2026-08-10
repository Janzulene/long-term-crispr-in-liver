import pandas as pd
from snakemake.script import snakemake

from src.detect_translocation import detect_translocation

def main(
        bam_file          : str,
        primers_table_file: str,
        primer_name       : str,
        output_stats_file : str,
        min_mapq          : int   = 0,
        min_complexity    : float = 1.5,
        min_map_length    : int   = 50
    ) -> None:

    primers_table = pd.read_csv(primers_table_file, dtype={"chr":str})
    _primer       = primers_table[primers_table["name"] == primer_name].reset_index(drop=True)
    primer        = _primer.iloc[0]

    translocation_stat = detect_translocation(
        primer_chr     = primer.chr,     # type:ignore
        primer_pos     = primer.pos,     # type:ignore
        primer_length  = primer.length,  # type:ignore
        primer_strand  = primer.strand,  # type:ignore
        bam_file_path  = bam_file,
        min_mapq       = min_mapq,
        min_complexity = min_complexity,
        min_map_length = min_map_length
    )

    translocation_stat.to_csv(output_stats_file, sep="\t", index=False)

main(
    bam_file           = snakemake.input.bam,                 # type:ignore
    primers_table_file = snakemake.input.primers_table,       # type:ignore
    primer_name        = snakemake.params.primer_name,        # type:ignore
    output_stats_file  = snakemake.output.translocation_stat, # type:ignore
    min_mapq           = snakemake.params.min_mapq,           # type:ignore
    min_complexity     = snakemake.params.min_complexity,     # type:ignore
    min_map_length     = snakemake.params.min_map_length      # type:ignore
)
    