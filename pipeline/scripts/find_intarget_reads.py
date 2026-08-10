from typing import Literal

import numpy as np
import pandas as pd
from snakemake.script import snakemake

from src.filter_intarget import find_intarget_reads

def main(
        bam_file          : str,
        primers_table_file: str,
        primer_name       : str,
        strategy          : Literal["only_perfect_read1", "only_perfect", "cover", "start_within4"],
        output_stats_file : str,
        output_qnames_file: str
    ) -> None:

    primers_table = pd.read_csv(primers_table_file, dtype={"chr":str})
    _primer       = primers_table[primers_table["name"] == primer_name].reset_index(drop=True)
    primer        = _primer.iloc[0]

    filter_stat, qnames = find_intarget_reads(
        primer_chr    = primer.chr,    # type:ignore
        primer_pos    = primer.pos,    # type:ignore
        primer_length = primer.length, # type:ignore
        primer_strand = primer.strand, # type:ignore
        bam_file_path = bam_file,
        strategy      = strategy
    )

    # save results
    filter_stat = pd.DataFrame([filter_stat])
    pd.concat([_primer, filter_stat], axis=1).to_csv(output_stats_file, sep="\t", index=False)

    qnames = np.unique(qnames)
    np.savetxt(output_qnames_file, qnames, fmt="%s")

main(
    bam_file           = snakemake.input.bam,           # type:ignore
    primers_table_file = snakemake.input.primers_table, # type:ignore
    primer_name        = snakemake.params.primer_name,  # type:ignore
    strategy           = snakemake.params.strategy,     # type:ignore
    output_stats_file  = snakemake.output.bam_stat,     # type:ignore
    output_qnames_file = snakemake.output.qname_list    # type:ignore
)