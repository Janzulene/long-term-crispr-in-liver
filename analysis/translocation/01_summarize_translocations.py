# ══════════════════════════════════════════════════════════════════════════════
# Summarize per-primer translocation calls into the final table.
#
# Reads the per-sample/per-primer translocation calls produced by
# target_analyse_2.smk (translocation_stat__cover.tsv), pools the
# forward and reverse primers of each target region, deduplicates
# read-level records to one row per unique breakpoint combination, and
# writes data/final/targetsequence_20240501/translocation_stat.tsv.
#
# Article references:
#   - Fig 2L (circos plots): upstream data table for
#     analysis/translocation/plot_translocation_circos.R
#
# Usage:
#   python analysis/translocation/01_summarize_translocations.py
# ══════════════════════════════════════════════════════════════════════════════

import sys
from pathlib import Path

import pandas as pd
import siuba as xb
from siuba import _
from siuba.dply.vector import row_number

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

SAMPLES_TABLE = "configs/target_sequence_240501/samples.tsv"
PRIMERS_TABLE = "configs/target_sequence_240501/target_primers.csv"

TRANS_DATA = "data/processed/targetsequence_20240501"
OUTPUT     = "data/final/targetsequence_20240501/translocation_stat.tsv"

samples = pd.read_table(SAMPLES_TABLE)
primers = pd.read_csv(PRIMERS_TABLE, dtype={"chr": str})

res_list = []
for sample_name in samples.sample_name:
    for primer_name in primers.name:
        trans_path = (
            f"{TRANS_DATA}/{sample_name}/{primer_name}/"
            "translocation_stat__cover.tsv"
        )
        try:
            trans_df = pd.read_table(trans_path)
        except pd.errors.EmptyDataError:
            continue

        # Grouped counts per unique breakpoint combination
        a = (
            trans_df
            >> xb.mutate(
                microhomo_len = _.first_query_breakpoint - _.second_query_breakpoint,
                sample        = sample_name,
                primer        = primer_name
            )
            >> xb.count(
                _.intarget_chr,
                _.intarget_strand,
                _.intarget_breakpoint_pos,
                _.intarget_breakpoint_loc,
                _.offtarget_chr,
                _.offtarget_strand,
                _.offtarget_breakpoint_pos,
                _.offtarget_breakpoint_loc,
                _.microhomo_len,
                _.breakpoint_location,
                _.intarget_first
            )
            >> xb.filter(_.microhomo_len <= 0, _.breakpoint_location != "between")
            >> xb.ungroup()
        )
        if len(a) == 0:
            continue

        # One representative record per breakpoint combination
        # (keeps sequence and strand annotation of the first read)
        b = (
            trans_df
            >> xb.mutate(
                microhomo_len = _.first_query_breakpoint - _.second_query_breakpoint,
                sample        = sample_name,
                primer        = primer_name
            )
            >> xb.group_by(
                _.intarget_chr,
                _.intarget_breakpoint_pos,
                _.offtarget_chr,
                _.offtarget_strand,
                _.offtarget_breakpoint_pos,
                _.microhomo_len,
                _.breakpoint_location,
                _.sample,
                _.primer
            )
            >> xb.mutate(row_num = row_number(_))
            >> xb.ungroup()
            >> xb.filter(_.row_num == 1)
            >> xb.select(-_.row_num, -_.qname, -_.intarget_first)
        )

        res_list.append(xb.left_join(a, b))

final_stat = pd.concat(res_list)
Path(OUTPUT).parent.mkdir(parents=True, exist_ok=True)
final_stat.to_csv(OUTPUT, index=False, sep="\t")
print(f"wrote {OUTPUT} ({len(final_stat)} rows)")
