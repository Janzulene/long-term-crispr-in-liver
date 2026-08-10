# ══════════════════════════════════════════════════════════════════════════════
# Summarize translocation calls into the final table.
#
# Collects the per-sample final tables produced by the
# filter_translocation rule of target_analyse_2.smk
# (translocation_stat__cover.filtered.tsv, one row per unique
# breakpoint combination with support count `n`) and pools them into a
# single cross-sample table with a `sample` column.
#
# Article references:
#   - Fig 2L (circos plots): upstream data table for
#     analysis/translocation/03_plot_translocation_circos.R
#   - Fig 2L (breakpoint validation): input for
#     analysis/translocation/02_trans_breakpoint_genbank.py
#
# Usage:
#   python analysis/translocation/01_summarize_translocations.py
# ══════════════════════════════════════════════════════════════════════════════

import pandas as pd

FINAL_DIR = Path("data/final/targetsequence")
OUTPUT    = FINAL_DIR / "translocation_stat.tsv"

filtered_tables = sorted(FINAL_DIR.glob("*/translocation_stat__cover.filtered.tsv"))

res_list = []
for table_path in filtered_tables:
    trans_df = pd.read_table(table_path)
    if len(trans_df) == 0:
        continue
    trans_df["sample"] = table_path.parent.name
    res_list.append(trans_df)

final_stat = pd.concat(res_list, ignore_index=True)
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
final_stat.to_csv(OUTPUT, index=False, sep="\t")
print(f"wrote {OUTPUT} ({len(final_stat)} rows)")
