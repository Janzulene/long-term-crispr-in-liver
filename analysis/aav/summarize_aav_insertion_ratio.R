# ══════════════════════════════════════════════════════════════════════════════
# Summarize AAV insertion ratios into the final tables.
#
# Collects the per-sample AAV insertion ratios produced by
# target_analyse_aav.smk and writes:
#   - data/final/targetsequence_20240501/aav_insertion_ratio.tsv
#     (overall ratio, using the "con" reference sgRNA)
#   - data/final/targetsequence_20240501/aav_insertion_ratio.2.tsv
#     (per-sgRNA ratio relative to the editing-evidence reads)
#
# Article references:
#   - Fig 2K (AAV integration detection): upstream data tables
#
# Usage:
#   Rscript analysis/aav/summarize_aav_insertion_ratio.R
# ══════════════════════════════════════════════════════════════════════════════

library(conflicted)

conflicts_prefer(dplyr::filter)

library(tidyverse)

AAV_RESPATH <- "data/processed/targetsequence_20240501_aav/mapping"
FINAL_DIR   <- "data/final/targetsequence_20240501"

# ---------- method 1: overall ratio ----------

aav_ratio_df <- map_dfr(
    list.files(AAV_RESPATH, pattern = "*.aav_insertion_ratio.txt", recursive = TRUE, full.names = TRUE),
    function(file_path) {
        sample_name <- dirname(file_path) |> basename()
        sgRNA_name  <- basename(file_path) |> str_remove(".aav_insertion_ratio.txt")
        res_df <- read_tsv(file_path, progress = FALSE, show_col_types = FALSE) |>
            mutate(
                sample_name    = sample_name,
                ref_sgRNA_name = sgRNA_name
            )
    }
)

# Using "con" as the reference sgRNA
dir.create(FINAL_DIR, showWarnings = FALSE, recursive = TRUE)
aav_ratio_df |>
    filter(ref_sgRNA_name == "con") |>
    select(-ref_sgRNA_name) |>
    write_tsv(file.path(FINAL_DIR, "aav_insertion_ratio.tsv"))

p1 <- (
    aav_ratio_df
    |> filter(ref_sgRNA_name == "con")
    |> ggplot(aes(x = sample_name, y = aav_insertion_ratio))
    + geom_col(fill = "steelblue", width = 0.6)
    + labs(y = "AAV Insertion Ratio", x = "Sample Name")
    + scale_y_continuous(expand = c(0, 0), labels = scales::percent_format())
    + theme_classic(base_size = 18)
)

# ---------- method 2: per-region ratio ----------

aav_ratio2_df <- map_dfr(
    list.files(AAV_RESPATH, pattern = "*.aav_insertion_ratio.2.txt", recursive = TRUE, full.names = TRUE),
    function(file_path) {
        sample_name <- dirname(file_path) |> basename()
        sgRNA_name  <- basename(file_path) |> str_remove(".aav_insertion_ratio.2.txt")
        res_df <- read_tsv(file_path, progress = FALSE, show_col_types = FALSE) |>
            mutate(
                sample_name = sample_name
            )
    }
)

aav_ratio2_df |> write_tsv(file.path(FINAL_DIR, "aav_insertion_ratio.2.tsv"))

p2 <- (
    aav_ratio2_df
    |> ggplot(aes(x = sample_name, y = aav_insertion_ratio))
    + geom_col(fill = "steelblue", width = 0.6)
    + labs(y = "AAV Insertion Ratio", x = "Sample Name")
    + scale_y_continuous(expand = c(0, 0, 0.05, 0), labels = scales::percent_format())
    + theme_bw(base_size = 18)
    + facet_wrap("~sgRNA")
    + theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
    )
)

dir.create("reports/figures/targetsequence_20240501", showWarnings = FALSE, recursive = TRUE)
ggsave("reports/figures/targetsequence_20240501/aav_insertion_ratio.svg",  p1, width = 8, height = 4)
ggsave("reports/figures/targetsequence_20240501/aav_insertion_ratio.2.svg", p2, width = 8, height = 4)
