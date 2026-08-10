# Aggregate per-primer translocation calls into the final table
# (snakemake script rule; paths are relative to the repository root).
library(conflicted)

conflicts_prefer(dplyr::filter)

library(tidyverse)

TARGET_DATA_PATH <- file.path("data/processed", snakemake@config[["project_name"]])

# Read translocation calls for one sample and one target region,
# pooling the forward and reverse primers of the region.
.read_trans <- function(sample_name, primer_prefix) {
    primer_names <- str_c(primer_prefix, c("-f", "-r"))

    trans_res <- map(
        primer_names,
        function(primer_name) {
            trans_res_path <- glue::glue("{TARGET_DATA_PATH}/{sample_name}/{primer_name}/translocation_stat__cover.tsv")
            if (file.size(trans_res_path) > 1) {
                read_tsv(
                    trans_res_path,
                    col_types = cols(offtarget_chr = col_character(), intarget_chr = col_character())
                )
            } else {
                tibble()
            }
        }
    ) |> bind_rows()

    if (nrow(trans_res) > 0) {
        sum_trans_res <- (
            trans_res
            |> filter(
                first_query_breakpoint - second_query_breakpoint <= 0,
                breakpoint_location != "between"
            )
            |> group_by(across(-c(qname, sequence, first_query_breakpoint, second_query_breakpoint)))
            |> summarise(n = n(), .groups = "drop")
            |> select(n, everything())
        )
    } else {
        sum_trans_res <- tibble()
    }

    sum_trans_res
}

region_list <- c("angptl3-g4", "pcsk9-g1", "pcsk9-g3")

sample_name <- snakemake@wildcards[["sample_name"]]
output_file <- snakemake@output[[1]]

trans_df <- (
    lapply(
        region_list,
        \(region) {
            .read_trans(sample_name, region) |> mutate(region = region)
        }
    )
    |> bind_rows()
)

trans_df |> write_tsv(output_file)
