# ══════════════════════════════════════════════════════════════════════════════
# WES somatic mutation analysis — step 2/7: somatic filtering
#
# Two-stage filtering of the Mutect2 calls:
#   1. Initial somatic filter (allele-frequency + quality + germline
#      removal), an exact reimplementation of the original
#      filter_somatic.vcf_df logic
#   2. Hard filter: remove mutations shared across multiple samples
#      (likely systematic sequencing artifacts), while protecting
#      mutations in the target regions (Angptl3 / Pcsk9)
#
# Article references:
#   - Table S5 / Fig 2I / lollipop figure: upstream filtering step
#     (produces samples.somatic_hard)
#
# Usage:
#   Rscript analysis/wes/02_somatic_filter.R
# ══════════════════════════════════════════════════════════════════════════════

library(conflicted)

conflicts_prefer(base::intersect)
conflicts_prefer(base::rank)
conflicts_prefer(base::setdiff)
conflicts_prefer(base::union)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::slice)
conflicts_prefer(magrittr::set_names)
conflicts_prefer(purrr::reduce)
conflicts_prefer(tidyr::unite)

library(glue)
library(magrittr)
library(tidyverse)

# Shared path conventions (see src/R/paths.R)
source("src/R/paths.R")

source("src/R/utils.R")  # .add_mut_id()

MIN_AF <- 0.1

# Helper: list of element sets -> binary matrix (rows = elements)
.list_to_BinaryMatrix <- function(upset_list) {
    unique_elements <- unique(unlist(upset_list))
    binary_matrix <- lapply(
        upset_list,
        \(x) { unique_elements %in% x }
    ) |> as.data.frame()
    rownames(binary_matrix) <- unique_elements
    binary_matrix
}

samples <- readRDS(file.path(WES_PROCESSED, "samples.rds"))

# ══════════════════════════════════════════════════════════════════════════════
# 2.1 Initial somatic filter
# ══════════════════════════════════════════════════════════════════════════════

.filter_somatic <- function(df, min_af = MIN_AF) {
    # Rescue parameters: keep low-AF somatic if the control AF is negligible
    RESCUE_RATIO       <- 4
    RESCUE_CON_AF_MAX  <- 0.01

    # --- below: exp_AF < min_af, rescue if signal >> noise ---
    below_muts <- df |>
        filter(exp_AF < min_af) |>
        mutate(
            con_AF = if_else(
                !is.na(con_AF) &
                    (exp_AF >= RESCUE_RATIO * con_AF) &
                    (con_AF <= RESCUE_CON_AF_MAX),
                NA,
                con_AF
            )
        ) |>
        filter(is.na(con_AF))

    # --- Germline positions (by CHROM+POS only, ignoring REF/ALT) ---
    germline_muts <- df |>
        filter(con_AF > min_af) |>
        select(CHROM, POS)

    # --- above: exp_AF >= min_af, keep where con_AF < min_af ---
    above_muts <- df |>
        mutate(
            exp_AF = if_else(exp_AF >= min_af, exp_AF, NA),
            con_AF = if_else(con_AF >= min_af, con_AF, NA)
        ) |>
        filter(
            is.na(con_AF),
            !is.na(exp_AF)
        )

    # --- Combine + quality filter + anti_join germline ---
    bind_rows(below_muts, above_muts) |>
        filter(
            exp_DP > 20,
            MBQ   > 20,
            MMQ   > 30
        ) |>
        anti_join(germline_muts, by = c("CHROM", "POS"))
}

samples.somatic <- map(samples, .filter_somatic)

# ══════════════════════════════════════════════════════════════════════════════
# 2.2 Hard filter (multi-sample shared removal)
# ══════════════════════════════════════════════════════════════════════════════

# Identify multi-sample shared mutations
get_uniq_mut_id <- function(df) {
    .add_mut_id(df) |> pull(mut_id) |> unique()
}

mut_ids <- map(samples.somatic, get_uniq_mut_id)
bm_df   <- .list_to_BinaryMatrix(mut_ids)

mut_occur_df <- bm_df |>
    mutate(n_occurs = rowSums(across(everything()))) |>
    filter(n_occurs > 1) |>
    rownames_to_column("mut_id")

# --- Hard filter per sample ---
.filter_hard <- function(df) {
    df <- .add_mut_id(df)

    angptl3_muts <- df |>
        filter(CHROM == 4, POS >= 98919191, POS <= 98934348) |>
        pull(mut_id)

    pcsk9_muts <- df |>
        filter(CHROM == 4, POS >= 106299526, POS <= 106321526) |>
        pull(mut_id)

    target_muts <- c(angptl3_muts, pcsk9_muts)

    bind_rows(
        # Non-target: remove shared + AF cap
        df |> filter(
            !mut_id %in% mut_occur_df$mut_id,
            !mut_id %in% target_muts,
            exp_AF < 0.3
        ),
        # Target: keep all (even if shared)
        df |> filter(mut_id %in% target_muts)
    )
}

samples.somatic_hard <- map(samples.somatic, .filter_hard)

# ══════════════════════════════════════════════════════════════════════════════
# Save intermediate results
# ══════════════════════════════════════════════════════════════════════════════

saveRDS(samples.somatic,      file.path(WES_PROCESSED, "samples.somatic.rds"))
saveRDS(samples.somatic_hard, file.path(WES_PROCESSED, "samples.somatic_hard.rds"))
