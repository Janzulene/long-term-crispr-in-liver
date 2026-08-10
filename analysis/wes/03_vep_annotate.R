# ══════════════════════════════════════════════════════════════════════════════
# WES2 somatic mutation analysis — step 3/6: VEP annotation
#
# Explodes the VEP CSQ field of the hard-filtered somatic mutations,
# keeping a single transcript per gene (prefer CANONICAL, then
# protein_coding biotype, then highest IMPACT).
#
# Article references:
#   - Table S5 / Fig 2I / lollipop figure: upstream annotation step
#
# Usage:
#   Rscript analysis/wes/03_vep_annotate.R
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

WES_PROCESSED <- "data/processed/WES_20250105_analysis"

IMPACT_ORDER <- c("HIGH", "MODERATE", "LOW", "MODIFIER")

samples.somatic_hard <- readRDS(file.path(WES_PROCESSED, "samples.somatic_hard.rds"))

# ══════════════════════════════════════════════════════════════════════════════
# explode_csq_by_gene() — one transcript per gene
# ══════════════════════════════════════════════════════════════════════════════

#' Explode VEP CSQ field, deduplicate to one transcript per gene.
#'
#' Steps:
#' 1. Add .row_id for mutation-level deduplication
#' 2. separate_rows(CSQ, sep = ",") to explode by transcript
#' 3. Parse CSQ pipe-delimited fields
#' 4. Add primary_consequence = first term before "&"
#' 5. Filter: keep only gene-level annotations (vep_gene_id != "")
#' 6. One transcript per gene: prefer CANONICAL=YES, then protein_coding, then highest IMPACT
#' 7. Remove CANONICAL column, keep .row_id
#' 8. Extract protein_pos as integer from Protein_position field
explode_csq_by_gene <- function(vcf_df) {
    vcf_df |>
        mutate(.row_id = row_number()) |>
        separate_rows(CSQ, sep = ",") |>
        mutate(
            vep_allele          = str_split_i(CSQ, fixed("|"), 1),
            consequence         = str_split_i(CSQ, fixed("|"), 2),
            IMPACT              = str_split_i(CSQ, fixed("|"), 3),
            vep_gene_name       = str_split_i(CSQ, fixed("|"), 4),
            vep_gene_id         = str_split_i(CSQ, fixed("|"), 5),
            feature_type        = str_split_i(CSQ, fixed("|"), 6),
            biotype             = str_split_i(CSQ, fixed("|"), 8),
            protein_pos_raw     = str_split_i(CSQ, fixed("|"), 15),
            vep_amino_acids     = str_split_i(CSQ, fixed("|"), 16),
            variant_class       = str_split_i(CSQ, fixed("|"), 22),
            CANONICAL           = str_split_i(CSQ, fixed("|"), 25),
            SIFT                = str_split_i(CSQ, fixed("|"), 37),
            primary_consequence = str_split_i(consequence, fixed("&"), 1)
        ) |>
        filter(vep_gene_id != "") |>
        # One transcript per gene per mutation
        group_by(.row_id, vep_gene_id) |>
        arrange(
            desc(CANONICAL == "YES"),
            desc(biotype == "protein_coding"),
            factor(IMPACT, levels = IMPACT_ORDER, ordered = TRUE)
        ) |>
        slice(1) |>
        ungroup() |>
        select(-CANONICAL) |>
        mutate(protein_pos = str_extract(protein_pos_raw, "^[0-9]+") |> as.integer()) |>
        select(-protein_pos_raw)
}

# ══════════════════════════════════════════════════════════════════════════════
# Apply to the hard-filtered somatic mutations
# ══════════════════════════════════════════════════════════════════════════════

samples.exploded <- map(samples.somatic_hard, explode_csq_by_gene)

saveRDS(samples.exploded, file.path(WES_PROCESSED, "samples.exploded.rds"))
