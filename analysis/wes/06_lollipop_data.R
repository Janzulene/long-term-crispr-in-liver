# ══════════════════════════════════════════════════════════════════════════════
# WES somatic mutation analysis — step 6/7: lollipop data assembly
#
# Assembles the gene-level and CDS-level mutation detail tables used by
# the standalone lollipop plot script (analysis/wes/07_plot_lollipop.R).
#
# Article references:
#   - Lollipop figure: genomic distribution and allele frequencies of
#     somatic mutations in representative genes (Angptl3, Pcsk9, Fancm,
#     Atm, Nsd1, Trp53bp1, Atr, Fancd2)
#
# Outputs:
#   - data/final/WES/mutation_detail.gene.tsv
#   - data/final/WES/mutation_detail.cds.tsv
#   - data/final/WES/gene_widths.tsv
#
# Usage:
#   Rscript analysis/wes/06_lollipop_data.R
#   Rscript analysis/wes/07_plot_lollipop.R   # standalone plotting
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

samples.exploded <- readRDS(file.path(WES_PROCESSED, "samples.exploded.rds"))
reference        <- readRDS(file.path(WES_PROCESSED, "reference.rds"))

grcm39_gene_table <- reference$grcm39_gene_table

# ══════════════════════════════════════════════════════════════════════════════
# 7.1 Data assembly
# ══════════════════════════════════════════════════════════════════════════════

goi_list <- c(
    "Angptl3", "Pcsk9",
    "Fancm", "Atm", "Nsd1", "Trp53bp1", "Atr", "Fancd2"
)

# --- Format amino acid label ---
.format_aa_label <- function(aa, pos) {
    if_else(
        is.na(aa) | aa == "" | is.na(pos),
        NA_character_,
        str_c("p.", str_sub(aa, 1, 1), pos, str_sub(aa, 3, 3))
    )
}

# --- Build gene-level data (strand-aware gene_pos) ---
all_gene_df <- map_dfr(names(samples.exploded), \(sn) {
    samples.exploded[[sn]] |>
        filter(vep_gene_name %in% goi_list) |>
        left_join(
            grcm39_gene_table |> select(gene_name, seqnames, start, end, strand, width),
            by = c("vep_gene_name" = "gene_name")
        ) |>
        mutate(
            gene_pos = if_else(strand == "+", POS - start + 1, end - POS + 1)
        ) |>
        filter(gene_pos >= 1, gene_pos <= width) |>
        mutate(
            aa_label = .format_aa_label(vep_amino_acids, protein_pos),
            sample_name = sn,
            RPA = as.character(RPA)  # keep RPA column type consistent across samples when binding rows
        )
})

# --- Build CDS-level data (non-NA protein_pos only) ---
all_cds_df <- all_gene_df |>
    filter(!is.na(protein_pos)) |>
    select(
        sample_name,
        gene_name = vep_gene_name,
        CHROM, POS, REF, ALT,
        allele_frequency = exp_AF,
        IMPACT,
        consequence,
        primary_consequence,
        protein_pos,
        vep_amino_acids,
        aa_label
    )

# ══════════════════════════════════════════════════════════════════════════════
# 7.2 Save TSV tables
# ══════════════════════════════════════════════════════════════════════════════

dir.create(WES_FINAL, showWarnings = FALSE, recursive = TRUE)

# Gene widths (used by 07_plot_lollipop.R)
grcm39_gene_table |>
    filter(gene_name %in% goi_list) |>
    distinct(gene_name, width) |>
    write_tsv(file.path(WES_FINAL, "gene_widths.tsv"))

all_gene_df |>
    select(
        sample_name, gene_name = vep_gene_name, CHROM, POS, REF, ALT,
        allele_frequency = exp_AF, IMPACT,
        consequence, primary_consequence,
        gene_pos, amino_acids = vep_amino_acids, aa_label
    ) |>
    write_tsv(file.path(WES_FINAL, "mutation_detail.gene.tsv"))

all_cds_df |>
    write_tsv(file.path(WES_FINAL, "mutation_detail.cds.tsv"))
