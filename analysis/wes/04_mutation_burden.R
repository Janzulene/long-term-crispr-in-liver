# ══════════════════════════════════════════════════════════════════════════════
# WES somatic mutation analysis — step 4/7: mutation burden
#
# Counts distinct mutations (CHROM/POS/REF/ALT) per sample for the raw
# calls and for the hard-filtered somatic calls, and saves the
# hard-filtered burden table.
#
# Article references:
#   - Table S4: mutation burden
#     (data/final/WES/mutation_burden.tsv)
#
# Usage:
#   Rscript analysis/wes/04_mutation_burden.R
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

library(RColorBrewer)
library(magrittr)
library(tidyverse)

# Shared path conventions (see src/R/paths.R)
source("src/R/paths.R")

SAMPLE_ORDER <- c(
    str_c("PBS_Rep", seq(3)),
    str_c("NT_Rep", seq(3)),
    str_c("Angptl3G4_Rep", seq(5)),
    str_c("Pcsk9G1_Rep", seq(5)),
    str_c("Pcsk9G3_Rep", seq(5))
)

SAMPLE_GROUP <- c(
    rep("PBS_Rep", 3),
    rep("NT_Rep", 3),
    rep("Angptl3G4_Rep", 5),
    rep("Pcsk9G1_Rep", 5),
    rep("Pcsk9G3_Rep", 5)
)
names(SAMPLE_GROUP) <- SAMPLE_ORDER

samples             <- readRDS(file.path(WES_PROCESSED, "samples.rds"))
samples.somatic_hard <- readRDS(file.path(WES_PROCESSED, "samples.somatic_hard.rds"))

# ══════════════════════════════════════════════════════════════════════════════
# Burden functions
# ══════════════════════════════════════════════════════════════════════════════

# Count distinct mutations per sample
calc_mut_burden.vcf_df <- function(vcf_df) {
    vcf_df |> n_distinct("CHROM", "POS", "REF", "ALT")
}

calc_mut_burden.samples <- function(samples) {
    map(samples, calc_mut_burden.vcf_df) |>
        unlist() |>
        enframe(name = "sample", value = "mut_burden")
}

plot_mut_burden <- function(samples, title = "Mutation Burden") {
    calc_mut_burden.samples(samples) |>
        mutate(
            group  = SAMPLE_GROUP[sample],
            sample = factor(sample, levels = SAMPLE_ORDER)
        ) |>
        ggplot() +
        aes(sample, mut_burden, fill = group) +
        geom_col(width = 0.6) +
        labs(x = NULL, y = "Mutation Burden", title = title) +
        scale_y_continuous(expand = c(0, 0, 0.1, 0)) +
        scale_fill_brewer(palette = "Accent") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

# ══════════════════════════════════════════════════════════════════════════════
# Burden plots + table
# ══════════════════════════════════════════════════════════════════════════════

dir.create(WES_FIGURES, showWarnings = FALSE, recursive = TRUE)

# All mutations (raw unfiltered)
ggsave(
    file.path(WES_FIGURES, "mutation_burden.raw.svg"),
    plot_mut_burden(samples),
    width = 8, height = 4
)

# Hard-filtered somatic only
ggsave(
    file.path(WES_FIGURES, "mutation_burden.somatic_hard.svg"),
    plot_mut_burden(samples.somatic_hard, title = "Mutation Burden (Somatic Hard)"),
    width = 8, height = 4
)

# Save burden table
dir.create(WES_FINAL, showWarnings = FALSE, recursive = TRUE)
mut_burden_hard <- calc_mut_burden.samples(samples.somatic_hard)
write_tsv(mut_burden_hard, file.path(WES_FINAL, "mutation_burden.tsv"))
