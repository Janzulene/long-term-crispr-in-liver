# ══════════════════════════════════════════════════════════════════════════════
# WES somatic mutation analysis — step 5/7: per-gene non-synonymous ranking
#
# Ranks genes by non-synonymous mutation count and by summed allele
# frequency (optionally scaled by gene length), computes mean-rank
# statistics across sample groups, and saves the ranking tables and
# plots (data/final/WES/rank_results/).
#
# Article references:
#   - Fig 2I: mutated genes ranked by normalized non-silent mutation
#     ratio in the Angptl3-G4 group (group "cancer" / "angptl3g4" below)
#   - Table S5A: NSgene.RatioRank.all.tsv (all edited samples vs control)
#   - Table S5B: NSgene.RatioRank.cancer.tsv (Angptl3-G4 group)
#
# Usage:
#   Rscript analysis/wes/05_gene_ranking.R
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
library(glue)
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

EXP_SAMPLES <- SAMPLE_ORDER[!str_detect(SAMPLE_ORDER, "PBS|NT")]

samples.exploded <- readRDS(file.path(WES_PROCESSED, "samples.exploded.rds"))
reference        <- readRDS(file.path(WES_PROCESSED, "reference.rds"))

grcm39_gene_table  <- reference$grcm39_gene_table
id2length          <- reference$id2length
cancer_gene_id_list <- reference$cancer_gene_id_list
cancer_gene_list   <- reference$cancer_gene_list

# Mark genes from the cancer gene list (mirrors the definition in 01_load_data.R)
annotate_cancer <- function(df) {
    df |> mutate(cancer_related = vep_gene_id %in% cancer_gene_id_list)
}

tilt_x_axis_text <- theme(
    axis.text.x = element_text(angle = 45, color = "black", hjust = 1, vjust = 1)
)

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

# ══════════════════════════════════════════════════════════════════════════════
# 6.1 Count & ratio functions
# ══════════════════════════════════════════════════════════════════════════════

#' Count non-synonymous mutations per gene (deduplicated by .row_id).
#' Only counts HIGH + MODERATE impact, mut_status = "acquired".
calc_mut_status_count <- function(df) {
    df |>
        filter(IMPACT %in% c("HIGH", "MODERATE")) |>
        annotate_cancer() |>
        group_by(vep_gene_id, vep_gene_name, cancer_related) |>
        summarise(
            n          = n_distinct(.row_id),
            mut_status = "acquired",
            .groups    = "drop"
        ) |>
        arrange(desc(n)) |>
        mutate(rank = rank(-n, ties.method = "min"))
}

#' Sum exp_AF per gene as ratio-based ranking.
#' Optionally scale by gene length.
calc_mut_status_ratio <- function(df, scale_by_gene_length = FALSE) {
    rank_df <- df |>
        filter(IMPACT %in% c("HIGH", "MODERATE")) |>
        annotate_cancer() |>
        distinct(.row_id, vep_gene_id, vep_gene_name, cancer_related, exp_AF) |>
        group_by(vep_gene_id, vep_gene_name, cancer_related) |>
        summarise(
            n          = sum(exp_AF),
            mut_status = "acquired",
            .groups    = "drop"
        )

    if (scale_by_gene_length) {
        rank_df <- rank_df |>
            left_join(id2length, by = c("vep_gene_id" = "gene_id")) |>
            mutate(n = n / width) |>
            select(-width)
    }

    rank_df |>
        arrange(desc(n)) |>
        mutate(rank = rank(-n, ties.method = "min"))
}

# ══════════════════════════════════════════════════════════════════════════════
# 6.2 Bar plot (ranking visualization)
# ══════════════════════════════════════════════════════════════════════════════

#' Bar plot of top N genes per mut_status, coloured by gene type.
count_plot_all <- function(count_df, n_top = 20, status_list = NULL) {
    if (!is.null(status_list)) {
        count_df <- count_df |> filter(mut_status %in% status_list)
    }

    plot_df <- count_df |>
        filter(mut_status != "constant") |>
        unite(".sort_gene_id", mut_status, vep_gene_id, sep = "_", remove = FALSE) |>
        mutate(
            mut_status = factor(mut_status, levels = c("acquired", "removed", "emerging", "disappearing")),
            gene_type  = case_when(
                vep_gene_name %in% c("Pcsk9", "Angptl3") ~ "target",
                cancer_related ~ "cancer related",
                TRUE ~ "normal"
            )
        ) |>
        group_by(mut_status) |>
        arrange(mut_status, desc(n)) |>
        mutate(.sort_gene_id = factor(.sort_gene_id, levels = .sort_gene_id)) |>
        slice_max(n, n = n_top) |>
        ungroup()

    ggplot(plot_df) +
        aes(.sort_gene_id, n, fill = gene_type) +
        scale_fill_manual(
            values = c("cancer related" = "red", "normal" = "gray70", "target" = "#7fc97f")
        ) +
        geom_col(width = 0.7) +
        facet_wrap("mut_status", scales = "free", ncol = 1) +
        xlab("Gene") +
        scale_x_discrete(breaks = plot_df$.sort_gene_id, labels = plot_df$vep_gene_name) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 90))
}

# ══════════════════════════════════════════════════════════════════════════════
# 6.3 Compute all rankings
# ══════════════════════════════════════════════════════════════════════════════

# Count-based
samples.hard_count <- map(samples.exploded, calc_mut_status_count)

# Ratio-based (raw + scaled)
samples.hard_ratio        <- map(samples.exploded, calc_mut_status_ratio)
samples.hard_ratio_scaled <- map(samples.exploded, \(x) calc_mut_status_ratio(x, scale_by_gene_length = TRUE))

# Example ranking plot for one sample
dir.create(WES_FIGURES, showWarnings = FALSE, recursive = TRUE)
TOP_N <- 50
sample_name <- "Angptl3G4_Rep5"

ggsave(
    file.path(WES_FIGURES, "rank_example.count.svg"),
    samples.hard_count[[sample_name]] |>
        count_plot_all(n_top = TOP_N) +
        labs(title = glue("{sample_name}: Top {TOP_N} (count)")),
    width = 8, height = 8
)

ggsave(
    file.path(WES_FIGURES, "rank_example.ratio.svg"),
    samples.hard_ratio[[sample_name]] |>
        count_plot_all(n_top = TOP_N) +
        labs(title = glue("{sample_name}: Top {TOP_N} (ratio)")),
    width = 8, height = 8
)

# ══════════════════════════════════════════════════════════════════════════════
# 6.4 Share gene analysis
# ══════════════════════════════════════════════════════════════════════════════

# Genes appearing in top N of multiple samples
.top_n <- 50

share_gene_df <- map(
    samples.hard_ratio,
    \(x) slice_max(x, order_by = n, n = .top_n) |> pull(vep_gene_name)
) |>
    .list_to_BinaryMatrix() |>
    mutate(n_share = rowSums(across(everything()))) |>
    filter(n_share > 1) |>
    rownames_to_column("gene_name") |>
    mutate(cancer_related = gene_name %in% cancer_gene_list)

share_gene_df |>
    select(gene_name, n_share, cancer_related) |>
    arrange(desc(n_share)) |>
    filter(cancer_related)

# ══════════════════════════════════════════════════════════════════════════════
# 6.5 Cancer gene ratio per sample
# ══════════════════════════════════════════════════════════════════════════════

cancer_ratio_df <- map2(
    names(samples.hard_ratio),
    samples.hard_ratio,
    \(.name, x) {
        x |>
            slice_max(n, n = 50) |>
            count(cancer_related) |>
            mutate(ratio = n / sum(n), sample = .name) |>
            filter(cancer_related)
    }
) |> bind_rows()

ggsave(
    file.path(WES_FIGURES, "cancer_gene_ratio.top50.svg"),
    cancer_ratio_df |>
        mutate(group = SAMPLE_GROUP[sample]) |>
        ggplot() +
        aes(sample, ratio, fill = group) +
        geom_col() +
        scale_fill_brewer(palette = "Accent") +
        labs(y = "Cancer Gene Ratio", title = "Cancer Gene Ratio in Top 50 Genes") +
        theme_classic() +
        tilt_x_axis_text,
    width = 8, height = 4
)

# ══════════════════════════════════════════════════════════════════════════════
# 6.6 Mean rank (MR) computation
# ══════════════════════════════════════════════════════════════════════════════

# --- Find genes present in > 1 sample ---
share_gene_ids <- map(
    samples.hard_ratio,
    \(x) pull(x, vep_gene_id)
) |>
    .list_to_BinaryMatrix() |>
    mutate(n_share = rowSums(across(everything()))) |>
    filter(n_share > 1) |>
    rownames_to_column("gene_id")

# --- Build unified rank dataframe ---
all_ratio_df <- map(
    names(samples.hard_ratio),
    \(x) {
        samples.hard_ratio[[x]] |>
            filter(vep_gene_id %in% share_gene_ids$gene_id) |>
            mutate(sample = x)
    }
) |> bind_rows()

all_ratio_df.scaled <- map(
    names(samples.hard_ratio_scaled),
    \(x) {
        samples.hard_ratio_scaled[[x]] |>
            filter(vep_gene_id %in% share_gene_ids$gene_id) |>
            mutate(sample = x)
    }
) |> bind_rows()

# --- MR helper: compute mean rank across specified samples ---
.get_mean_rank <- function(.all_rank_df, do_filter = TRUE, filtered_samples = c()) {
    if (do_filter) {
        stopifnot(length(filtered_samples) > 0)
        .all_rank_df <- .all_rank_df |> filter(sample %in% filtered_samples)
    }

    .all_rank_df |>
        mutate(is_exp = sample %in% EXP_SAMPLES) |>
        group_by(vep_gene_id, vep_gene_name, cancer_related) |>
        summarise(
            MR            = mean(rank),
            MRR           = mean(1 / rank),
            n_occurs      = n(),
            only_exp      = all(is_exp),
            mean_sum_ratio = mean(n),
            .groups       = "drop"
        ) |>
        arrange(MR)
}

# --- Plot ranked genes ---
.plot_rank <- function(samples_rank, rank_key = "mean_sum_ratio") {
    samples_rank |>
        filter(n_occurs >= 3) |>
        mutate(
            is_target = vep_gene_name %in% c("Angptl3", "Pcsk9"),
            gene_anno = case_when(
                is_target ~ "target",
                cancer_related ~ "cancer related",
                TRUE ~ "others"
            )
        ) |>
        ggplot() +
        aes(reorder(vep_gene_name, .data[[rank_key]]), .data[[rank_key]], fill = gene_anno) +
        geom_col(width = 0.7) +
        scale_fill_manual(values = c("target" = "lightgreen", "cancer related" = "red", "others" = "gray70")) +
        coord_flip() +
        labs(x = NULL, y = rank_key)
}

# ══════════════════════════════════════════════════════════════════════════════
# 6.7 Batch save rank results
# ══════════════════════════════════════════════════════════════════════════════

res_dir <- file.path(WES_FINAL, "rank_results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

.save_result <- function(group_name, group_samples) {
    .all_rank        <- .get_mean_rank(all_ratio_df,        filtered_samples = group_samples)
    .all_rank_scaled <- .get_mean_rank(all_ratio_df.scaled, filtered_samples = group_samples)

    write_tsv(.all_rank,        file.path(res_dir, glue("NSgene.RatioRank.{group_name}.tsv")))
    write_tsv(.all_rank_scaled, file.path(res_dir, glue("NSgene.ScaledRatioRank.{group_name}.tsv")))

    p1 <- .plot_rank(.all_rank, "MRR")            + labs(title = glue("{group_name}: MRR"))
    p2 <- .plot_rank(.all_rank, "mean_sum_ratio") + labs(title = glue("{group_name}: Mean Sum Ratio"))
    p3 <- .plot_rank(.all_rank_scaled, "MRR")            + labs(title = glue("{group_name}: MRR (scaled)"))
    p4 <- .plot_rank(.all_rank_scaled, "mean_sum_ratio") + labs(title = glue("{group_name}: Mean Sum Ratio (scaled)"))

    ggsave(file.path(res_dir, glue("NSgene.RatioRank.{group_name}.MRR.png")),          p1, width = 6, height = 12)
    ggsave(file.path(res_dir, glue("NSgene.RatioRank.{group_name}.MeanSumRatio.png")), p2, width = 6, height = 12)
    ggsave(file.path(res_dir, glue("NSgene.ScaledRatioRank.{group_name}.MRR.png")),          p3, width = 6, height = 12)
    ggsave(file.path(res_dir, glue("NSgene.ScaledRatioRank.{group_name}.MeanSumRatio.png")), p4, width = 6, height = 12)
}

# Sample groupings
high_samples      <- c("Angptl3G4_Rep3", "Angptl3G4_Rep4", "Angptl3G4_Rep5", "Pcsk9G3_Rep3", "Pcsk9G3_Rep4")
angptl3g4_samples <- SAMPLE_ORDER[str_detect(SAMPLE_ORDER, "Angptl3G4")]
pcsk9g3_samples   <- SAMPLE_ORDER[str_detect(SAMPLE_ORDER, "Pcsk9G3")]
pcsk9g1_samples   <- SAMPLE_ORDER[str_detect(SAMPLE_ORDER, "Pcsk9G1")]

.save_result("all",       SAMPLE_ORDER)
.save_result("exp",       EXP_SAMPLES)
.save_result("cancer",    high_samples)
.save_result("angptl3g4", angptl3g4_samples)
.save_result("pcsk9g3",   pcsk9g3_samples)
.save_result("pcsk9g1",   pcsk9g1_samples)
