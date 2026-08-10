# ══════════════════════════════════════════════════════════════════════════════
# WES2 somatic mutation analysis — step 1/6: load data
#
# Loads sequencing depth, reference annotations (GRCm39 GTF, cancer gene
# list) and the Mutect2 + VEP annotated variant tables, then saves them
# as RDS files consumed by the downstream analysis steps.
#
# Article references:
#   - Table S4: sequencing depth (data/final/WES_20250105/sequencing_depth.tsv)
#   - Fig 2I / Table S5 / lollipop figure: upstream data (samples.rds)
#
# External data (not shipped with the repository):
#   - GRCm39 gene annotation GTF (GENCODE vM33 / Ensembl 110) — set
#     GRCm39_GTF_PATH below
#   - cancer gene list: data/raw/WES_20240323/cancer_list.csv
#   - Mutect2 + VEP annotated VCF tables (sarek 3.4.2 output):
#     data/processed/WES2_paired_annotate/annotation/mutect2/
#   - per-sample mean depth tables (samtools depth):
#     data/processed/WES2_paired/preprocessing/markduplicates/
#
# Usage:
#   Rscript analysis/wes/01_load_data.R
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
library(assertthat)
library(glue)
library(magrittr)
library(patchwork)
library(rtracklayer)
library(tidyverse)

# ══════════════════════════════════════════════════════════════════════════════
# Parameters
# ══════════════════════════════════════════════════════════════════════════════

# GRCm39 gene annotation GTF (GENCODE vM33 / Ensembl release 110).
# Download from https://www.ensembl.org (or GENCODE) and set the path here.
GRCm39_GTF_PATH <- "path/to/GRCm39_ensembl.gtf.gz"

# Sample groups: PBS controls, NT (untreated) controls and three edited groups
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

EXP_SAMPLES  <- SAMPLE_ORDER[!str_detect(SAMPLE_ORDER, "PBS|NT")]
IMPACT_ORDER <- c("HIGH", "MODERATE", "LOW", "MODIFIER")

MIN_AF <- 0.1

# Output directories
WES_FINAL      <- "data/final/WES_20250105"
WES_PROCESSED  <- "data/processed/WES_20250105_analysis"
FIG_PATH       <- "reports/figures/WES_20250105"

DEPTH_TABLE <- file.path("data/processed/WES2_paired/preprocessing/markduplicates")
VCF_TABLE   <- "data/processed/WES2_paired_annotate/annotation/mutect2"

# ══════════════════════════════════════════════════════════════════════════════
# 1. Sequencing depth
# ══════════════════════════════════════════════════════════════════════════════

.read_samtools_stat <- function(sample_name) {
    stat_path <- str_glue("{DEPTH_TABLE}/{sample_name}/{sample_name}.depth.txt")
    dt <- data.table::fread(
        stat_path,
        col.names = c("chr", "pos", "depth"),
        colClasses = c("character", "integer", "integer")
    )
    mean_depth <- dt[, mean(depth)]
    tibble(sample = sample_name, mean_depth = mean_depth)
}

stat_df <- map(SAMPLE_ORDER, .read_samtools_stat) |> bind_rows()

dir.create(WES_FINAL, showWarnings = FALSE, recursive = TRUE)
stat_df |> write_tsv(file.path(WES_FINAL, "sequencing_depth.tsv"))

dir.create(FIG_PATH, showWarnings = FALSE, recursive = TRUE)
p_depth <- (
    stat_df |>
        mutate(sample = factor(sample, levels = SAMPLE_ORDER)) |>
        ggplot() +
        aes(sample, mean_depth) +
        geom_col(width = 0.6, fill = "gray70", color = "gray40", linewidth = 0.3) +
        geom_hline(yintercept = 125, color = "red", linetype = "dashed", linewidth = 0.6) +
        annotate("text", x = Inf, y = 125, label = "125×", color = "red",
                 hjust = 1.1, vjust = -0.5, size = 3.5, fontface = "bold") +
        scale_y_continuous(expand = c(0, 0, 0.05, 0)) +
        labs(x = NULL, y = "Mean Depth") +
        theme_classic(base_size = 12) +
        theme(
            axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.text    = element_text(color = "gray20"),
            axis.title.y = element_text(margin = margin(r = 8)),
            axis.ticks   = element_line(color = "gray60"),
            axis.line    = element_line(color = "gray60"),
            plot.margin  = margin(10, 20, 10, 10)
        )
)
ggsave(file.path(FIG_PATH, "sequencing_depth.svg"), p_depth, width = 8, height = 4)

# ══════════════════════════════════════════════════════════════════════════════
# 2. Reference data
# ══════════════════════════════════════════════════════════════════════════════

# --- GTF gene table ---
grcm39_gtf_table <- rtracklayer::import(GRCm39_GTF_PATH) |> as_tibble()

grcm39_gene_table <- grcm39_gtf_table |>
    filter(type == "gene") |>
    select(seqnames, start, end, width, strand, gene_id, gene_type, gene_name)

# --- Cancer gene list ---
cancer_gene <- read_csv(
    "data/raw/WES_20240323/cancer_list.csv",
    col_names = "oncogene"
)

cancer_gene_table <- grcm39_gene_table |>
    mutate(
        cancer_related = str_to_upper(gene_name) %in% c(
            cancer_gene$oncogene,
            "ABRAXAS1", "BABAM2", "ELOC", "FH1", "H3C2", "H3C3",
            "MYCL", "NECTIN4", "NSD2", "NSD3", "POT1A", "PRKN",
            "TENT5C", "TRP53", "TRP53BP1", "WWTR1", "ZFP217"
        )
    )

cancer_gene_id_list <- cancer_gene_table |>
    filter(cancer_related) |>
    pull(gene_id) |>
    str_split_i(fixed("."), 1)

cancer_gene_list <- cancer_gene_table |>
    filter(cancer_related) |>
    pull(gene_name) |>
    unique()

# --- Gene length for scaling ---
id2length <- grcm39_gene_table |>
    mutate(gene_id = str_split_i(gene_id, fixed("."), 1)) |>
    group_by(gene_id) |>
    summarise(width = sum(width))

# ══════════════════════════════════════════════════════════════════════════════
# 3. Read VCF data
# ══════════════════════════════════════════════════════════════════════════════

.read_paired <- function(sample_name, con_name, res_path, line_name) {
    vcf_path <- glue::glue(
        "{res_path}/{sample_name}/{sample_name}_vs_{con_name}.mutect2.biallelic_VEP.ann.vcf.tsv"
    )
    data.table::fread(vcf_path) |>
        rename_with(
            \(x) str_replace(x, str_c(line_name, "_", sample_name, "."), "exp_")
        ) |>
        rename_with(
            \(x) str_replace(x, str_c(line_name, "_", con_name, "."), "con_")
        ) |>
        mutate(
            exp_AF = as.numeric(exp_AF),
            con_AF = as.numeric(con_AF)
        )
}

strain1_sample_names <- c(
    str_c("PBS_Rep", seq(1, 2)),
    str_c("Angptl3G4_Rep", seq(5)),
    str_c("Pcsk9G1_Rep", seq(1, 2)),
    str_c("Pcsk9G3_Rep", seq(5))
)
strain1_samples <- map(
    strain1_sample_names,
    \(x) .read_paired(x, con_name = "PBS_Rep3", res_path = VCF_TABLE, line_name = "line1")
) |> set_names(strain1_sample_names)

strain2_sample_names <- c(
    str_c("NT_Rep", seq(1, 2)),
    str_c("Pcsk9G1_Rep", seq(3, 5))
)
strain2_samples <- map(
    strain2_sample_names,
    \(x) .read_paired(x, con_name = "NT_Rep3", res_path = VCF_TABLE, line_name = "line2")
) |> set_names(strain2_sample_names)

samples <- c(strain1_samples, strain2_samples)
rm(strain1_samples, strain2_samples)

# ══════════════════════════════════════════════════════════════════════════════
# Save intermediate results
# ══════════════════════════════════════════════════════════════════════════════

dir.create(WES_PROCESSED, showWarnings = FALSE, recursive = TRUE)

saveRDS(
    list(
        grcm39_gene_table  = grcm39_gene_table,
        id2length          = id2length,
        cancer_gene_id_list = cancer_gene_id_list,
        cancer_gene_list   = cancer_gene_list
    ),
    file.path(WES_PROCESSED, "reference.rds")
)

saveRDS(samples, file.path(WES_PROCESSED, "samples.rds"))
