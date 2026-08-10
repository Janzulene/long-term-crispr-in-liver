# ══════════════════════════════════════════════════════════════════════════════
# Lollipop Plot — Somatic Mutations (ggplot2, standalone)
#
# Reads mutation_detail.{gene,cds}.tsv and produces publication-ready lollipop
# figures.  No external R files needed beyond tidyverse.
#
# Usage:
#   Rscript analysis/wes/plot_lollipop.R
#   or source in RStudio for interactive use.
# ══════════════════════════════════════════════════════════════════════════════

library(tidyverse)

# ══════════════════════════════════════════════════════════════════════════════
# Parameters — adjust these before running
# ══════════════════════════════════════════════════════════════════════════════

# Input files
DATA_DIR       <- "data/final/WES_20250105"
GENE_TSV       <- file.path(DATA_DIR, "mutation_detail.gene.tsv")
CDS_TSV        <- file.path(DATA_DIR, "mutation_detail.cds.tsv")
GENE_WIDTHS    <- file.path(DATA_DIR,  "gene_widths.tsv")

# Output
OUTPUT_DIR     <- file.path(DATA_DIR, "lollipop_plot")
GENE_PNG       <- file.path(OUTPUT_DIR, "lollipop_gene.png")
CDS_PNG        <- file.path(OUTPUT_DIR, "lollipop_cds.png")

# Genes to plot
GOI_LIST <- c(
    "Angptl3", "Pcsk9",
    "Fancm", "Atm", "Nsd1", "Trp53bp1", "Atr", "Fancd2"
)

# Impact levels to include
IMPACT_KEEP <- c("HIGH", "MODERATE")

# Samples to exclude (regex)
SAMPLE_EXCLUDE <- "NT|PBS"

# Consequence colour mapping (unmatched consequences get grey)
CONSEQUENCE_COLORS <- c(
    "missense_variant"          = "#E41A1C",
    "protein_altering_variant"  = "#377EB8",
    "frameshift_variant"        = "#FF7F00",
    "splice_region_variant"     = "#984EA3",
    "inframe_deletion"          = "#A65628",
    "inframe_insertion"         = "#F781BF",
    "stop_gained"               = "#4DAF4A"
)
FALLBACK_COLOR <- "#999999"

# Plot dimensions
PNG_WIDTH_IN  <- 12
PNG_HEIGHT_IN <- 20
PNG_DPI       <- 300

# ══════════════════════════════════════════════════════════════════════════════
# Functions
# ══════════════════════════════════════════════════════════════════════════════

#' Build a fake-boundary dataframe using the actual gene widths so every
#' gene panel spans [1, width] even if no mutations fall near the ends.
.build_fake_boundary <- function(gene_widths, gene_col = "gene_name",
                                 width_col = "width", pos_col = "gene_pos") {
    gene_widths |>
        tidyr::uncount(2) |>
        group_by(.data[[gene_col]]) |>
        mutate(
            !!pos_col := if_else(row_number() == 1L, 1L, .data[[width_col]]),
            allele_frequency = 0
        ) |>
        ungroup()
}

#' Core lollipop plot (shared by gene and CDS variants).
.plot_lollipop <- function(plot_df, fake_boundary, gene_col, pos_col,
                           consequence_col = "primary_consequence") {
    ggplot(plot_df, aes(.data[[pos_col]], allele_frequency,
                        color = .data[[consequence_col]])) +
        # Gene body line
        geom_hline(yintercept = 0, linewidth = 1.8, color = "grey50") +
        # Stem
        geom_segment(
            aes(xend = .data[[pos_col]], y = 0, yend = allele_frequency),
            linewidth = 0.5, alpha = 0.35
        ) +
        # Head
        geom_point(size = 2.5, alpha = 0.85) +
        # Force full range
        geom_blank(
            aes(x = .data[[pos_col]], y = allele_frequency),
            data = fake_boundary,
            inherit.aes = FALSE
        ) +
        # Colour
        scale_color_manual(
            values = CONSEQUENCE_COLORS,
            na.value = FALLBACK_COLOR
        ) +
        # Facet: sample (row) × gene (column)
        facet_grid(
            as.formula(paste("sample_name", "~", gene_col)),
            scales = "free"
        ) +
        labs(
            x     = "Position",
            y     = "Allele Frequency",
            color = "Consequence"
        ) +
        theme_classic(base_size = 11) +
        theme(
            axis.text.x       = element_text(angle = 45, hjust = 1, vjust = 1),
            strip.text.y      = element_text(angle = 0),
            strip.background  = element_rect(fill = "grey92", color = NA),
            panel.spacing     = unit(0.8, "lines"),
            legend.position   = "bottom"
        )
}

#' Gene-level lollipop (genomic position).
plot_lollipop_gene <- function(gene_tsv = GENE_TSV, gene_widths_path = GENE_WIDTHS,
                               goi = GOI_LIST) {

    df <- read_tsv(gene_tsv, show_col_types = FALSE) |>
        filter(gene_name %in% goi)

    plot_df <- df |>
        filter(IMPACT %in% IMPACT_KEEP) |>
        filter(!str_starts(sample_name, SAMPLE_EXCLUDE))

    widths <- read_tsv(gene_widths_path, show_col_types = FALSE) |>
        filter(gene_name %in% goi)

    fake_boundary <- .build_fake_boundary(widths, "gene_name", "width", "gene_pos")

    p <- .plot_lollipop(plot_df, fake_boundary, "gene_name", "gene_pos")
    p + labs(x = "Gene Position")
}


#' CDS-level lollipop (protein position).
plot_lollipop_cds <- function(cds_tsv = CDS_TSV, goi = GOI_LIST) {

    df <- read_tsv(cds_tsv, show_col_types = FALSE) |>
        filter(gene_name %in% goi)

    plot_df <- df |>
        filter(IMPACT %in% IMPACT_KEEP) |>
        filter(!str_starts(sample_name, SAMPLE_EXCLUDE))

    # CDS: derive fake boundary from data (no genomic width equivalent)
    fake_boundary <- plot_df |>
        group_by(gene_name) |>
        summarise(protein_pos = max(protein_pos, na.rm = TRUE), .groups = "drop") |>
        tidyr::uncount(2) |>
        group_by(gene_name) |>
        mutate(
            protein_pos      = if_else(row_number() == 1L, 1L, protein_pos),
            allele_frequency = 0
        ) |>
        ungroup()

    p <- .plot_lollipop(plot_df, fake_boundary, "gene_name", "protein_pos")
    p + labs(x = "Protein Position")
}

# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Gene-level lollipop ...\n")
p_gene <- plot_lollipop_gene()
ggsave(GENE_PNG, p_gene, width = PNG_WIDTH_IN, height = PNG_HEIGHT_IN, dpi = PNG_DPI)

cat("CDS-level lollipop ...\n")
p_cds <- plot_lollipop_cds()
ggsave(CDS_PNG, p_cds, width = PNG_WIDTH_IN, height = PNG_HEIGHT_IN, dpi = PNG_DPI)

cat(sprintf("Done. → %s/\n", OUTPUT_DIR))
