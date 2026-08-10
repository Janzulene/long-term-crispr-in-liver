# ============================================================
# AAV integration plots (Fig. 2K)
#
# Reads per-sample genome and AAV-vector coverage files produced
# by the AAV pipeline (target_analyse_aav.smk) and draws:
#   - per-chromosome genome coverage of AAV-genome junction reads
#     (plot_genome_col_break, y axis break)
#   - coverage across the AAV2 vector genome
#
# Usage (from the repository root):
#   Rscript analysis/aav/plot_aav.R
#
# Outputs SVG files into OUTPUT_FIGS_PATH/{sample}/.
# ============================================================

library(conflicted)
conflicts_prefer(dplyr::filter)

library(tidyverse)

# Plotting helpers (same directory): plot_genome_col(), plot_genome_col_break(),
# scale_y_log10_plus1(), and the global chrom_sizes / CHROM_ORDER.
source("analysis/aav/plot_genome_col.R")

# ============================================================
# Configuration
# ============================================================

# All paths are relative to the repository root.
OUTPUT_FIGS_PATH <- "results/aav"
INPUT_DIR <- "data/processed/targetsequence_20240501_aav/mapping"

# Shared y-axis-break position for the genome coverage plots.
Y_BREAK_VALS <- c(9000, 45000)
# Relative heights of the top and bottom panels (bottom 1/4).
CUT_RATIO_VALS <- c(3, 1)

# Number of parallel workers.
MC_CORES <- 8

# ============================================================
# Plotting
# ============================================================

dir.create(OUTPUT_FIGS_PATH, recursive = TRUE, showWarnings = FALSE)

sample_names <- c("con_1", "con_2", "con_3", "exp_1", "exp_2", "exp_3", "tail", "pbs-1", "pbs-2")
sgRNA_list <- c("con", "angptl3_g4")

all_condition <- expand_grid(
    sample_name = sample_names,
    sgRNA       = sgRNA_list
)

parallel::mclapply(
    1:nrow(all_condition),
    \(i) {
        sample_name <- all_condition$sample_name[i]
        sgRNA_name  <- all_condition$sgRNA[i]

        genome_cov_path <- str_glue("{INPUT_DIR}/{sample_name}/{sgRNA_name}.genome.filtered.coverage.txt")
        aav_cov_path    <- str_glue("{INPUT_DIR}/{sample_name}/{sgRNA_name}.coverage.txt")

        .output_dir <- str_glue("{OUTPUT_FIGS_PATH}/{sample_name}")
        if (!dir.exists(.output_dir)) {
            dir.create(.output_dir, recursive = TRUE)
        }

        # Genome-wide distribution of AAV-genome junction reads
        cov_data <- read_tsv(genome_cov_path, col_names = c("rname", "pos", "coverage"))
        if (nrow(cov_data) == 0) {
            cov_data <- (
                chrom_sizes
                |> filter(rname %in% CHROM_ORDER)
                |> mutate(
                    rname    = factor(rname, levels = CHROM_ORDER),
                    pos      = 1,
                    coverage = 0
                )
                |> select(rname, pos, coverage)
            )
        }

        p1 <- plot_genome_col_break(
            cov_data,
            cov_key   = "coverage",
            bin_size  = 3e6,
            y_break   = Y_BREAK_VALS,
            y_lab     = "Reads Count",
            cut_ratio = CUT_RATIO_VALS
        )

        out_path1 <- str_glue("{.output_dir}/{sgRNA_name}.aav.insertion.svg")
        svg(out_path1, width = 15, height = 4)
        print(p1)
        dev.off()

        # Coverage across the AAV2 vector genome
        aav_cov_data <- read_tsv(aav_cov_path, col_names = c("rname", "pos", "coverage"))
        if (nrow(aav_cov_data) == 0) {
            aav_cov_data <- tibble(
                rname = "AAV2",
                pos = 1:4849,
                coverage = 0
            )
        }

        p2 <- (
            ggplot(aav_cov_data)
            + aes(x = pos, y = coverage)
            + geom_col(fill = "steelblue", width = 1)
            + scale_y_continuous(expand = c(0, 0))
            + theme_minimal(base_size = 18)
            + labs(y = "Reads Count", x = "AAV2 Genome Position")
        )

        out_path2 <- str_glue("{.output_dir}/{sgRNA_name}.aav.coverage.svg")
        svg(out_path2, width = 8, height = 4)
        print(p2)
        dev.off()
    },
    mc.cores = MC_CORES
)
