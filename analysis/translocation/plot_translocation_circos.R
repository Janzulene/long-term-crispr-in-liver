# ============================================================
# Translocation circos plots (Fig. 2L)
#
# Reads per-sample translocation calls and on-target sequencing
# depth from the target sequencing pipeline output, and draws
# circos plots (circlize) showing translocation events from the
# targeted loci to other genomic regions, with the outer track
# representing on-target read depth.
#
# Usage:
#   Rscript analysis/translocation/plot_translocation_circos.R
#
# Outputs SVG files into OUTPUT_FIGS_PATH:
#   - trans_circle_{sample}_{region}.svg         (per sample x region)
#   - sum_trans_circle_{sample}.svg              (regions pooled per sample)
#   - sum_trans_circle_{region}.svg              (samples pooled per region)
# ============================================================

library(conflicted)
conflicts_prefer(dplyr::filter)

library(circlize)
library(tidyverse)

# ============================================================
# Configuration
# ============================================================

# Path to the GRCm39 chromosome sizes file (external reference
# data, not shipped with this repository; see data/README.md).
CHROM_SIZES_PATH <- "path/to/GRCm39.chrom.sizes"

# All other paths are relative to the repository root.
TARGET_DATA_PATH <- "data/processed/targetsequence_20240501"
OUTPUT_FIGS_PATH <- "results/translocation"

# Only report events supported by more than this many read pairs.
MIN_N_TRANS <- 2
# Number of parallel workers.
MC_CORES <- 8

# ============================================================
# Load reference data
# ============================================================

VALID_CHROMOSOMES <- c(c(1:19), "X", "Y")

chrom_sizes <- read_table(CHROM_SIZES_PATH, col_names = c("chrom", "size")) |>
    filter(str_to_upper(chrom) %in% str_to_upper(VALID_CHROMOSOMES)) |>
    mutate(chrom = factor(chrom, levels = VALID_CHROMOSOMES))

plot_chrom <- chrom_sizes |> mutate(start = 0) |> arrange(chrom)

# ============================================================
# Helper functions
# ============================================================

# Read on-target depth data for one sample and one target region.
.read_depth <- function(sample_name, primer_prefix) {
    depth_path <- glue::glue("{TARGET_DATA_PATH}/{sample_name}/.dirty/{primer_prefix}/intarget__cover.depth")
    read_tsv(
        depth_path,
        col_names = c("chr", "pos", "depth"),
        col_types = cols(
            chr   = col_character(),
            pos   = col_integer(),
            depth = col_integer()
        )
    )
}

# Read translocation calls for one sample and one target region.
# The forward and reverse primers of the region are pooled, and
# reads are aggregated into unique translocation events.
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
        sum_trans_res <- tibble(
            n                        = integer(),
            intarget_chr             = character(),
            intarget_breakpoint_pos  = integer(),
            offtarget_chr            = character(),
            offtarget_breakpoint_pos = integer()
        )
    }

    sum_trans_res
}

# Draw one circos plot.
# Args:
#   depth_df:      per-position depth (columns: chr / pos / depth)
#   sum_trans_res: translocation events (columns: n / intarget_chr /
#                  intarget_breakpoint_pos / offtarget_chr / offtarget_breakpoint_pos)
#   col_fun:       color scale for event counts
#   title_name:    plot title
#   max_value:     count mapped to the darkest color
# NOTE: relies on the global `plot_chrom`.
plot_trans <- function(
    depth_df,
    sum_trans_res,
    col_fun    = colorRamp2(c(0, log(max_value + 1)), c("white", "red")),
    title_name = "trans plot",
    max_value  = 21590
) {

    # Initialize the circos layout
    circos.clear()
    circos.par(
        start.degree = 90,
        track.margin = c(0, 0.02)
    )

    xlim_mat <- plot_chrom |> select(start, size) |> as.matrix()
    circos.initialize(factors = plot_chrom$chrom, xlim = xlim_mat)

    # Track 1: sequencing depth histogram (outermost ring)
    .max_depth <- max(depth_df$depth, na.rm = TRUE) + 1
    circos.track(
        factors      = plot_chrom$chrom,
        track.height = 0.25,
        ylim         = c(0, .max_depth |> log()),
        bg.border    = NA,
        panel.fun    = function(x, y) {
            current_chr <- CELL_META$sector.index
            chr_data <- depth_df %>% filter(chr == current_chr)

            # Bin the depth profile (5 Mbp bins)
            bin_size <- 5e6
            current_xlim <- CELL_META$xlim
            breaks <- seq(0, current_xlim[2], by = bin_size)
            if (max(breaks) < current_xlim[2]) breaks <- c(breaks, current_xlim[2])

            depth_means <- tapply(chr_data$depth, cut(chr_data$pos, breaks = breaks, include.lowest = TRUE), mean, na.rm = TRUE)

            for (i in seq_along(depth_means)) {
                if (!is.na(depth_means[i])) {
                    circos.rect(
                        xleft   = breaks[i],
                        xright  = breaks[i + 1],
                        ybottom = 0,
                        ytop    = log(depth_means[i] + 1),
                        col     = "#2885f7",
                        border  = NA
                    )
                }
            }
        }
    )

    # Track 2: chromosome labels (middle ring)
    circos.track(
        factors      = plot_chrom$chrom,
        track.height = 0.08,
        ylim         = c(0, 1),
        panel.fun = function(x, y) {
            circos.rect(
                CELL_META$cell.xlim[1],
                CELL_META$cell.ylim[1],
                CELL_META$cell.xlim[2],
                CELL_META$cell.ylim[2],
                col = "#DDDDDD",
                border = NA
            )
            circos.text(
                CELL_META$xcenter,
                CELL_META$ycenter,
                CELL_META$sector.index,
                col = "white",
                cex = 0.9,
                facing = "clockwise",
                niceFacing = TRUE
            )
        },
        bg.border = NA
    )

    # Track 3: translocation links (drawn lightest first)
    if (nrow(sum_trans_res) > 0) {
        sum_trans_res %>%
            arrange(n) %>%
            mutate(log_value = log(n + 1)) %>%
            pmap(function(...) {
                data <- list(...)
                circos.link(
                    sector.index1 = data$intarget_chr,
                    point1        = data$intarget_breakpoint_pos,
                    sector.index2 = data$offtarget_chr,
                    point2        = data$offtarget_breakpoint_pos,
                    col           = col_fun(data$log_value),
                    border        = NA,
                    lwd           = 2
                )
            })
    }

    title(title_name)
}

# ============================================================
# Plotting
# ============================================================

dir.create(OUTPUT_FIGS_PATH, recursive = TRUE, showWarnings = FALSE)

sample_names <- c("con_1", "con_2", "con_3", "exp_1", "exp_2", "exp_3", "tail", "pbs-1", "pbs-2")
region_list  <- c("angptl3-g4", "pcsk9-g1", "pcsk9-g3")

# Mode 1: one plot per sample x region
parallel::mclapply(
    cross2(sample_names, region_list),
    \(x) {
        sample_name <- x[[1]]
        region <- x[[2]]

        depth_df <- .read_depth(sample_name, region)
        trans_df <- .read_trans(sample_name, region) |> filter(n > MIN_N_TRANS)

        output_path <- glue::glue("{OUTPUT_FIGS_PATH}/trans_circle_{x[[1]]}_{x[[2]]}.svg")
        svg(output_path)
        if (nrow(depth_df) == 0 || nrow(trans_df) == 0) {
            # Draw an empty placeholder
            plot.new()
            plot.window(xlim = c(0, 1), ylim = c(0, 1))

            reason <- if (nrow(depth_df) == 0) "No depth data" else "No trans data"
            text(0.5, 0.5,
                paste0(sample_name, ":", region, "\n", reason),
                cex = 1.5, col = "gray50"
            )

            message(glue::glue("{reason} for `trans_circle_{x[[1]]}_{x[[2]]}`, plot placeholder"))
        } else {
            plot_trans(depth_df, trans_df, title_name = glue::glue("{sample_name}_{region}"))
        }
        dev.off()
        circos.clear()

        message(glue::glue("plot circlize plot for {sample_name} in {region}"))
        output_path
    },
    mc.cores = MC_CORES
)

# Mode 2: one plot per sample, regions pooled
parallel::mclapply(
    sample_names,
    \(sample_name) {
        depth_df <- (
            lapply(region_list, \(region) .read_depth(sample_name, region))
            |> bind_rows()
            |> group_by(chr, pos)
            |> summarise(depth = sum(depth), .groups = "drop")
        )

        trans_df <- (
            lapply(region_list, \(region) .read_trans(sample_name, region))
            |> bind_rows()
            |> group_by(across(-c(n)))
            |> summarise(n = sum(n), .groups = "drop")
            |> select(n, everything())
            |> filter(n > MIN_N_TRANS)
        )

        output_path <- glue::glue("{OUTPUT_FIGS_PATH}/sum_trans_circle_{sample_name}.svg")
        svg(output_path)
        if (nrow(depth_df) == 0 || nrow(trans_df) == 0) {
            plot.new()
            plot.window(xlim = c(0, 1), ylim = c(0, 1))

            reason <- if (nrow(depth_df) == 0) "No depth data" else "No trans data"
            text(0.5, 0.5,
                paste0(sample_name, "\n", reason),
                cex = 1.5, col = "gray50"
            )

            message(glue::glue("{reason} for `sum of {sample_name}`, plot placeholder"))
        } else {
            plot_trans(depth_df, trans_df, title_name = glue::glue("{sample_name}_all"))
        }
        dev.off()
        circos.clear()

        message(glue::glue("plot sum circlize plot for {sample_name}"))
        output_path
    },
    mc.cores = MC_CORES
)

# Mode 3: one plot per region, samples pooled
parallel::mclapply(
    region_list,
    \(region) {
        depth_df <- (
            lapply(sample_names, \(sample_name) .read_depth(sample_name, region))
            |> bind_rows()
            |> group_by(chr, pos)
            |> summarise(depth = sum(depth), .groups = "drop")
        )

        trans_df <- (
            lapply(sample_names, \(sample_name) .read_trans(sample_name, region))
            |> bind_rows()
            |> group_by(across(-c(n)))
            |> summarise(n = sum(n), .groups = "drop")
            |> select(n, everything())
            |> filter(n > MIN_N_TRANS)
        )

        output_path <- glue::glue("{OUTPUT_FIGS_PATH}/sum_trans_circle_{region}.svg")
        svg(output_path)

        if (nrow(depth_df) == 0 || nrow(trans_df) == 0) {
            plot.new()
            plot.window(xlim = c(0, 1), ylim = c(0, 1))

            reason <- if (nrow(depth_df) == 0) "No depth data" else "No trans data"
            text(0.5, 0.5,
                paste0(region, "\n", reason),
                cex = 1.5, col = "gray50"
            )

            message(glue::glue("{reason} for {region}, plot placeholder"))
        } else {
            plot_trans(depth_df, trans_df, title_name = glue::glue("{region}_all"))
            message(glue::glue("plot sum circlize plot for {region}"))
        }
        dev.off()
        circos.clear()

        output_path
    },
    mc.cores = MC_CORES
)
