# ============================================================
# Genome coverage plotting library (used by plot_aav.R, Fig. 2K)
#
# Provides:
#   - scale_y_log10_plus1()      log10(x + 1) y transform for ggplot2
#   - plot_genome_col()          per-chromosome coverage bar plot
#   - plot_genome_col_break()    same, with a y-axis break (two panels)
#
# Depends on: grid, gtable, tidyverse, patchwork, ggplotify
# ============================================================

library(grid)
library(gtable)
library(tidyverse)
library(patchwork)

# Custom log10(x + 1) transformation for ggplot2
transform_log10_plus1 <- scales::new_transform(
    name = "log10_plus1",
    transform = function(x) log10(x + 1),
    inverse = function(x) 10^x - 1,
    domain = c(0, Inf)
)

scale_y_log10_plus1 <- function(...) {
    ggplot2::scale_y_continuous(..., transform = transform_log10_plus1)
}

CHROM_ORDER <- c(1:19, "X", "Y")

# Path to the GRCm39 chromosome sizes file (external reference data,
# not shipped with this repository; see data/README.md).
CHROM_SIZES_PATH <- "path/to/GRCm39.chrom.sizes"

chrom_sizes <- read_table(
    CHROM_SIZES_PATH,
    col_names = c("rname", "size")
)

# Per-chromosome coverage bar plot (manual binning).
# Args:
#   cov_data:    coverage data frame with columns rname / pos / <cov_key>
#   cov_key:     name of the coverage column
#   bin_size:    bin width in bp
#   scale_log_y: apply log10(x + 1) transform to the y axis
plot_genome_col <- function(cov_data, cov_key = "cov_rate", bin_size = 3e6, scale_log_y = TRUE) {
    binned_data <- (
        cov_data
        |> filter(rname %in% CHROM_ORDER)
        |> group_by(rname)
        |> mutate(bin = floor(pos / bin_size) * bin_size)
        |> group_by(rname, bin)
        |> summarise(
            pos     = mean(pos),
            value   = mean(!!sym(cov_key), na.rm = TRUE),
            .groups = "drop"
        )
        |> mutate(rname = factor(rname, levels = CHROM_ORDER))
    )

    chrom_placeholder <- (
        chrom_sizes
        |> filter(rname %in% CHROM_ORDER)
        |> mutate(
            rname = factor(rname, levels = CHROM_ORDER),
            xmin = 0
        )
        |> pivot_longer(-rname, names_to = "names", values_to = "pos")
        |> mutate(value = 0)
        |> select(rname, pos, value)
    )

    p <- (
        ggplot(binned_data)
        + aes(x = pos, y = value)
        + geom_col(data = chrom_placeholder, aes(x = pos, y = value), alpha = 0)
        + geom_col(fill = "black", width = bin_size * 0.9)
        + facet_grid(. ~ rname, scales = "free_x", space = "free_x", switch = "x")
        + labs(y = "Coverage Rate", x = "Chromosome")
        + theme_classic(base_size = 18)
        + theme(
            strip.text.y.right = element_text(angle = 0),
            strip.background  = element_blank(),
            strip.text.x = element_text(
                color    = "white",
                size     = 16,
                face     = "bold"
            ),
            strip.placement    = "outside",
            axis.line          = element_line(),
            axis.text.x        = element_blank()
        )
    )
    if (scale_log_y) {
        p <- p + scale_y_log10_plus1(expand = c(0, 0))
    } else {
        p <- p + scale_y_continuous(expand = c(0, 0))
    }

    # Add rounded-rectangle strip backgrounds
    g <- ggplotGrob(p)
    g <- apply_rounded_strip(g)
    p_final <- ggplotify::as.ggplot(g)
    p_final
}

# Helper: apply rounded-rectangle strip backgrounds to a ggplotGrob
apply_rounded_strip <- function(g) {
    strip_indices <- which(grepl("strip", g$layout$name))

    for (i in strip_indices) {
        strip_grob <- g$grobs[[i]]
        rounded_rect <- roundrectGrob(
            x      = unit(0.5, "npc"),
            y      = unit(0.5, "npc"),
            width  = unit(1, "npc"),
            height = unit(1, "npc"),
            r      = unit(0.15, "inches"),
            gp     = gpar(fill = "#003366", col = NA)
        )
        g$grobs[[i]] <- gTree(
            children = gList(rounded_rect, strip_grob),
            vp = strip_grob$vp
        )
    }
    return(g)
}

# Helper: generate a slash grob (used for the y-axis break)
create_slash_grob <- function(y_pos_npc) {
    segmentsGrob(
        x0 = unit(0, "npc") - unit(2, "mm"),
        x1 = unit(1, "npc") + unit(2, "mm"),
        y0 = unit(y_pos_npc, "npc") - unit(2, "mm"),
        y1 = unit(y_pos_npc, "npc") + unit(2, "mm"),
        gp = gpar(lwd = 1.5, col = "black")
    )
}

# Per-chromosome coverage bar plot with a y-axis break (two panels).
# Args:
#   cov_data: coverage data frame with columns rname / pos / <cov_key>
#   cov_key:  name of the coverage column
#   bin_size: bin width in bp
#   y_break:  numeric vector of length 2 giving the omitted range, e.g. c(9000, 45000)
#   y_lab:    y axis title
#   cut_ratio: relative heights of the top and bottom panels, e.g. c(3, 1)
plot_genome_col_break <- function(cov_data, cov_key = "cov_rate", bin_size = 3e6, y_break = NULL, y_lab = "Reads Count", cut_ratio = c(3, 1)) {
    # --- 1. Data preparation ---
    binned_data <- (
        cov_data
        |> filter(rname %in% CHROM_ORDER)
        |> group_by(rname)
        |> mutate(bin = floor(pos / bin_size) * bin_size)
        |> group_by(rname, bin)
        |> summarise(
            pos     = mean(pos),
            value   = mean(!!sym(cov_key), na.rm = TRUE),
            .groups = "drop"
        )
        |> mutate(rname = factor(rname, levels = CHROM_ORDER))
    )

    chrom_placeholder <- (
        chrom_sizes
        |> filter(rname %in% CHROM_ORDER)
        |> mutate(
            rname = factor(rname, levels = CHROM_ORDER),
            xmin = 0
        )
        |> pivot_longer(-rname, names_to = "names", values_to = "pos")
        |> mutate(value = 0)
        |> select(rname, pos, value)
    )

    # --- 2. Base plot ---
    build_base_plot <- function(data_subset) {
        p <- (
            ggplot(data_subset)
            + aes(x = pos, y = value)
            + geom_col(data = chrom_placeholder, aes(x = pos, y = value), alpha = 0)
            + geom_col(fill = "black", width = bin_size * 0.9)
            # switch = "x" places the strip at the bottom
            + facet_grid(. ~ rname, scales = "free_x", space = "free_x", switch = "x")
            + theme_classic(base_size = 18)
            + theme(
                strip.background = element_blank(),
                strip.placement  = "outside",
                axis.line        = element_line(),
                axis.text.x      = element_blank(), # chromosome labels come from the strip
            )
        )
        return(p)
    }

    # --- 3. Axis-break layout ---
    y_break   <- sort(y_break)
    gap_start <- y_break[1]
    gap_end   <- y_break[2]
    max_val   <- max(binned_data$value, na.rm = TRUE) * 1.15 # headroom at the top

    p_top <- (
        build_base_plot(
            data_subset = binned_data |> filter(value >= gap_end)
        )
        + labs(x = NULL, y = NULL)
        + theme(
            strip.text   = element_blank(),
            axis.line.x  = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin  = margin(b = 0, t = 25),
        )
        + coord_cartesian(ylim = c(gap_end, max_val), expand = FALSE)
    )

    p_bottom <- (
        build_base_plot(
            data_subset = binned_data
        )
        + labs(y = NULL, x = "Chromosome")
        + theme(
            strip.text.x = element_text(
                color    = "white",
                size     = 16,
                face     = "bold"
            ),
            plot.margin = margin(t = 4),
        )
        + coord_cartesian(ylim = c(0, gap_start), expand = FALSE)
    )

    g_top <- ggplotGrob(p_top)
    g_bottom <- apply_rounded_strip(ggplotGrob(p_bottom))

    # Align the widths of the two panels
    max_widths <- unit.pmax(g_top$widths, g_bottom$widths)
    g_top$widths <- max_widths
    g_bottom$widths <- max_widths

    g_final <- rbind(g_top, g_bottom, size = "last")

    # Set the relative heights of the two panels
    panel_heights <- g_final$layout$t[grepl("panel-1-1", g_final$layout$name)]
    g_final$heights[panel_heights[1]] <- unit(cut_ratio[1], "null")
    g_final$heights[panel_heights[2]] <- unit(cut_ratio[2], "null")

    # --- 4. Add slash marks on the y axis ---
    axis_layout <- g_final$layout[grepl("axis-l", g_final$layout$name), ]
    axis_indices <- sort(unique(axis_layout$t))

    top_axis_row <- axis_indices[1]
    bottom_axis_row <- axis_indices[2]

    axis_col_idx <- axis_layout$l[1]

    # Slash parameters
    d <- unit(3, "mm") # half width
    h <- unit(2, "mm") # half height

    slash_top <- segmentsGrob(
        x0 = unit(1, "npc") - d, x1 = unit(1, "npc") + d,
        y0 = unit(0, "npc") - h, y1 = unit(0, "npc") + h,
        gp = gpar(lwd = 2.327953, col = "black")
    )
    slash_bottom <- segmentsGrob(
        x0 = unit(1, "npc") - d, x1 = unit(1, "npc") + d,
        y0 = unit(1, "npc") - h, y1 = unit(1, "npc") + h,
        gp = gpar(lwd = 2.327953, col = "black")
    )

    g_final2 <- gtable_add_grob(
        g_final,
        slash_top,
        t    = top_axis_row,
        l    = axis_col_idx,
        z    = Inf,
        clip = "off",
        name = "axis-slash-top"
    )
    g_final3 <- gtable_add_grob(
        g_final2,
        slash_bottom,
        t    = bottom_axis_row,
        l    = axis_col_idx,
        z    = Inf,
        clip = "off",
        name = "axis-slash-bottom"
    )

    # --- 5. Add the global y axis title ---
    y_lab_grob <- textGrob(
        y_lab,
        rot = 90,
        gp  = gpar(fontsize = 18, col = "black")
    )
    lab_width <- grobWidth(y_lab_grob) + unit(4, "mm")

    g_final4 <- gtable_add_cols(
        g_final3,
        width = lab_width,
        pos = 0
    )
    g_final4 <- gtable_add_grob(
        g_final4,
        y_lab_grob,
        t    = 1,
        l    = 1,
        b    = nrow(g_final4),
        r    = 1,
        z    = Inf,
        clip = "off",
        name = "y-axis-label-global"
    )

    p_final <- ggplotify::as.ggplot(g_final4)
    p_final
}
