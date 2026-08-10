# Functions for early-stage target-sequencing analysis (per-region
# translocation rate estimation from BAM files).
#
# NOTE: the final paper pipeline uses the Python translocation detector
# (src/detect_translocation.py) driven by pipeline/target_analyse_2.smk.
# These R helpers were used during method development and for the
# per-region summary tables; they are kept for reference and for
# reproducing the summary statistics.

.read_bam <- function(bam_file) {
    bam_file <- Rsamtools::BamFile(bam_file)

    bam_data <- Rsamtools::scanBam(bam_file)[[1]]
    bam_df <- as.data.frame(bam_data) |> dplyr::as_tibble()
    bam_df
}

.analyse_region <- function(
    bam_df, region_name, target_chr = 4,
    region_start = NULL, region_end = NULL, only_primary = TRUE, return_trans_reads = FALSE) {
    # if region_start is null, will use the whole chromosome to analyse
    total_mapped <- .count_reads(bam_df)

    # get region_df
    region_df <- bam_df |>
        filter(rname == target_chr)
    if (!is.null(region_start)) {
        region_df <- region_df |>
            filter(region_start <= pos, pos <= region_end)
    }
    region_df <- region_df |>
        mutate(n_match = .count_match_cigar(cigar))

    if (only_primary) {
        # only consider about primary reads
        region_df <- region_df |> filter(bitwAnd(flag, 2048 + 256) == 0)
    }
    n_mapped <- .count_reads(region_df)
    trans_reads <- region_df |>
        filter(mrnm != target_chr) |>
        filter(30 < n_match, n_match < 150)
    n_trans <- .count_reads(trans_reads)

    region_report <- tidyr::tibble(
        region = region_name,
        n_mapped = n_mapped,
        total_mapeed = total_mapped,
        mapping_rate = n_mapped / total_mapped,
        n_trans = n_trans,
        trans_rate = n_trans / n_mapped
    )

    if (return_trans_reads) {
        trans_reads_df <- bam_df |> filter(qname %in% trans_reads$qname)
        return(list(region_report, trans_reads_df))
    }
    region_report
}

# Revised region analysis
.analyse_region_2 <- function(
    bam_df,
    target_chr,
    left_region_start,
    left_region_end,
    right_region_start,
    right_region_end,
    region_name,
    only_primary = FALSE,
    return_trans_reads = FALSE,
    return_region_df = FALSE) {
    # Fetch the region df early to reduce downstream work
    region_df <- bam_df |>
        filter(rname == target_chr) |>
        filter(
            pos >= left_region_start,
            pos <= right_region_end
        ) |>
        mutate(right_pos = pos + qwidth)

    # Region dfs for the two primer directions
    left_region_df <- region_df |>
        filter(
            pos >= left_region_start,
            pos <= left_region_end
        )
    right_region_df <- region_df |>
        filter(
            right_pos >= right_region_start + 1,
            right_pos <= right_region_end + 1
        )

    # Combine the region dfs from both sides
    region_df <- bind_rows(left_region_df, right_region_df) |>
        mutate(n_match = .count_match_cigar(cigar)) |>
        filter(n_match >= 40)

    if (only_primary) {
        # only consider about primary reads
        region_df <- region_df |> filter(bitwAnd(flag, 2048 + 256) == 0)
    }
    n_mapped <- .count_reads(region_df)
    n_left_primer <- region_df |>
        filter((pos == left_region_start) | (right_pos == right_region_start + 1)) |>
        distinct(qname) |>
        nrow()
    n_right_primer <- region_df |>
        filter((pos == left_region_end) | (right_pos == right_region_end + 1)) |>
        distinct(qname) |>
        nrow()
    trans_reads <- region_df |>
        filter(mrnm != target_chr) |>
        filter(n_match < 150) # a real translocation should not match the target perfectly
    n_trans <- .count_reads(trans_reads)

    region_report <- tidyr::tibble(
        region = region_name,
        n_mapped = n_mapped,
        n_trans = n_trans,
        trans_rate = n_trans / n_mapped,
        n_left_primer = n_left_primer,
        n_right_primer = n_right_primer
    )

    if (return_trans_reads) {
        trans_reads_qnames <- unique(trans_reads$qname)
        return(list(region_report, trans_reads_qnames))
    }

    if (return_region_df) {
        return(list(region_report, region_df))
    }

    region_report
}

# Further improved analysis, less manual work
#' f_primer_pos: forward primer position
#' r_primer_pos: reverse primer end position
.analyse_region_3 <- function(
    bam_df,
    target_chr,
    f_primer_pos,
    r_primer_pos,
    region_name,
    only_primary = FALSE,
    return_other = FALSE
) {
    if (f_primer_pos > r_primer_pos) {
        # f > r: reverse strand
        strand <- "-"

        f_primer_rightpos <- f_primer_pos - 150
        r_primer_rightpos <- r_primer_pos + 150

        range_start <- min(f_primer_rightpos, r_primer_pos)
        range_end   <- max(f_primer_pos, r_primer_rightpos)

        l_primer_start <- r_primer_pos # l for low
        l_primer_end   <- r_primer_rightpos
        h_primer_start <- f_primer_rightpos # h for high
        h_primer_end   <- f_primer_pos

    } else {
        # f < r: forward strand
        strand <- "+"

        f_primer_rightpos <- f_primer_pos + 150
        r_primer_rightpos <- r_primer_pos - 150

        range_start <- min(f_primer_pos, r_primer_rightpos)
        range_end   <- max(f_primer_rightpos, r_primer_pos)

        l_primer_start <- f_primer_pos # l for low
        l_primer_end   <- f_primer_rightpos
        h_primer_start <- r_primer_rightpos # h for high
        h_primer_end   <- r_primer_pos
    }

    # Widen the range a little
    range_start <- range_start - 3
    range_end <- range_end + 3

    # Fetch the region df early to reduce downstream work
    region_df <- bam_df |>
        filter(rname == target_chr) |>
        mutate(right_pos = pos + qwidth) |>
        filter(
            pos <= range_end,
            right_pos >= range_start
        ) |>
        mutate(n_match = .count_match_cigar(cigar)) |>
        filter(n_match >= 40)

    if (only_primary) {
        # only consider about primary reads
        region_df <- region_df |> filter(bitwAnd(flag, 2048 + 256) == 0)
    }

    # Compute the metrics
    n_mapped <- .count_reads(region_df)
    n_l_primer <- region_df |>
        filter(
            (between(pos, l_primer_start - 3, l_primer_start + 3))
            | (between(right_pos, l_primer_end - 3, l_primer_end + 3))
        ) |>
        distinct(qname) |>
        nrow()
    n_h_primer <- region_df |>
        filter(
            (between(pos, h_primer_start - 3, h_primer_start + 3))
            | (between(right_pos, h_primer_end - 3, h_primer_end + 3))
        ) |>
        distinct(qname) |>
        nrow()

    if (strand == "+") {
        n_f_primer <- n_l_primer
        n_r_primer <- n_h_primer
    } else {
        n_f_primer <- n_h_primer
        n_r_primer <- n_l_primer
    }

    trans_reads <- region_df |>
        filter(mrnm != target_chr) |>
        filter(n_match < 150) # a real translocation should not match the target perfectly
    n_trans <- .count_reads(trans_reads)

    region_report <- tidyr::tibble(
        region = region_name,
        n_mapped = n_mapped,
        n_trans = n_trans,
        trans_rate = n_trans / n_mapped,
        n_forward_primer = n_f_primer,
        n_reverse_primer = n_r_primer
    )

    if (return_other) {
        trans_reads_qnames <- unique(trans_reads$qname)
        region_df <- region_df |> mutate(region = region_name)
        return(list(region_report, trans_reads_qnames, region_df))
    }

    region_report
}

# Stricter analysis conditions:
# 1. keep only reads anchored at the correct primer start/end (+/- 3 bp)
# 2. the other segment must have mapping quality > 10
# 3. primary mappings only
# 4. mapping length should be < 120
.analyse_region_strict <- function(
    bam_df, # bam_df must be the complete validated_df
    target_chr,
    f_primer_pos,
    r_primer_pos,
    region_name,
    return_other = FALSE
) {
    if (f_primer_pos > r_primer_pos) {
        # f > r: reverse strand
        strand <- "-"

        f_primer_rightpos <- f_primer_pos - 150
        r_primer_rightpos <- r_primer_pos + 150

        range_start <- min(f_primer_rightpos, r_primer_pos)
        range_end   <- max(f_primer_pos, r_primer_rightpos)

        l_primer_start <- r_primer_pos # l for low
        l_primer_end   <- r_primer_rightpos
        h_primer_start <- f_primer_rightpos # h for high
        h_primer_end   <- f_primer_pos

    } else {
        # f < r: forward strand
        strand <- "+"

        f_primer_rightpos <- f_primer_pos + 150
        r_primer_rightpos <- r_primer_pos - 150

        range_start <- min(f_primer_pos, r_primer_rightpos)
        range_end   <- max(f_primer_rightpos, r_primer_pos)

        l_primer_start <- f_primer_pos # l for low
        l_primer_end   <- f_primer_rightpos
        h_primer_start <- r_primer_rightpos # h for high
        h_primer_end   <- r_primer_pos
    }
    # Widen the range a little
    range_start <- range_start - 3
    range_end <- range_end + 3

    # Fetch the region df early to reduce downstream work
    region_df <- bam_df |>
        filter(rname == target_chr) |>
        mutate(right_pos = pos + qwidth) |>
        filter(
            pos <= range_end,
            right_pos >= range_start
        ) |>
        filter(bitwAnd(flag, 2048 + 256) == 0) |> # only primary
        mutate(n_match = .count_match_cigar(cigar)) |>
        filter(n_match >= 40)

    # Keep only reads anchored at the primers
    low_primer_df <- region_df |>
        filter(
            (between(pos, l_primer_start - 3, l_primer_start + 3))
            | (between(right_pos, l_primer_end - 3, l_primer_end + 3))
        )
    high_primer_df <- region_df |>
        filter(
            (between(pos, h_primer_start - 3, h_primer_start + 3))
            | (between(right_pos, h_primer_end - 3, h_primer_end + 3))
        )

    n_l_mapped <- .count_reads(low_primer_df)
    n_h_mapped <- .count_reads(high_primer_df)

    low_trans_reads <- low_primer_df |>
        filter(mrnm != target_chr) |>
        filter(n_match < 150)
    high_trans_reads <- high_primer_df |>
        filter(mrnm != target_chr) |>
        filter(n_match < 150)

    # Analyse the other half of the reads
    l_tr_qname <- bam_df |>
        filter(rname != target_chr) |>
        filter(qname %in% low_trans_reads$qname) |>
        filter(bitwAnd(flag, 2048 + 256) == 0) |> # also should be primary
        filter(mapq > 10) |> # further mapq filtering
        distinct(qname)

    h_tr_qname <- bam_df |>
        filter(rname != target_chr) |>
        filter(qname %in% high_trans_reads$qname) |>
        filter(bitwAnd(flag, 2048 + 256) == 0) |> # also should be primary
        filter(mapq > 10) |> # further mapq filtering
        distinct(qname)

    n_l_trans <- l_tr_qname |> nrow()
    n_h_trans <- h_tr_qname |> nrow()

    if (strand == "+") {
        h_label <- "reverse"
        l_label <- "forward"
    } else {
        h_label <- "forward"
        l_label <- "reverse"
    }

    region_report <- tidyr::tibble(
        region     = region_name,
        n_mapped   = n_l_mapped + n_h_mapped,
        n_trans    = n_l_trans + n_h_trans,
        trans_rate = (n_l_trans + n_h_trans) / (n_l_mapped + n_h_mapped),
    )

    if (return_other) {

        low_trans_reads <- bam_df |>
            filter(qname %in% l_tr_qname$qname) |>
            mutate(region = region_name, primer = l_label)

        high_trans_reads <- bam_df |>
            filter(qname %in% h_tr_qname$qname) |>
            mutate(region = region_name, primer = h_label)

        trans_reads <- bind_rows(low_trans_reads, high_trans_reads)

        return(list(region_report, trans_reads))
    }

    region_report
}

# Count reads via the read1/read2 flags
.count_reads <- function(bam_df) {
    n_read1 <- bam_df %>%
        filter(bitwAnd(flag, 64) != 0) %>%
        distinct(qname) %>%
        nrow()
    n_read2 <- bam_df %>%
        filter(bitwAnd(flag, 128) != 0) %>%
        distinct(qname) %>%
        nrow()
    n_read1 + n_read2
}

# Optimized fast count of M operations in CIGAR strings
.count_match_cigar <- function(cigar) {
    matchs <- str_replace_all(cigar, r"{(\d+(?![M\d]).)}", "")

    lapply(
        str_match_all(matchs, r"{\d+(?=[M])}"),
        \(x) {
            as.numeric(x) %>%
                sum()
        }
    ) %>% unlist()
}

.get_chr_bincount <- function(bam_df, chr, n_bin = 1000, sample_size = NULL) {
    bam_df <- bam_df |>
        filter(rname == chr)

    if (!is.null(sample_size)) {
        bam_df <- bam_df |>
            dplyr::sample_n(size = sample_size)
    }
    bam_df |>
        dplyr::mutate(bin = cut(pos, n_bin)) |>
        dplyr::group_by(bin) |>
        mutate(bin_reads = dplyr::n())
}

.get_chr_slidecount <- function(
    bam_df, chr, sample_size = NULL, slide_window = 2e+05) {
    # Sliding-window read counts
    # sample_size: if NULL, use the whole chromosome
    bam_df <- bam_df |>
        dplyr::filter(rname == chr)

    if (!is.null(sample_size)) {
        bam_df <- bam_df |>
            dplyr::sample_n(size = sample_size)
    }

    bam_df |>
        dplyr::arrange(pos) |>
        dplyr::mutate(
            slide_count = slider::slide_index_dbl(
                .x = qname,
                .i = pos,
                .f = \(x) length(x),
                .before = slide_window / 2,
                .after = slide_window / 2
            )
        ) |>
        dplyr::arrange(-slide_count)
}

# Analysis for specific regions
.analyse_chr4 <- function(bam_df) {
    .analyse_region(bam_df, "chr_4", target_chr = 4)
}

.analyse_gene_pcsk9 <- function(bam_df) {
    .analyse_region(
        bam_df,
        "gene_pcsk9",
        target_chr = 4,
        region_start = 106299526,
        region_end = 106321526
    )
}

.analyse_gene_angptl3 <- function(bam_df) {
    .analyse_region(
        bam_df,
        "gene_angptl3",
        target_chr = 4,
        region_start = 98919191,
        region_end = 98934348
    )
}

.analyse_pcsk9_g1 <- function(bam_df) {
    .analyse_region(
        bam_df,
        "pcsk9_g1",
        target_chr = 4,
        region_start = 106320987,
        region_end = 106320987 + 150
    )
}

.analyse_pcsk9_g3 <- function(bam_df) {
    .analyse_region(
        bam_df,
        "pcsk9_g3",
        target_chr = 4,
        region_start = 106314006,
        region_end = 106314194
    )
}

.analyse_angptl3_g4 <- function(bam_df) {
    .analyse_region(
        bam_df,
        "angptl3_g4",
        target_chr = 4,
        region_start = 98920295,
        region_end = 98920449
    )
}
