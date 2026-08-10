# QC plots for the AAV pipeline (snakemake script rule):
#   - full genome coverage (all reads)
#   - genome coverage of AAV-mapped reads (insertion sites)
#   - coverage across the AAV2 vector
library(tidyverse)

source("analysis/aav/plot_genome_col.R")

full_genome_cov_path <- snakemake@input[["full_genome_cov"]]
genome_cov_path      <- snakemake@input[["genome_cov"]]
aav_cov_path         <- snakemake@input[["aav_cov"]]

# Full genome coverage
full_cov_data <- read_tsv(full_genome_cov_path, col_names = c("rname", "pos", "coverage"))
p0 <- plot_genome_col(full_cov_data, cov_key = "coverage", bin_size = 3e6) + labs(y = "Reads Count", x = "Chromosome")
svg(snakemake@output[["full_genome_plot"]], width = 15, height = 4)
print(p0)
dev.off()

# Genome coverage of AAV-mapped reads
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

p1 <- plot_genome_col(cov_data, cov_key = "coverage", bin_size = 3e6) + labs(y = "Reads Count", x = "Chromosome")
svg(snakemake@output[["genome_plot"]], width = 15, height = 4)
print(p1)
dev.off()

# AAV2 vector coverage
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

svg(snakemake@output[["aav_plot"]], width = 8, height = 4)
print(p2)
dev.off()
