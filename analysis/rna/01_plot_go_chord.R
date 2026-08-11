# ══════════════════════════════════════════════════════════════════════════════
# GO chord diagram of DAVID enrichment for key DEGs (Fig. 2G)
#
# Reads the limma differential-expression results and the DAVID GO
# enrichment table, and draws a GO chord diagram (GOplot) connecting
# key DEGs (left) to the enriched pathways (right), with the ribbon
# color indicating pathway membership and the bar color the mean
# log2 fold change of each gene.
#
# Article references:
#   - Fig. 2G: DAVID enrichment analysis of significant pathways
#     associated with key DEGs
#
# Input (data/raw/rna/, shipped with this repository):
#   - Angptl3G4_vs_NT_Limma_results.csv, Pcsk9G1_vs_NT_Limma_results.csv
#     (limma DEG results of the two edited groups vs NT control)
#   - GO_terms_from_DAVID.csv (DAVID functional annotation of key DEGs)
#
# Usage:
#   Rscript analysis/rna/01_plot_go_chord.R
#
# Output:
#   reports/figures/rna/GO_chord.svg
# ══════════════════════════════════════════════════════════════════════════════

library(GOplot)
library(dplyr)
library(tidyr)
library(stringr)
library(export)
library(conflicted)

# Specify function preferences to avoid namespace conflicts
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")

# ------------------------------------------------------------------------------
# 1. Data import
# ------------------------------------------------------------------------------
# Differential expression results for the two edited groups
DEGP <- read.csv(file = "data/raw/rna/Pcsk9G1_vs_NT_Limma_results.csv", row.names = 1)
DEGP$gene <- row.names(DEGP)

DEGA <- read.csv(file = "data/raw/rna/Angptl3G4_vs_NT_Limma_results.csv", row.names = 1)
DEGA$gene <- row.names(DEGA)

# DAVID GO enrichment results
go_data <- read.csv(file = "data/raw/rna/GO_terms_from_DAVID.csv")

# ------------------------------------------------------------------------------
# 2. Data processing and integration
# ------------------------------------------------------------------------------
# Standardize gene symbols to uppercase for consistent matching
DEGP <- DEGP %>% mutate(gene_upper = toupper(gene))
DEGA <- DEGA %>% mutate(gene_upper = toupper(gene))

# Merge datasets and calculate mean log-fold change
de_combined <- DEGP %>%
  full_join(DEGA, by = "gene_upper", suffix = c("_group1", "_group2")) %>%
  mutate(avg_logFC = (as.numeric(logFC_group1) + as.numeric(logFC_group2)) / 2) %>%
  select(gene_upper, avg_logFC)

# Extract and parse gene lists from GO results
all_genes <- go_data$Genes
genes_vector <- unlist(strsplit(toupper(all_genes), ",\\s*"))
unique_genes <- unique(genes_vector)

# Filter for genes present in the GO analysis
diffgene <- de_combined[which(de_combined$gene_upper %in% unique_genes), ]
diffgene$avg_logFC <- as.numeric(diffgene$avg_logFC)

# ------------------------------------------------------------------------------
# 3. Visualization data preparation
# ------------------------------------------------------------------------------
# Format data for GOplot compatibility (columns: ID, logFC)
diffgene$ID <- diffgene$gene_upper
diffgene$logFC <- diffgene$avg_logFC
diffgene <- diffgene[, c("ID", "logFC")]

# Construct the GOplot input format from the raw GO data
# Select necessary columns and map P.Value to adj_pval
goenrichment <- go_data %>%
  select(Category, ID, Term, Genes, adj_pval = P.Value)

# Generate the circle data matrix
circ <- circle_dat(goenrichment, diffgene)

# Convert gene symbols to title case for publication-quality display
circ$genes <- str_to_title(circ$genes)
diffgene$ID <- str_to_title(diffgene$ID)

# Generate chord data matrix
chord <- chord_dat(circ, diffgene, goenrichment$Term)

# ------------------------------------------------------------------------------
# 4. Plotting
# ------------------------------------------------------------------------------
# Define color palette for the chord diagram
chordcol <- c("#529b61", "#ffec85", "#A14462", "#c2e4ff")

# Generate the GO Chord plot
p <- GOChord(chord,
             space = 0.01,
             gene.order = "logFC",
             ribbon.col = chordcol,
             lfc.col = c("#E53A40", "white", "#30A9DE"),
             lfc.min = -4,
             lfc.max = 3,
             gene.space = 0.2,
             gene.size = 2,
             process.label = 2)

print(p)

# Save the figure
dir.create("reports/figures/rna", showWarnings = FALSE, recursive = TRUE)
svg("reports/figures/rna/GO_chord.svg", width = 10, height = 8)
print(p)
dev.off()
