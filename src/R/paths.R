# Shared path conventions for all R analyses.
#
# All paths are relative to the repository root: run every script with
# `Rscript analysis/...` from the repository root. The Python side keeps
# an equivalent configuration in src/config.py.
#
# Directory names deliberately carry no batch date: the published
# pipeline targets a single final analysis round.

# --- data roots ---
DATA_RAW       <- "data/raw"
DATA_PROCESSED <- "data/processed"
DATA_FINAL     <- "data/final"
FIGURES        <- "reports/figures"

# --- targeted-sequencing module (translocation + AAV) ---
TS_PROCESSED     <- file.path(DATA_PROCESSED, "targetsequence")
TS_PROCESSED_AAV <- file.path(DATA_PROCESSED, "targetsequence_aav")
TS_FINAL         <- file.path(DATA_FINAL, "targetsequence")
TS_FIGURES       <- file.path(FIGURES, "translocation")
AAV_FIGURES      <- file.path(FIGURES, "aav")

# --- WES module ---
WES_PROCESSED <- file.path(DATA_PROCESSED, "WES_analysis")
WES_FINAL     <- file.path(DATA_FINAL, "WES")
WES_FIGURES   <- file.path(FIGURES, "WES")

# --- external inputs (not shipped with the repository; see data/README.md) ---
# Mutect2 + VEP annotated variant tables (sarek 3.4.2 output)
SAREK_VCF_TABLE   <- file.path(DATA_PROCESSED, "WES2_paired_annotate", "annotation", "mutect2")
# per-sample mean depth tables (samtools depth, sarek preprocessing output)
SAREK_DEPTH_TABLE <- file.path(DATA_PROCESSED, "WES2_paired", "preprocessing", "markduplicates")
# cancer gene list
CANCER_GENE_LIST  <- file.path(DATA_RAW, "WES", "cancer_list.csv")
# AAV vector references (one fasta per sgRNA)
AAV_REFERENCE_DIR <- file.path(DATA_RAW, "aav")
