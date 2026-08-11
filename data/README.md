# Data

Raw sequencing data for this study are deposited in the CNGB Sequence
Archive (CNSA) of the China National GeneBank DataBase (CNGBdb):

- **Accession**: `CNP0008752`
- **Contents**: raw FASTQ files for whole-exome sequencing (WES) and
  transposon-mediated targeted enrichment sequencing libraries

See the paper's Data Availability statement for details.

## Expected local layout

All paths referenced by the pipelines are relative to this repository
and follow the layout below (raw FASTQ and pipeline outputs are
git-ignored):

```
data/
├── raw/
│   ├── <fastq from CNSA>           # targeted-sequencing libraries (downloaded)
│   ├── WES/
│   │   └── cancer_list.csv         # cancer gene list
│   ├── aav/
│   │   └── aav_{sgRNA}.fa          # AAV vector reference per sgRNA
│   └── rna/
│       ├── Angptl3G4_vs_NT_Limma_results.csv  # limma DEG results
│       ├── Pcsk9G1_vs_NT_Limma_results.csv    # limma DEG results
│       └── GO_terms_from_DAVID.csv            # DAVID GO enrichment
├── processed/                      # intermediate outputs (fastp, BWA, filtering, ...)
│   ├── targetsequence/             # translocation pipeline output
│   ├── targetsequence_aav/         # AAV pipeline output
│   ├── WES_analysis/               # WES analysis intermediates (RDS)
│   ├── WES2_paired/                # sarek 3.4.2 preprocessing output (external)
│   └── WES2_paired_annotate/       # sarek Mutect2 + VEP annotated VCF tables
│                                   # (annotation/mutect2/*.tsv, external)
├── final/                          # final result tables
│   ├── targetsequence/             # translocation / AAV summary tables
│   └── WES/                        # WES result tables
└── (figures are written to reports/figures/)
```

The `WES2_paired*` directories are produced by the external
[sarek 3.4.2](https://nf-co.re/sarek) pipeline (fastp → BWA-MEM →
MarkDuplicates → Mutect2 paired tumor-normal → VEP v110) and are not
part of this repository; see the paper's methods for details.

