# Data

Raw sequencing data for this study are deposited in the CNGB Sequence
Archive (CNSA) of the China National GeneBank DataBase (CNGBdb):

- **Accession**: `CNP0008752`
- **Contents**: raw FASTQ files for whole-exome sequencing (WES) and
  transposon-mediated targeted enrichment sequencing libraries

See the paper's Data Availability statement for details.

## Expected local layout

All paths referenced by the pipelines are relative to this repository
and follow the layout below (git-ignored, not shipped):

```
data/
├── raw/                            # input data, downloaded / copied by the user
│   ├── <fastq from CNSA>           # targeted-sequencing libraries
│   ├── WES/
│   │   └── cancer_list.csv         # cancer gene list used by analysis/wes/01_load_data.R
│   └── aav/
│       └── aav_{sgRNA}.fa          # AAV vector reference per sgRNA
│                                   # (target_analyse_aav.smk)
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

## External reference files (not shipped)

- GRCm39 reference genome (Ensembl 110): fasta + BWA index
  (`bwa index` prefix, set via the `genome_index` config key) and
  chromosome sizes file (for the circos plots)
- GRCm39 gene annotation GTF (GENCODE vM33 / Ensembl 110)
  (`GRCm39_GTF_PATH` in `analysis/wes/01_load_data.R`)
- AAV vector sequences (`data/raw/aav/aav_{sgRNA}.fa`)

Paths to these files must be set in the corresponding config / script
headers; see the README for each pipeline.
