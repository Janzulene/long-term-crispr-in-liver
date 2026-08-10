# long-term-crispr-in-liver

Analysis code for the study:

> **Long-term in vivo CRISPR/Cas9 gene editing in liver triggers host immunity and clonal mutations in DNA damage response genes beyond off-target effects**

This repository contains the analysis code for three modules:

| Module | Description |
|---|---|
| **Translocation detection** | Detection of CRISPR/Cas9 editing-induced chromosomal translocations from transposon-mediated targeted enrichment sequencing (Tn5 tagmentation + target-specific primer enrichment) |
| **AAV integration detection** | Detection of AAV vector genome integrations from the same targeted sequencing data |
| **WES analysis** | Whole-exome sequencing somatic mutation analysis: filtering of Mutect2/VEP output, mutational burden, gene ranking, and visualization (e.g. lollipop plots) |

## Repository layout

```
analysis/
  aav/            AAV integration figures and ratio tables
  translocation/  translocation summary tables, breakpoint GenBank export, circos plots
  wes/            WES analysis pipeline (6 steps) + lollipop plot
configs/          example samplesheet, primer table and Snakemake config
pipeline/         Snakemake workflows + scripts (targeted sequencing, AAV)
src/              shared modules (Python detector + R utilities)
data/README.md    data layout and input data expectations
environment.yml   conda environment (single-file dependency definition)
```

## Environment

All dependencies are defined in a single conda environment file:

```bash
mamba env create -f environment.yml -n long-term-crispr-in-liver
```

## Pipeline: targeted sequencing (translocations + indels)

```bash
snakemake -s pipeline/target_analyse_2.smk \
    --configfile configs/target_sequence_240501/snakemake_config.yaml \
    --cores 32 --use-conda
```

Steps: fastp trimming → BWA-MEM alignment to GRCm39 → on-target read
filtering (primer match, strand, complexity) → translocation detection →
short-indel counting → final translocation table
(`data/final/targetsequence_20240501/translocation_stat__cover.filtered.tsv`).

The GRCm39 reference (Ensembl 110) is external data: set the BWA index
path via the `genome_index` config key (see
`configs/target_sequence_240501/snakemake_config.yaml`).

## Pipeline: AAV integration detection

```bash
snakemake -s pipeline/target_analyse_aav.smk \
    --configfile configs/target_sequence_240501/snakemake_config.yaml \
    --cores 8 --use-conda
```

Realigns the target-sequencing reads to per-sgRNA AAV vector references
(`data/raw/aav/aav_{sgRNA}.fa`, external data) and quantifies AAV
integrations at the CRISPR target sites.

## Analysis: WES somatic mutations

The six steps are run in order from the repository root:

```bash
Rscript analysis/wes/01_load_data.R   # depth + reference + VCF loading
Rscript analysis/wes/02_somatic_filter.R  # somatic + hard filtering
Rscript analysis/wes/03_vep_annotate.R    # VEP annotation explosion
Rscript analysis/wes/04_mutation_burden.R # burden table (Table S4)
Rscript analysis/wes/05_gene_ranking.R    # per-gene ranking (Fig 2I, Table S5)
Rscript analysis/wes/06_lollipop_data.R   # lollipop data tables
Rscript analysis/wes/plot_lollipop.R      # lollipop figure
```

Intermediate results are stored as RDS files under
`data/processed/WES_20250105_analysis/`; final tables and figures are
written to `data/final/WES_20250105/` and `reports/figures/`.

## Analysis: translocation / AAV result tables

```bash
python analysis/translocation/01_summarize_translocations.py
python analysis/translocation/02_trans_breakpoint_genbank.py
Rscript analysis/translocation/plot_translocation_circos.R   # Fig 2L
Rscript analysis/aav/summarize_aav_insertion_ratio.R          # Fig 2K tables
Rscript analysis/aav/plot_aav.R                               # Fig 2K
```

## Data availability

Raw sequencing data have been deposited in the CNGB Sequence Archive
(CNSA) of the China National GeneBank DataBase (CNGBdb) under accession
**CNP0008752**. See [`data/README.md`](data/README.md) for the expected
local data layout and external reference files (GRCm39 genome, GTF,
assembly chain, AAV vector sequences).

## Figure-to-script mapping

*Coming soon.* A table mapping each figure/table in the paper to the
script that produced it will be added here upon publication.

## License

MIT — see [LICENSE](LICENSE).
