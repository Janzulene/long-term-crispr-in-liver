# long-term-crispr-in-liver

Analysis code for the study:

> **Long-term in vivo CRISPR/Cas9 gene editing in liver triggers host immunity and clonal mutations in DNA damage response genes beyond off-target effects**

## Repository layout

```
analysis/
  aav/            AAV integration figures and ratio tables (steps 01-02)
  rna/            GO chord diagram of DAVID enrichment (step 01)
  translocation/  translocation summary tables, breakpoint GenBank export,
                  circos plots (steps 01-03)
  wes/            WES analysis pipeline (steps 01-07)
configs/          example samplesheet, primer table and Snakemake config
pipeline/         Snakemake workflows + scripts (targeted sequencing, AAV)
src/              shared modules (Python detector + R utilities and plotting helpers)
data/README.md    data layout and input data expectations
environment.yml   conda environment (single-file dependency definition)
```

Numbered scripts in `analysis/*/` run in ascending order; each step
consumes the outputs of the previous one.

## Environment

```bash
mamba env create -f environment.yml -n long-term-crispr-in-liver
```

## Analyses

### Targeted sequencing (translocations + indels)

Runs the targeted-sequencing pipeline to detect translocations and short indels at
the CRISPR target sites.

```bash
snakemake -s pipeline/target_analyse_2.smk \
    --configfile configs/target_sequence/snakemake_config.yaml \
    --cores 32 --use-conda
```

- Input: raw FASTQ (listed in `configs/target_sequence/samples.tsv`)
- Output: intermediates under `data/processed/targetsequence/`; final
  translocation table at
  `data/final/targetsequence/translocation_stat__cover.filtered.tsv`
- Downstream analysis code: `analysis/translocation/` and `analysis/aav/`

### AAV integration detection

Realigns the targeted-sequencing reads to per-sgRNA AAV vector
references and quantifies AAV integrations at the CRISPR target sites.

```bash
snakemake -s pipeline/target_analyse_aav.smk \
    --configfile configs/target_sequence/snakemake_config.yaml \
    --cores 8 --use-conda
```

- Input: trimmed reads from the targeted-sequencing pipeline; AAV
  vector references
- Output: `data/processed/targetsequence_aav/` and `reports/figures/aav/` (QC plots)
- Downstream analysis code: `analysis/aav/`

### WES somatic mutations

The WES analysis starts from Mutect2 + VEP annotated variant tables,
produced from the raw WES FASTQ by
[nf-core/sarek](https://nf-co.re/sarek) 3.4.2, using
`PBS_Rep3` / `NT_Rep3` as germline references. Arrange the sarek output
as:

```
data/processed/WES2_paired/preprocessing/markduplicates/{sample}/{sample}.depth.txt      # per-sample mean depth (samtools depth)
data/processed/WES2_paired_annotate/annotation/mutect2/{sample}_{control}.mutect2.biallelic_VEP.ann.vcf.tsv   # Mutect2 + VEP tables
```

Example (see the sarek documentation for the full command):

```bash
nextflow run nf-core/sarek -r 3.4.2 --input samplesheet.csv \
    --genome GRCm39 --tools mutect2 --annotate_tools vep --outdir results
```

Downstream analysis code: `analysis/wes/`

### RNA analysis

The RNA-seq raw data were processed with
[nf-core/rnaseq](https://nf-co.re/rnaseq) v3.14.0. Differential
expression was analysed with limma, and key DEGs were functionally
annotated with the DAVID web tool. This script draws the GO chord
diagram from the limma / DAVID result tables shipped in
`data/raw/rna/`.

## Data availability

Raw sequencing data have been deposited in the CNGB Sequence Archive
(CNSA) of the China National GeneBank DataBase (CNGBdb) under accession
**CNP0008752**. See [`data/README.md`](data/README.md) for the expected
local data layout and reference files.

## License

MIT — see [LICENSE](LICENSE).
