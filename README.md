# long-term-crispr-in-liver

Analysis code for the study:

> **Long-term in vivo CRISPR/Cas9 gene editing in liver triggers host immunity and clonal mutations in DNA damage response genes beyond off-target effects**

This repository contains the analysis code for three modules:

| Module | Description |
|---|---|
| **Translocation detection** | Detection of CRISPR/Cas9 editing-induced chromosomal translocations from transposon-mediated targeted enrichment sequencing (Tn5 tagmentation + target-specific primer enrichment) |
| **AAV integration detection** | Detection of AAV vector genome integrations from the same targeted sequencing data |
| **WES analysis** | Whole-exome sequencing somatic mutation analysis: filtering of Mutect2/VEP output, mutational burden, gene ranking, and visualization (e.g. lollipop plots) |

## Status

**This repository is being finalized.** The complete reproducible pipeline and data access instructions will be available upon publication.

## Data availability

Raw sequencing data have been deposited in the CNGB Sequence Archive (CNSA) of the China National GeneBank DataBase (CNGBdb) under accession **CNP0008752**. See [`data/README.md`](data/README.md) for details.

## Figure-to-script mapping

*Coming soon.* A table mapping each figure/table in the paper to the script that produced it will be added here upon publication.

## License

MIT — see [LICENSE](LICENSE).
