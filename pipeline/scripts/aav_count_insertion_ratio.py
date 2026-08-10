"""
Count the overlap between AAV-related reads and the f/r editing
evidence reads, and compute the AAV insertion ratio.

Called by Snakemake through the `script:` directive; the global
`snakemake` object provides input / output / wildcards.

Output: a four-column TSV
    sgRNA  aav_reads_count  fr_reads_count  intersect_reads_count  aav_insertion_ratio

aav_insertion_ratio = intersect_reads_count / fr_reads_count
"""


def read_set(path: str) -> set[str]:
    """Read a read-name file into a deduplicated set."""
    with open(path) as f:
        return {line.strip() for line in f if line.strip()}


def main() -> None:
    aav_reads_path    = snakemake.input.aav_reads
    f_edit_reads_path = snakemake.input.f_edit_reads
    r_edit_reads_path = snakemake.input.r_edit_reads
    output_path       = snakemake.output.ratio_file
    sgRNA             = snakemake.wildcards.sgRNA

    aav = read_set(aav_reads_path)
    f_set = read_set(f_edit_reads_path)
    r_set = read_set(r_edit_reads_path)

    fr = f_set | r_set
    intersect = aav & fr

    aav_reads_count = len(aav)
    fr_reads_count = len(fr)
    intersect_reads_count = len(intersect)
    ratio = intersect_reads_count / fr_reads_count if fr_reads_count else 0.0

    with open(output_path, "w") as out:
        out.write(
            "sgRNA\taav_reads_count\tfr_reads_count\t"
            "intersect_reads_count\taav_insertion_ratio\n"
        )
        out.write(
            f"{sgRNA}\t{aav_reads_count}\t{fr_reads_count}\t"
            f"{intersect_reads_count}\t{ratio}\n"
        )


if __name__ == "__main__":
    main()
