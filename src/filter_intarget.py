from dataclasses import dataclass
from typing import Literal, cast

import numpy as np
import pysam

from src.transform_coordinates import transform_primer_coordinates


@dataclass
class FilterStat:
    n_total: int
    n_same_strand: int
    n_perfect: int
    n_cover: int
    n_left_sub: int
    n_right_sub: int
    n_inner: int
    n_start_within4: int


def find_intarget_reads(
    primer_chr: str,
    primer_pos: int,
    primer_length: int,
    primer_strand: Literal["+", "-"],
    bam_file_path: str,
    strategy: Literal["only_perfect_read1", "only_perfect", "cover", "start_within4"] = "only_perfect_read1",
    read_length: int = 150
) -> tuple[FilterStat, np.ndarray]:
    """Find reads that are anchored at the target (primer) region.

    Args:
        primer_chr: chromosome of the primer.
        primer_pos: leftmost primer position in the genome, 1-based (NCBI BLAST result).
        primer_length: primer length in bp.
        primer_strand: "+" or "-".
        bam_file_path: path to the indexed BAM file.
        strategy: which reads to keep:
            - "only_perfect_read1": perfect primer match, read 1 only
            - "only_perfect": perfect primer match
            - "cover": reads fully covering the primer
            - "start_within4": reads starting within +/- 4 bp of the primer
        read_length: read length in bp.

    Returns:
        (FilterStat, np.ndarray): per-category counts and the unique names of
        the retained reads.
    """

    # Input validation
    assert isinstance(primer_chr, str), f"chromosome should be string, not {type(primer_chr)}!"
    assert primer_strand in ("+", "-"), f"strand of primer should be `+` or `-`, not {primer_strand}!"

    bam_file = pysam.AlignmentFile(bam_file_path, "rb")
    bam_file.check_index()

    region_start, region_end, primer_start, primer_end = transform_primer_coordinates(
        primer_pos, primer_length, primer_strand, read_length
    )

    n_total         = 0
    n_same_strand   = 0
    n_perfect       = 0
    n_cover         = 0
    n_left_sub      = 0
    n_right_sub     = 0
    n_inner         = 0
    n_start_within4 = 0

    read_name_list = []
    for read in bam_file.fetch(primer_chr, region_start, region_end):
        # Skip unmapped reads
        if read.is_unmapped:
            continue
        n_total += 1

        # Keep reads whose strand matches the primer strand
        read_strand = "+" if read.is_forward else "-"
        if read_strand != primer_strand:
            continue
        n_same_strand += 1

        # Express the alignment on the primer strand
        if primer_strand == "+":
            align_start = read.reference_start
            align_end   = cast(int, read.reference_end)  # mapped reads have an end position
        else:
            align_start = -cast(int, read.reference_end)
            align_end   = -read.reference_start

        match align_start, align_end:
            case (align_start, align_end) if (align_start == primer_start) & (align_end >= primer_end):
                # Read starts exactly at the primer and covers it fully
                align_type = "perfect match"
                n_perfect += 1
            case (align_start, align_end) if (align_start < primer_start) & (align_end >= primer_end):
                # Read covers the primer completely
                align_type = "cover match"
                n_cover += 1
            case (align_start, align_end) if (align_start > primer_start) & (align_end >= primer_end):
                # Read overlaps the left part of the primer
                align_type = "left sub match"
                n_left_sub += 1
            case (align_start, align_end) if (align_start < primer_start) & (align_end < primer_end):
                # Read overlaps the right part of the primer
                align_type = "right sub match"
                n_right_sub += 1
            case (align_start, align_end) if (align_start >= primer_start) & (align_end < primer_end):
                # Read is fully inside the primer
                align_type = "inner match"
                n_inner += 1
            case _:
                raise AssertionError(
                    f"Unexpected read/primer alignment: read:{align_start}-{align_end}, primer:{primer_start}-{primer_end}"
                )

        start_within4 = False
        if primer_start - 4 <= align_start <= primer_start + 4:
            start_within4 = True
            n_start_within4 += 1

        match strategy:
            case "only_perfect_read1":
                if align_type != "perfect match":
                    continue
                if read.is_read2:
                    continue
            case "only_perfect":
                if align_type != "perfect match":
                    continue
            case "cover":
                if align_type not in ["cover match", "perfect match"]:
                    continue
            case "start_within4":
                if not start_within4:
                    continue
            case _:
                raise AssertionError(
                    f"strategy should be one of ['only_perfect_read1', 'only_perfect', 'cover', 'start_within4'], not {strategy}"
                )

        read_name_list.append(read.query_name)

    bam_file.close()

    read_name_list = np.unique(read_name_list)
    filter_stat = FilterStat(
        n_total         = n_total,
        n_same_strand   = n_same_strand,
        n_perfect       = n_perfect,
        n_cover         = n_cover,
        n_left_sub      = n_left_sub,
        n_right_sub     = n_right_sub,
        n_inner         = n_inner,
        n_start_within4 = n_start_within4
    )
    return filter_stat, read_name_list


# After filtering, extract the private BAM and sort it by name:
# samtools view -@ 16 -N names.txt input.bam | samtools sort -n -@ 16 -o sorted.bam -
