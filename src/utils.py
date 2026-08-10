"""Utility functions for read-level analysis."""
from typing import TypeVar, Literal
from operator import attrgetter
from collections.abc import Callable

import pysam

T = TypeVar("T")


def partition(predicate: Callable[[T], bool], values: list[T]) -> tuple[list[T], list[T]]:
    """Partition values into two lists based on the predicate function.

    Args:
        predicate (Callable[[T],bool]): a function to test each value.
        values (list[T]): a list of values.

    Returns:
        tuple[list[T], list[T]]: two lists, the first contains the values
        for which the predicate returns True, the second the rest.
    """
    a: list[T] = []
    b: list[T] = []
    for item in values:
        if predicate(item):
            a.append(item)
        else:
            b.append(item)
    return a, b


def print_read_info(read: pysam.AlignedSegment) -> None:
    print(
        f"read name: {read.query_name}, "
        f"read strand: {'+' if read.is_forward else '-'}, "
        f"read order: {1 if read.is_read1 else 2}, "
        f"read length: {read.query_length}, "
        f"reference name: {read.reference_name}, "
        f"reference start: {read.reference_start}, "
        f"reference end: {read.reference_end}, "
        f"flag: {read.flag}, "
        f"is secondary: {read.is_secondary}, "
        f"is supplementary: {read.is_supplementary}, "
        f"cigar: {read.cigarstring}, "
        f"is mapped: {read.is_mapped}, "
        f"mapping quality: {read.mapping_quality}"
    )


def is_overlap(read: pysam.AlignedSegment, region_chr: str, region_start: int, region_end: int) -> bool:
    if read.reference_name != region_chr:
        return False
    overlap = read.get_overlap(region_start, region_end)

    if overlap is None:
        return False

    return overlap > 0


def partition_pair_reads(pair: list[pysam.AlignedSegment]) -> tuple[list[pysam.AlignedSegment], list[pysam.AlignedSegment]]:
    return partition(attrgetter('is_read1'), pair)


def partition_reads_by_primary(reads: list[pysam.AlignedSegment]) -> tuple[pysam.AlignedSegment, list[pysam.AlignedSegment]]:
    primary_reads, other_reads = partition(
        lambda read: read.is_secondary or read.is_supplementary, reads
    )

    assert len(primary_reads) == 1, "There should be exactly one primary read"
    primary_read = primary_reads[0]
    return primary_read, other_reads


# Both arguments should be mapping positions of the same read.
def calc_intra_overlap(read_map1: pysam.AlignedSegment, read_map2: pysam.AlignedSegment) -> int:
    """Calculate the overlap between two mapping positions of one read.

    Args:
        read_map1 (pysam.AlignedSegment): the first mapping position of a read.
        read_map2 (pysam.AlignedSegment): the second mapping position of the same read.

    Returns:
        int: the overlap value between the two mapping positions.
    """
    assert read_map1.query_name == read_map2.query_name and read_map1.is_read1 == read_map2.is_read1, "Should be the same read."

    read1_start = read_map1.query_alignment_start
    read1_end   = read_map1.query_alignment_end

    read2_start = read_map2.query_alignment_start
    read2_end   = read_map2.query_alignment_end

    return max(0, min(read1_end, read2_end) - max(read1_start, read2_start))


def calc_intra_overlap_rate(read_map1: pysam.AlignedSegment, read_map2: pysam.AlignedSegment) -> tuple[float, float]:
    """Calculate the overlap rate between two mapping positions of one read.

    Args:
        read_map1 (pysam.AlignedSegment): the first mapping position of a read.
        read_map2 (pysam.AlignedSegment): the second mapping position of the same read.

    Returns:
        float: the overlap rate between the two mapping positions.
    """
    intra_overlap = calc_intra_overlap(read_map1, read_map2)
    return intra_overlap / read_map1.query_alignment_length, intra_overlap / read_map2.query_alignment_length


def count_total_reads(bam_file_path: str, min_mapq: int) -> tuple[int, int]:
    read1_count = 0
    read2_count = 0
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        for read in bam_file:
            if read.is_unmapped:
                continue

            if read.is_secondary or read.is_supplementary:
                continue

            # Skip low-quality mappings
            if read.mapping_quality <= min_mapq:
                continue

            if read.is_read1:
                read1_count += 1
            else:
                read2_count += 1
    return read1_count, read2_count


def get_read_strand(read: pysam.AlignedSegment) -> Literal["+", "-"]:
    return "+" if read.is_forward else "-"


def calc_inter_overlap(read1: pysam.AlignedSegment, read2: pysam.AlignedSegment) -> int:
    """Calculate the overlap between two reads (mate pair).

    Args:
        read1 (pysam.AlignedSegment): a read.
        read2 (pysam.AlignedSegment): the other read of the same pair.

    Returns:
        int: the overlap between the two reads.
    """
    assert read1.query_name == read2.query_name, "Should be the same read."

    read1_start = read1.reference_start
    read1_end   = read1.reference_end
    assert read1_end is not None, "Read 1 should be mapped."

    read2_start = read2.reference_start
    read2_end   = read2.reference_end
    assert read2_end is not None, "Read 2 should be mapped."

    return max(0, min(read1_end, read2_end) - max(read1_start, read2_start))
