"""Short-indel detection on targeted sequencing reads.

Simple version: judge indels directly from the CIGAR strings of the
on-target reads.
Advanced version (not implemented): use mpileup results to cross-check.
"""

from typing import Literal, cast

import pysam

from src.transform_coordinates import transform_primer_coordinates
from src.kolmogorov_complexity import calc_kolmogorov_complexity
from src.read_pair import get_paired_reads, ReadFilter, partition_read_maps
from src.utils import is_overlap, partition_pair_reads, partition


def count_short_indel(
    primer_chr: str,
    primer_pos: int,
    primer_length: int,
    primer_strand: Literal["+", "-"],
    bam_file_path: str,
    read_length: int = 150,
    min_mapq: int = 0,
    min_complexity: float = 1.5,
    min_map_length: int = 50
) -> int:

    region_start, region_end, primer_start, primer_end = transform_primer_coordinates(
        primer_pos, primer_length, primer_strand, read_length
    )

    # Temporary filter functions
    is_intarget      = lambda read: is_overlap(read, primer_chr, region_start - 1000, region_end + 1000)
    is_translocation = lambda read: not is_overlap(read, primer_chr, region_start - 1000, region_end + 1000)
    has_complexity   = lambda read: calc_kolmogorov_complexity(read.query_alignment_sequence) > min_complexity

    read_filters: dict[str, ReadFilter] = {}
    read_filters["is_intarget"]      = is_intarget
    read_filters["has_complexity"]   = has_complexity
    read_filters["is_translocation"] = is_translocation

    event = 0

    for pair_reads in get_paired_reads(bam_file_path, min_mapq):
        read1_reads, read2_reads = partition_pair_reads(pair_reads)

        # TODO [enhancement: handle more cases]: reads containing Ns
        if (len(read1_reads) > 2) or (len(read2_reads) > 2):
            "contains N"
            continue

        if read1_reads:
            read1_intarget_map, read1_offtarget_map = partition_read_maps(read1_reads, read_filters)
        else:
            read1_intarget_map, read1_offtarget_map = None, None

        if read2_reads:
            read2_intarget_map, read2_offtarget_map = partition_read_maps(read2_reads, read_filters)
        else:
            read2_intarget_map, read2_offtarget_map = None, None

        _maps = (read1_intarget_map, read1_offtarget_map, read2_intarget_map, read2_offtarget_map)

        _map_res = tuple(map(lambda x: x is not None, _maps))

        if _map_res == (False, True, True, False):
            # Undefined behaviour: read1 maps off-target while read2 maps
            # on-target. Most likely a Tn5-induced artefact, but could also
            # be a real translocation amplified from background noise.
            # Skip.
            continue

        # TODO: choose the primary alignment, or the longest one
        if isinstance(read1_intarget_map, list):
            read1_intarget_map = choose_read_from_group(read1_intarget_map)

        if isinstance(read2_intarget_map, list):
            read2_intarget_map = choose_read_from_group(read2_intarget_map)

        if (read1_intarget_map and read1_intarget_map.reference_length < min_map_length):  # type: ignore
            read1_intarget_map = None
        if (read2_intarget_map and read2_intarget_map.reference_length < min_map_length):  # type: ignore
            read2_intarget_map = None

        check_1 = is_short_indel(read1_intarget_map) if read1_intarget_map is not None else False
        check_2 = is_short_indel(read2_intarget_map) if read2_intarget_map is not None else False

        if check_1 or check_2:
            event += 1

    return event


def choose_read_from_group(reads: list[pysam.AlignedSegment]) -> pysam.AlignedSegment:
    primary_reads, other_reads = partition(
        lambda read: read.is_secondary or read.is_supplementary, reads
    )

    # Return the primary alignment if there is one
    if primary_reads:
        return primary_reads[0]

    # Otherwise return the longest alignment
    max_length_read = max(
        other_reads,
        key=lambda x: x.reference_length  # type: ignore
    )
    return max_length_read


def count_ins_len(read: pysam.AlignedSegment) -> int:
    """Total insertion length in the read's CIGAR."""
    ins_len = 0
    for op, length in read.cigartuples:  # type: ignore
        if op == 1:  # insertion
            ins_len += length
    return ins_len


def count_del_len(read: pysam.AlignedSegment) -> int:
    """Total deletion length in the read's CIGAR."""
    del_len = 0
    for op, length in read.cigartuples:  # type: ignore
        if op == 2:  # deletion
            del_len += length
    return del_len


def count_ins_and_del_len(read: pysam.AlignedSegment) -> tuple[int, int]:
    """Insertion and deletion lengths from the read's CIGAR."""
    ins_len = 0
    del_len = 0
    for op, length in read.cigartuples:  # type: ignore
        if op == 1:  # insertion
            ins_len += length
        elif op == 2:  # deletion
            del_len += length
    return ins_len, del_len


def is_short_indel(read: pysam.AlignedSegment, max_len: int = 10) -> bool:
    ins_len, del_len = count_ins_and_del_len(read)
    is_short_ins = 0 < ins_len <= max_len
    is_short_del = 0 < del_len <= max_len
    return is_short_ins or is_short_del
