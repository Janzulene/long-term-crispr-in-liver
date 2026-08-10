from dataclasses import dataclass
from typing import Literal, cast

import pandas as pd
import pysam

from src.kolmogorov_complexity import calc_kolmogorov_complexity
from src.read_pair import ReadFilter, get_paired_reads, partition_read_maps
from src.transform_coordinates import transform_primer_coordinates
from src.utils import get_read_strand, is_overlap, partition_pair_reads


@dataclass
class TransResult:
    """A validated translocation event.

    qname: read name.
    intarget_chr / intarget_strand / intarget_breakpoint_pos:
        on-target locus of the breakpoint.
    intarget_breakpoint_loc: "start" or "end" — which side of the on-target
        alignment the breakpoint lies on.
    offtarget_*: same fields for the off-target locus.
    first_query_breakpoint / second_query_breakpoint:
        breakpoint positions on the read, ordered as in the FASTQ file.
    sequence: forward sequence of the read (same orientation as FASTQ).
    breakpoint_location: "read1" / "read2" — the read carrying the split,
        "between" — the breakpoint lies between the two reads, or "both" —
        both reads carry it.
    intarget_first: whether the on-target segment comes first on the read;
        for "between" events, whether read1 is the forward read.
    """
    qname                   : str
    intarget_chr            : str
    intarget_strand         : Literal["+", "-"]
    intarget_breakpoint_pos : int
    intarget_breakpoint_loc : Literal["start", "end"]
    offtarget_chr           : str
    offtarget_strand        : Literal["+", "-"]
    offtarget_breakpoint_pos: int
    offtarget_breakpoint_loc: Literal["start", "end"]
    first_query_breakpoint  : int
    second_query_breakpoint : int
    sequence                : str
    breakpoint_location     : Literal["read1", "read2", "between", "both"]
    intarget_first          : bool


@dataclass
class BreakPoint:
    first_query_breakpoint        : int
    second_query_breakpoint       : int
    intarget_reference_breakpoint : int
    offtarget_reference_breakpoint: int
    intarget_first                : bool
    intarget_breakpoint_loc       : Literal["start", "end"]
    offtarget_breakpoint_loc      : Literal["start", "end"]


# Detect translocations per primer region.
def detect_translocation(
    primer_chr    : str,
    primer_pos    : int,
    primer_length : int,
    primer_strand : Literal["+", "-"],
    bam_file_path : str,
    read_length   : int   = 150,
    min_mapq      : int   = 0,
    min_complexity: float = 1.5,
    min_map_length: int   = 50
) -> pd.DataFrame:
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

    res_list = []

    for pair_reads in get_paired_reads(bam_file_path, min_mapq):
        read1_reads, read2_reads = partition_pair_reads(pair_reads)

        # TODO [enhancement: handle more cases]: singleton reads
        if (not read1_reads) or (not read2_reads):
            "singleton"
            continue

        # TODO [enhancement: handle more cases]: reads containing Ns
        if (len(read1_reads) > 2) or (len(read2_reads) > 2):
            "contains N"
            continue

        # Classify the mappings of each read
        read1_intarget_map, read1_offtarget_map = partition_read_maps(read1_reads, read_filters)
        read2_intarget_map, read2_offtarget_map = partition_read_maps(read2_reads, read_filters)

        _maps = (read1_intarget_map, read1_offtarget_map, read2_intarget_map, read2_offtarget_map)

        # TODO [enhancement: handle more cases]: inverted reads (two intarget maps)
        if any(map(lambda x: isinstance(x, list), _maps)):
            "invert reads"
            continue

        # Analyse each mapping pattern
        match tuple(map(lambda x: x is not None, _maps)):
            # Both reads are split (two-two pair)
            # --------------------------------
            case True, True, True, True:
                read1_intarget_map  = cast(pysam.AlignedSegment, read1_intarget_map)
                read1_offtarget_map = cast(pysam.AlignedSegment, read1_offtarget_map)
                read2_offtarget_map = cast(pysam.AlignedSegment, read2_offtarget_map)
                read2_intarget_map  = cast(pysam.AlignedSegment, read2_intarget_map)

                read1_breakpoint = find_inner_breakpoint(read1_intarget_map, read1_offtarget_map)
                read2_breakpoint = find_inner_breakpoint(read2_intarget_map, read2_offtarget_map)

                # TODO [enhancement: tolerate small shifts]: accept breakpoints
                # within +/- 3 bp of each other
                # Check that read1 and read2 support the same breakpoint
                same_intarget  = read1_breakpoint.intarget_reference_breakpoint    == read2_breakpoint.intarget_reference_breakpoint
                same_offtarget = read1_breakpoint.offtarget_reference_breakpoint == read2_breakpoint.offtarget_reference_breakpoint
                if not (same_intarget and same_offtarget):
                    continue

                trans_res = TransResult(
                    qname                    = cast(str, read1_intarget_map.query_name),
                    intarget_chr             = cast(str, read1_intarget_map.reference_name),
                    intarget_strand          = get_read_strand(read1_intarget_map),
                    intarget_breakpoint_pos  = read1_breakpoint.intarget_reference_breakpoint,
                    intarget_breakpoint_loc  = read1_breakpoint.intarget_breakpoint_loc,
                    offtarget_chr            = cast(str, read1_offtarget_map.reference_name),
                    offtarget_strand         = get_read_strand(read1_offtarget_map),
                    offtarget_breakpoint_pos = read1_breakpoint.offtarget_reference_breakpoint,
                    offtarget_breakpoint_loc = read1_breakpoint.offtarget_breakpoint_loc,
                    first_query_breakpoint   = read1_breakpoint.first_query_breakpoint,
                    second_query_breakpoint  = read1_breakpoint.second_query_breakpoint,
                    sequence                 = cast(str, read1_intarget_map.get_forward_sequence()),
                    breakpoint_location      = "both",
                    intarget_first           = read1_breakpoint.intarget_first
                )

            # Possible in two-one pairs
            # --------------------------------
            case True, True, False, True:
                # read1 is split, read2 maps fully on-target
                read1_intarget_map  = cast(pysam.AlignedSegment, read1_intarget_map)
                read1_offtarget_map = cast(pysam.AlignedSegment, read1_offtarget_map)
                read2_offtarget_map = cast(pysam.AlignedSegment, read2_offtarget_map)

                # Check that read2's off-target position overlaps read1's
                if not is_overlap(
                    read2_offtarget_map,                  # type: ignore
                    read1_offtarget_map.reference_name,   # type: ignore
                    read1_offtarget_map.reference_start - 1000,  # type: ignore
                    read1_offtarget_map.reference_end   + 1000   # type: ignore
                ):
                    # Mismatched positions: likely an artefact
                    continue

                read1_breakpoint = find_inner_breakpoint(read1_intarget_map, read1_offtarget_map)

                trans_res = TransResult(
                    qname                    = cast(str, read1_intarget_map.query_name),
                    intarget_chr             = cast(str, read1_intarget_map.reference_name),
                    intarget_strand          = get_read_strand(read1_intarget_map),
                    intarget_breakpoint_pos  = read1_breakpoint.intarget_reference_breakpoint,
                    intarget_breakpoint_loc  = read1_breakpoint.intarget_breakpoint_loc,
                    offtarget_chr            = cast(str, read1_offtarget_map.reference_name),
                    offtarget_strand         = get_read_strand(read1_offtarget_map),
                    offtarget_breakpoint_pos = read1_breakpoint.offtarget_reference_breakpoint,
                    offtarget_breakpoint_loc = read1_breakpoint.offtarget_breakpoint_loc,
                    first_query_breakpoint   = read1_breakpoint.first_query_breakpoint,
                    second_query_breakpoint  = read1_breakpoint.second_query_breakpoint,
                    sequence                 = cast(str, read1_intarget_map.get_forward_sequence()),
                    breakpoint_location      = "read1",
                    intarget_first           = read1_breakpoint.intarget_first
                )

            case True, True, True, False:
                # read1 is split, read2 maps fully on-target.
                # Could be background noise amplified onto an existing
                # translocation, or a bridge-PCR error. Skip.
                continue

            # Possible in one-two pairs
            # --------------------------------
            case True, False, True, True:
                # read2 is split
                read1_intarget_map  = cast(pysam.AlignedSegment, read1_intarget_map)
                read2_intarget_map  = cast(pysam.AlignedSegment, read2_intarget_map)
                read2_offtarget_map = cast(pysam.AlignedSegment, read2_offtarget_map)

                # Determine the breakpoint from read2's alignments
                read2_breakpoint = find_inner_breakpoint(read2_intarget_map, read2_offtarget_map)

                # Sanity check
                assert read1_intarget_map.is_forward != read2_intarget_map.is_forward, "Read1 and Read2 should be on different strands."

                trans_res = TransResult(
                    qname                    = cast(str, read2_intarget_map.query_name),
                    intarget_chr             = cast(str, read2_intarget_map.reference_name),
                    intarget_strand          = get_read_strand(read2_intarget_map),
                    intarget_breakpoint_pos  = read2_breakpoint.intarget_reference_breakpoint,
                    intarget_breakpoint_loc  = read2_breakpoint.intarget_breakpoint_loc,
                    offtarget_chr            = cast(str, read2_offtarget_map.reference_name),
                    offtarget_strand         = get_read_strand(read2_offtarget_map),
                    offtarget_breakpoint_pos = read2_breakpoint.offtarget_reference_breakpoint,
                    offtarget_breakpoint_loc = read2_breakpoint.offtarget_breakpoint_loc,
                    first_query_breakpoint   = read2_breakpoint.first_query_breakpoint,
                    second_query_breakpoint  = read2_breakpoint.second_query_breakpoint,
                    sequence                 = cast(str, read2_intarget_map.get_forward_sequence()),
                    breakpoint_location      = "read2",
                    intarget_first           = read2_breakpoint.intarget_first
                )

            case False, True, True, True:
                # read2 is split, read1 maps off-target.
                # Could be background noise amplified onto an existing
                # translocation. Undefined behaviour: skip.
                continue

            # Possible in one-one pairs
            # --------------------------------
            case True, False, True, False:
                # Both reads map fully on-target; no translocation.
                # TODO [enhancement: indel analysis]: this pattern can be
                # used for indel analysis.
                continue

            case True, False, False, True:
                # Breakpoint between the two reads
                read1_intarget_map  = cast(pysam.AlignedSegment, read1_intarget_map)
                read2_offtarget_map = cast(pysam.AlignedSegment, read2_offtarget_map)

                # Reject short off-target mappings
                if read2_offtarget_map.reference_length < min_map_length:  # type: ignore
                    continue

                # Determine the breakpoint from the strand
                if read1_intarget_map.is_forward:
                    target_breakpoint_pos = read1_intarget_map.reference_end
                    breakpoint_query_pos  = read1_intarget_map.query_alignment_end
                else:
                    target_breakpoint_pos = read1_intarget_map.reference_start
                    breakpoint_query_pos  = read1_intarget_map.query_alignment_end

                if read2_offtarget_map.is_forward:
                    offtarget_breakpoint_pos = read2_offtarget_map.reference_start
                else:
                    offtarget_breakpoint_pos = read2_offtarget_map.reference_end

                trans_res = TransResult(
                    qname                    = cast(str, read1_intarget_map.query_name),
                    intarget_chr             = cast(str, read1_intarget_map.reference_name),
                    intarget_strand          = get_read_strand(read1_intarget_map),
                    intarget_breakpoint_pos  = cast(int, target_breakpoint_pos),
                    intarget_breakpoint_loc  = "end" if read1_intarget_map.is_forward else "start",
                    offtarget_chr            = cast(str, read2_offtarget_map.reference_name),
                    offtarget_strand         = get_read_strand(read2_offtarget_map),
                    offtarget_breakpoint_pos = cast(int, offtarget_breakpoint_pos),
                    offtarget_breakpoint_loc = "start" if read2_offtarget_map.is_forward else "end",
                    first_query_breakpoint   = breakpoint_query_pos,
                    second_query_breakpoint  = breakpoint_query_pos,
                    sequence                 = read1_intarget_map.get_forward_sequence(),  # type: ignore
                    breakpoint_location      = "between",
                    intarget_first           = read1_intarget_map.is_forward  # for "between", indicates whether the strand is correct
                )

            case False, True, True, False:
                # Undefined behaviour: read1 maps off-target while read2 maps
                # on-target. Most likely a Tn5-induced artefact, but could
                # also be a real translocation amplified from background
                # noise. Skip.
                continue

            # Singleton-like patterns
            # --------------------------------
            case False, False, _, _:
                # read1 has no usable mapping after filtering.
                # Skip (could potentially be analysed as a single read).
                continue
            case _, _, False, False:
                # read2 has no usable mapping after filtering.
                # Skip (could potentially be analysed as a single read).
                continue
            case _:
                raise ValueError(f"Unexpected situation occurred. {_maps}")

        res_list.append(trans_res)

    trans_res = pd.DataFrame(res_list)
    return trans_res


# Return a sequence oriented with the read (use `get_forward_sequence`).
def find_inner_breakpoint(
    intarget_map: pysam.AlignedSegment,
    offtarget_map: pysam.AlignedSegment
) -> BreakPoint:

    intarget_strand = intarget_map.is_forward
    same_strand     = intarget_map.is_forward == offtarget_map.is_forward

    # Express both alignments on the on-target strand
    if same_strand:
        offtarget_map_query_start = offtarget_map.query_alignment_start
        offtarget_map_query_end   = offtarget_map.query_alignment_end
    else:
        offtarget_map_query_start = offtarget_map.query_length - offtarget_map.query_alignment_end
        offtarget_map_query_end   = offtarget_map.query_length - offtarget_map.query_alignment_start

    intarget_first = intarget_map.query_alignment_start < offtarget_map_query_start

    # The strand combination and the order of the two segments determine
    # which side of each alignment the external breakpoint lies on.
    match same_strand, intarget_first:
        case True, True:
            intarget_reference_breakpoint  = intarget_map.reference_end
            offtarget_reference_breakpoint = offtarget_map.reference_start
            intarget_breakpoint_loc        = "end"
            offtarget_breakpoint_loc       = "start"
        case True, False:
            intarget_reference_breakpoint  = intarget_map.reference_start
            offtarget_reference_breakpoint = offtarget_map.reference_end
            intarget_breakpoint_loc        = "start"
            offtarget_breakpoint_loc       = "end"
        case False, True:
            intarget_reference_breakpoint  = intarget_map.reference_end
            offtarget_reference_breakpoint = offtarget_map.reference_end
            intarget_breakpoint_loc        = "end"
            offtarget_breakpoint_loc       = "end"
        case False, False:
            intarget_reference_breakpoint  = intarget_map.reference_start
            offtarget_reference_breakpoint = offtarget_map.reference_start
            intarget_breakpoint_loc        = "start"
            offtarget_breakpoint_loc       = "start"

    # Determine the internal breakpoint positions on the query
    # ================================
    if intarget_first:
        target_query_breakpoint    = intarget_map.query_alignment_end
        offtarget_query_breakpoint = offtarget_map_query_start
    else:
        target_query_breakpoint    = intarget_map.query_alignment_start
        offtarget_query_breakpoint = offtarget_map_query_end

    # Orient all breakpoints towards the forward read
    if not intarget_strand:
        target_query_breakpoint    = intarget_map.query_length - target_query_breakpoint
        offtarget_query_breakpoint = offtarget_map.query_length - offtarget_query_breakpoint

    # ================================
    if intarget_first:
        first_query_breakpoint  = target_query_breakpoint
        second_query_breakpoint = offtarget_query_breakpoint
    else:
        first_query_breakpoint  = offtarget_query_breakpoint
        second_query_breakpoint = target_query_breakpoint

    # Swap when the read is reverse-oriented
    if not intarget_strand:
        first_query_breakpoint, second_query_breakpoint = second_query_breakpoint, first_query_breakpoint

    breakpoint = BreakPoint(
        first_query_breakpoint         = first_query_breakpoint,
        second_query_breakpoint        = second_query_breakpoint,
        intarget_reference_breakpoint  = cast(int, intarget_reference_breakpoint),
        offtarget_reference_breakpoint = cast(int, offtarget_reference_breakpoint),
        intarget_first                 = intarget_first,
        intarget_breakpoint_loc        = intarget_breakpoint_loc,
        offtarget_breakpoint_loc       = offtarget_breakpoint_loc
    )
    return breakpoint
