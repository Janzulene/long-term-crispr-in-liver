"""Read-pair handling for translocation detection."""
from typing import Optional, TypeAlias
from collections.abc import Iterator, Callable

import pysam

from src.utils import partition

ReadFilter: TypeAlias = Callable[[pysam.AlignedSegment], bool]


def get_paired_reads(bam_file_path: str, min_mapq: int) -> Iterator[list[pysam.AlignedSegment]]:
    with pysam.AlignmentFile(bam_file_path, "rb") as bam_file:
        pair_reads: list[pysam.AlignedSegment] = []
        current_qname = None  # drop the group of the first read

        for read in bam_file:
            if read.query_name != current_qname:
                if pair_reads:
                    yield pair_reads

                # Reset the group and update current_qname
                pair_reads = []
                current_qname = read.query_name

            # Skip unmapped reads
            if read.is_unmapped:
                continue

            # Skip low-quality mappings
            if read.mapping_quality <= min_mapq:
                continue

            pair_reads.append(read)

        if pair_reads:
            yield pair_reads


# Analyse the mappings of one read: identify the on-target (intarget) and
# off-target maps.
# TODO [enhancement: stricter classification]: use the intarget map to decide
# whether an off-target map is a secondary or supplementary alignment.
def partition_read_maps(read_maps: list[pysam.AlignedSegment], read_filters: dict[str, ReadFilter]) -> tuple[Optional[pysam.AlignedSegment], Optional[pysam.AlignedSegment]]:
    intarget_maps, offtarget_maps = partition(read_filters["is_intarget"], read_maps)

    # No complexity test for intarget maps: the target region is complex by design
    if len(intarget_maps) == 0:
        intarget_map = None
    elif len(intarget_maps) == 1:
        intarget_map = intarget_maps[0]
    else:
        # Usually a long deletion
        # TODO [enhancement: handle more cases]: inverted reads (two intarget
        # maps) need a dedicated flag; for now return the list
        intarget_map = intarget_maps  # type: ignore

    # Skip complexity calculation for non-offtarget maps to save time
    offtarget_maps = list(filter(read_filters["has_complexity"], offtarget_maps))

    if len(offtarget_maps) == 0:
        offtarget_map = None
    elif len(offtarget_maps) == 1:
        offtarget_map = offtarget_maps[0]
    else:
        offtarget_map = offtarget_maps  # type: ignore

    return intarget_map, offtarget_map  # type: ignore
