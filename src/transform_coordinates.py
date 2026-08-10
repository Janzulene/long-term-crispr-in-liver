from typing import Literal


def transform_primer_coordinates(
    primer_pos: int,
    primer_length: int,
    primer_strand: Literal["+", "-"],
    read_length: int
) -> tuple[int, int, int, int]:
    """Transform primer coordinates to region coordinates.

    Args:
        primer_pos: leftmost primer position in the genome, 1-based (NCBI BLAST result).
        primer_length: primer length in bp.
        primer_strand: "+" or "-".
        read_length: read length in bp.

    Returns:
        (region_start, region_end, primer_start, primer_end):
        the query region covered by reads starting at the primer, and the
        true start/end of the primer within that region (on the primer strand).
    """

    # 1-based to 0-based
    primer_pos = primer_pos - 1

    if primer_strand == "+":
        region_start = primer_pos
        region_end   = primer_pos + read_length

        primer_start = primer_pos  # true start of the primer
        primer_end   = primer_pos + primer_length

    elif primer_strand == "-":
        region_start = primer_pos + primer_length - read_length
        region_end   = primer_pos + primer_length

        primer_start = -primer_pos - primer_length  # true start of the primer
        primer_end   = -primer_pos

    return region_start, region_end, primer_start, primer_end
