"""Build SnapGene-readable GenBank (.gb) records with annotated features.

This module provides small helpers for constructing ``SeqFeature`` objects
(genes / misc features / primer binding sites) that SnapGene and similar
GenBank viewers can render with labels, colors, and directions.

Typical use case: visualize CRISPR translocation breakpoints — assemble the
junction sequence (on-target flank + off-target flank) into a ``SeqRecord``,
attach colored features marking each flank and the breakpoint, and write the
record to a ``.gb`` file for manual inspection in SnapGene.

Dependencies: Biopython (``Bio``).

Example
-------
>>> from Bio.Seq import Seq
>>> from Bio.SeqIO import SeqRecord
>>> from snapgene_gb import create_snapgene_feature
>>> record = SeqRecord(Seq("ACGT" * 25), name="junction")
>>> record.features.append(
...     create_snapgene_feature(0, 50, "intarget_chr1:1000-1050", color="#b6ffb6")
... )
>>> record.features.append(
...     create_snapgene_feature(49, 51, "breakpoint", color="#ff0000")
... )
>>> # SeqIO.write(record, "junction.gb", "genbank")
"""

from typing import Literal, TypeAlias, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

_valid_type: TypeAlias = Literal[
    "misc_feature",
    "primer_bind",
]


def create_snapgene_feature(
    start: int,
    end: int,
    feature_name: str,
    strand: int | Literal["+", "-"] | None = None,
    feature_type: _valid_type = "misc_feature",
    **kwargs,
) -> SeqFeature:
    """Create a ``SeqFeature`` rendered by SnapGene.

    Args:
        start: 0-based start position on the sequence.
        end: 0-based end position (exclusive).
        feature_name: label shown in SnapGene.
        strand: "+" / "-" / 1 / -1, or None for unknown strand.
        feature_type: GenBank feature type (e.g. ``misc_feature``, ``primer_bind``).
        **kwargs: extra qualifiers (e.g. ``color="#ff0000"``, ``direction=1``).

    Returns:
        A ``SeqFeature`` suitable for appending to a ``SeqRecord``.
    """
    # Normalize strand to the Biopython convention (1 / -1 / 0).
    match strand:
        case "+":
            strand = 1
        case "-":
            strand = -1
        case None:
            strand = 0
        case _:
            pass

    # SnapGene qualifiers
    qualifiers = {
        "label": feature_name,
        "note": "; ".join(f"{k}: {v}" for k, v in kwargs.items()),
    }

    return SeqFeature(
        location=FeatureLocation(start, end, strand),
        type=feature_type,
        qualifiers=qualifiers,
    )


def create_snapgene_primer(
    start: int,
    end: int,
    name: str,
    sequence: str | Seq,
    strand: int | Literal["+", "-"] | None = None,
    **kwargs,
) -> SeqFeature:
    """Create a ``primer_bind`` feature annotated with its sequence.

    Args:
        start: 0-based start position on the sequence.
        end: 0-based end position (exclusive).
        name: primer name (label).
        sequence: primer sequence (stored in the ``sequence`` qualifier).
        strand: "+" / "-" / 1 / -1, or None.
        **kwargs: extra qualifiers (e.g. ``color``).

    Returns:
        A ``SeqFeature`` of type ``primer_bind``.
    """
    return create_snapgene_feature(
        start=start,
        end=end,
        feature_name=name,
        feature_type="primer_bind",
        color="black",
        sequence=sequence,
        strand=strand,
        **kwargs,
    )


def _create_snapgene_primer(
    primer_seq: str | Seq,
    record: SeqRecord,
    name: str,
    strand: Literal["+", "-"],
) -> SeqFeature:
    """Locate a primer sequence within ``record`` and annotate it.

    The primer is found by exact substring search; for the "-" strand the
    reverse complement is searched. Raises ``ValueError`` if ``strand`` is
    not "+" or "-".
    """
    primer_seq = Seq(primer_seq)

    if strand == "+":
        search_seq = primer_seq
        strand = None  # type: ignore
    elif strand == "-":
        search_seq = primer_seq.reverse_complement()
    else:
        raise ValueError("strand must be either '+' or '-'")

    start = record.seq.find(search_seq)  # type: ignore
    end = start + len(search_seq)

    return create_snapgene_primer(
        start=start,
        end=end,
        name=name,
        sequence=primer_seq,
        strand=strand,
    )


if __name__ == "__main__":
    # Minimal example: build a junction record with two flank features and a
    # breakpoint marker, then print the GenBank text.
    from Bio.Seq import Seq
    from Bio.SeqIO import SeqRecord

    seq = Seq("ACGT" * 50)
    record = SeqRecord(
        seq,
        name="junction_example",
        annotations={"molecule_type": "DNA", "topology": "linear"},
    )
    record.features.append(
        create_snapgene_feature(0, 50, "intarget_flank", color="#b6ffb6")
    )
    record.features.append(
        create_snapgene_feature(49, 51, "breakpoint", color="#ff0000")
    )
    record.features.append(
        create_snapgene_feature(50, 100, "offtarget_flank", color="#ff6600")
    )
    print(record.format("genbank"))
