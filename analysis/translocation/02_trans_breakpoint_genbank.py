# ══════════════════════════════════════════════════════════════════════════════
# Export SnapGene-readable GenBank records for translocation breakpoints.
#
# For each translocation event, extract the flanking sequences around the
# on-target and off-target breakpoints from the GRCm39 reference, join
# them into a junction sequence, and write a .gb record with colored
# features (on-target flank, off-target flank, breakpoint) for manual
# inspection in SnapGene.
#
# Article references:
#   - Fig 2L (circos plots): breakpoint validation workflow; the .gb
#     records are used to design PCR validation primers around the
#     detected breakpoints (input: data/final/targetsequence/translocation_stat.tsv,
#     produced by 01_summarize_translocations.py).
#
# External data (not shipped with the repository):
#   - GRCm39 reference fasta (Ensembl 110) — set GRCm39_FASTA below
#
# Usage:
#   python analysis/translocation/02_trans_breakpoint_genbank.py
# ══════════════════════════════════════════════════════════════════════════════

import sys
from pathlib import Path
from typing import Literal

import pandas as pd
import pysam
import siuba as xb
from Bio.Seq import Seq
from Bio.SeqIO import SeqRecord, SeqIO
from Bio.SeqFeature import Reference
from siuba import _

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.snapgene_gb import create_snapgene_feature

# GRCm39 reference fasta (Ensembl 110 primary assembly).
# Download from https://www.ensembl.org and set the path here.
GRCm39_FASTA = "path/to/GRCm39_ensembl.fa"

FINAL_STAT = "data/final/targetsequence/translocation_stat.tsv"
OUTPUT_DIR = "data/processed/translocation_genbank"

mm39_fa = pysam.FastaFile(GRCm39_FASTA)


def get_trans_seq(
        chr                  : str,
        breakpoint_loc       : Literal["start", "end"],
        breakpoint_pos       : int,
        flank_length         : int,
        output_breakpoint_loc: Literal["start", "end"]
    ) -> tuple[Seq, int, int]:
    """Fetch flank_length bp around a breakpoint and orient the sequence
    so that the breakpoint sits at the requested end."""

    if breakpoint_loc == "start":
        start = breakpoint_pos
        end   = breakpoint_pos + flank_length
    elif breakpoint_loc == "end":
        start = breakpoint_pos - flank_length
        end   = breakpoint_pos

    intarget_seq = mm39_fa.fetch(
        reference = chr,
        start     = start,
        end       = end
    )
    seq = Seq(intarget_seq)

    if breakpoint_loc == output_breakpoint_loc:
        return seq, start, end
    else:
        return seq.reverse_complement(), start, end


def trans_to_snapgene(
        trans_record,
        flank_length   : int = 500,
        intarget_color : str = "#b6ffb6",
        offtarget_color: str = "#ff6600"
    ) -> SeqRecord:
    """Build a SnapGene .gb record for one translocation event:
    on-target flank (green) + off-target flank (orange), joined at the
    breakpoint, with a red breakpoint feature in between."""

    intarget_seq, intarget_start, intarget_end = get_trans_seq(
        str(trans_record.intarget_chr),
        trans_record.intarget_breakpoint_loc,
        trans_record.intarget_breakpoint_pos,
        flank_length          = flank_length,
        output_breakpoint_loc = "end"
    )

    offtarget_seq, offtarget_start, offtarget_end = get_trans_seq(
        str(trans_record.offtarget_chr),
        trans_record.offtarget_breakpoint_loc,
        trans_record.offtarget_breakpoint_pos,
        flank_length,
        "start"
    )

    intarget_name = (
        "intarget_"
        f"{trans_record.intarget_chr}_"
        f"{intarget_start}-{intarget_end}"
        f"({trans_record.intarget_strand})"
    )
    intarget_feature = create_snapgene_feature(
        start        = 0,
        end          = flank_length,
        feature_name = intarget_name,
        color        = intarget_color
    )

    offtarget_name = (
        "offtarget_"
        f"{trans_record.offtarget_chr}_"
        f"{offtarget_start}-{offtarget_end}"
        f"({trans_record.offtarget_strand})"
    )
    offtarget_feature = create_snapgene_feature(
        start        = flank_length,
        end          = 2 * flank_length,
        feature_name = offtarget_name,
        color        = offtarget_color
    )

    breakpoint_feature = create_snapgene_feature(
        start        = flank_length - 1,
        end          = flank_length + 1,
        feature_name = "breakpoint",
        color        = "#ff0000"
    )

    reference         = Reference()
    reference.authors = "azulene"
    reference.title   = "Direct Submission"
    reference.journal = "Exported Jan 25, 2025 from SnapGene 4.2.4\nhttp://www.snapgene.com"

    trans_name = (
        f"({trans_record['sample']})_"
        f"({trans_record.region})_"
        f"({trans_record.intarget_chr}:{trans_record.intarget_breakpoint_pos}"
        "->"
        f"{trans_record.offtarget_chr}:{trans_record.offtarget_breakpoint_pos})"
    )
    seq_record = SeqRecord(
        name        = "Exported",
        seq         = intarget_seq + offtarget_seq,
        features    = [intarget_feature, offtarget_feature, breakpoint_feature],
        annotations = {
            "molecule_type": "DNA",
            "topology"     : "linear",
            "references"   : [reference],
            "keywords"     : [trans_name]
        }
    )
    return seq_record


if __name__ == "__main__":
    final_stat = pd.read_table(FINAL_STAT)

    # Example: the second most frequent event in exp_1
    trans_record = (
        final_stat
        >> xb.filter(_["sample"] == "exp_1")
        >> xb.arrange(-_.n)
        >> _.iloc[1]
    )

    out_path = Path(OUTPUT_DIR) / "example_translocation.gb"
    out_path.parent.mkdir(parents = True, exist_ok = True)
    SeqIO.write(trans_to_snapgene(trans_record, 500), out_path, "gb")
    print(f"wrote {out_path}")
