"""Repository-local path conventions.

All paths are relative to the repository root (the parent of this
file), so the code works regardless of where the repository is
checked out. Data lives under ``data/`` (raw / processed / final).
"""

from pathlib import Path


class Location(object):
    root:           Path    = Path(__file__).parents[1]
    raw_data:       Path    = root / "data" / "raw"
    processed_data: Path    = root / "data" / "processed"
    final_data:     Path    = root / "data" / "final"
    figures:        Path    = root / "reports" / "figures"
