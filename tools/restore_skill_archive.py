#!/usr/bin/env python
"""Restore the full v4 skill archive from base64 chunk files.

Usage:
    python tools/restore_skill_archive.py

This writes:
    dist/scvelo_cellrank2_fate_mapping_skill_v4_hygiene_revised.zip

Then unpack it with:
    unzip dist/scvelo_cellrank2_fate_mapping_skill_v4_hygiene_revised.zip
"""

from __future__ import annotations

import base64
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHUNK_DIR = ROOT / "dist" / "archive_chunks"
OUT = ROOT / "dist" / "scvelo_cellrank2_fate_mapping_skill_v4_hygiene_revised.zip"


def main() -> None:
    chunks = sorted(CHUNK_DIR.glob("skill_v4.zip.b64.part*"))
    if not chunks:
        raise SystemExit(f"No chunks found in {CHUNK_DIR}")
    data = "".join(p.read_text(encoding="utf-8").strip() for p in chunks)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_bytes(base64.b64decode(data))
    print(f"Wrote {OUT} ({OUT.stat().st_size} bytes) from {len(chunks)} chunks")


if __name__ == "__main__":
    main()
