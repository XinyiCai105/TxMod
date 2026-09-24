"""
txmod.variants
===================

Variant input: read substitutions from a VCF (optionally gzipped) or a plain
TSV/CSV, and pair each variant with every transcript whose 3'UTR contains it.

The pairing step is what makes TxMod transcript-resolved: one genomic
variant typically maps into the 3'UTR of several isoforms, and its predicted
consequence is evaluated separately in each.
"""

from __future__ import annotations

import csv
import gzip
from dataclasses import dataclass
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

from .annotation import Transcript, build_utr3_index, transcripts_overlapping


@dataclass
class Variant:
    """A single-nucleotide variant in genomic coordinates."""

    chrom: str
    pos: int          # 1-based
    ref: str
    alt: str
    variant_id: str = ""

    @property
    def is_snv(self) -> bool:
        return len(self.ref) == 1 and len(self.alt) == 1 and self.ref != self.alt


def _open_text(path: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def read_vcf(path: str, snv_only: bool = True) -> Iterator[Variant]:
    """Yield variants from a VCF.

    Multi-allelic ALT fields are split on commas. Non-substitutions are skipped
    when ``snv_only`` is set (the default) because length-changing variants break
    position-matched REF/MUT comparison.
    """
    with _open_text(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            chrom, pos_s, vid, ref, alt_field = f[0], f[1], f[2], f[3], f[4]
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            for alt in alt_field.split(","):
                v = Variant(chrom=chrom, pos=pos, ref=ref.upper(), alt=alt.upper(),
                            variant_id="" if vid == "." else vid)
                if snv_only and not v.is_snv:
                    continue
                yield v


def read_table(
    path: str,
    chrom_col: str = "chrom",
    pos_col: str = "pos",
    ref_col: str = "ref",
    alt_col: str = "alt",
    id_col: Optional[str] = None,
    snv_only: bool = True,
) -> Iterator[Variant]:
    """Yield variants from a delimited table (TSV or CSV, auto-detected).

    Column names are configurable so the reader accepts COSMIC-style exports and
    other in-house formats without reshaping them first.
    """
    with _open_text(path) as fh:
        sample = fh.read(8192)
        fh.seek(0)
        delim = "\t" if "\t" in sample.splitlines()[0] else ","
        rdr = csv.DictReader(fh, delimiter=delim)
        missing = [c for c in (chrom_col, pos_col, ref_col, alt_col)
                   if c not in (rdr.fieldnames or [])]
        if missing:
            raise ValueError(
                f"{path}: missing required column(s) {missing}; found {rdr.fieldnames}"
            )
        for row in rdr:
            try:
                pos = int(str(row[pos_col]).strip())
            except (ValueError, TypeError):
                continue
            v = Variant(
                chrom=str(row[chrom_col]).strip(),
                pos=pos,
                ref=str(row[ref_col]).strip().upper(),
                alt=str(row[alt_col]).strip().upper(),
                variant_id=str(row.get(id_col, "")).strip() if id_col else "",
            )
            if snv_only and not v.is_snv:
                continue
            yield v


def read_variants(path: str, snv_only: bool = True, **table_kwargs) -> Iterator[Variant]:
    """Dispatch to :func:`read_vcf` or :func:`read_table` by file extension."""
    p = str(path).lower()
    if p.endswith(".vcf") or p.endswith(".vcf.gz") or p.endswith(".bcf"):
        return read_vcf(path, snv_only=snv_only)
    return read_table(path, snv_only=snv_only, **table_kwargs)


def pair_variants_with_transcripts(
    variants: Sequence[Variant],
    transcripts: Dict[str, Transcript],
) -> List[Tuple[Variant, str]]:
    """Pair each variant with every transcript whose 3'UTR contains it.

    Returns a list of ``(variant, transcript_id)`` tuples — the transcript-resolved
    unit of analysis.
    """
    index = build_utr3_index(transcripts)
    pairs: List[Tuple[Variant, str]] = []
    for v in variants:
        for tid in transcripts_overlapping(index, v.chrom, v.pos):
            pairs.append((v, tid))
    return pairs
