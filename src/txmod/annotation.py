"""
txmod.annotation
=====================

Transcript annotation handling: parse a GTF/GFF3, build transcript models, and
derive 3'UTR coordinates and spliced sequences.

This is the layer that lets TxMod run from standard reference files
(GTF + genome FASTA) rather than pre-built per-variant FASTA.

A transcript's 3'UTR is the spliced mRNA region downstream of the stop codon:
for a ``+`` strand transcript everything after the last CDS base; for a ``-``
strand transcript everything before the first CDS base in genomic coordinates
(which is downstream in transcript orientation).

Because a 3'UTR may span several exons, transcript-level (spliced) offsets are
maintained alongside genomic coordinates so a genomic variant can be mapped to
its 1-based offset within the spliced 3'UTR.
"""

from __future__ import annotations

import gzip
import re
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

_ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')


def _open_text(path: str):
    """Open plain or gzipped text."""
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def parse_attributes(attr: str) -> Dict[str, str]:
    """Parse a GTF/GFF3 attribute column into a dict.

    Handles the GTF ``key "value";`` convention and the GFF3 ``key=value;``
    convention.
    """
    out: Dict[str, str] = {}
    if "=" in attr and '"' not in attr:
        for part in attr.strip().strip(";").split(";"):
            if not part.strip():
                continue
            if "=" in part:
                k, v = part.split("=", 1)
                out[k.strip()] = v.strip()
        return out
    for k, v in _ATTR_RE.findall(attr):
        out[k] = v
    return out


@dataclass
class Exon:
    """A genomic block, 1-based inclusive."""

    start: int
    end: int


@dataclass
class Transcript:
    """A transcript model with the pieces needed for 3'UTR reconstruction."""

    transcript_id: str
    gene_name: str
    chrom: str
    strand: str
    exons: List[Exon] = field(default_factory=list)
    cds: List[Exon] = field(default_factory=list)

    def sort_features(self) -> None:
        self.exons.sort(key=lambda e: e.start)
        self.cds.sort(key=lambda e: e.start)

    @property
    def has_cds(self) -> bool:
        return len(self.cds) > 0

    def utr3_blocks(self) -> List[Exon]:
        """Genomic blocks of the 3'UTR, sorted by genomic start.

        Empty when the transcript has no CDS (non-coding) or nothing downstream
        of the stop codon.
        """
        if not self.has_cds or not self.exons:
            return []
        self.sort_features()
        blocks: List[Exon] = []
        if self.strand == "+":
            cds_end = self.cds[-1].end
            for e in self.exons:
                if e.end <= cds_end:
                    continue
                start = max(e.start, cds_end + 1)
                if start <= e.end:
                    blocks.append(Exon(start, e.end))
        else:
            cds_start = self.cds[0].start
            for e in self.exons:
                if e.start >= cds_start:
                    continue
                end = min(e.end, cds_start - 1)
                if e.start <= end:
                    blocks.append(Exon(e.start, end))
        return blocks

    def utr3_length(self) -> int:
        return sum(b.end - b.start + 1 for b in self.utr3_blocks())

    def genomic_to_utr3_offset(self, pos: int) -> Optional[int]:
        """Map a 1-based genomic position to a 1-based spliced-3'UTR offset.

        Returns ``None`` when the position is outside the 3'UTR. Offsets run in
        transcript orientation, so on the ``-`` strand offset 1 is the highest
        genomic 3'UTR coordinate.
        """
        blocks = self.utr3_blocks()
        if not blocks:
            return None
        consumed = 0
        if self.strand == "+":
            for b in blocks:  # ascending genomic == transcript order
                if b.start <= pos <= b.end:
                    return consumed + (pos - b.start) + 1
                consumed += b.end - b.start + 1
        else:
            for b in reversed(blocks):  # descending genomic == transcript order
                if b.start <= pos <= b.end:
                    return consumed + (b.end - pos) + 1
                consumed += b.end - b.start + 1
        return None


_REVCOMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    """Reverse complement, preserving case and passing N through."""
    return seq.translate(_REVCOMP)[::-1]


def load_transcripts(
    gtf_path: str,
    transcript_ids: Optional[Iterable[str]] = None,
    coding_only: bool = True,
) -> Dict[str, Transcript]:
    """Parse exon and CDS features from a GTF/GFF3 into transcript models.

    Parameters
    ----------
    gtf_path
        GTF or GFF3 file (optionally gzipped).
    transcript_ids
        If given, keep only these transcript IDs (saves memory on large files).
    coding_only
        Drop transcripts without a CDS or without a 3'UTR.
    """
    keep = set(transcript_ids) if transcript_ids is not None else None
    tx: Dict[str, Transcript] = {}
    with _open_text(gtf_path) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            feature = f[2]
            if feature not in ("exon", "CDS"):
                continue
            attrs = parse_attributes(f[8])
            tid = attrs.get("transcript_id") or attrs.get("Parent") or ""
            if tid.startswith("transcript:"):
                tid = tid.split(":", 1)[1]
            if not tid:
                continue
            if keep is not None and tid not in keep:
                continue
            rec = tx.get(tid)
            if rec is None:
                rec = Transcript(
                    transcript_id=tid,
                    gene_name=attrs.get("gene_name") or attrs.get("gene_id") or "",
                    chrom=f[0],
                    strand=f[6],
                )
                tx[tid] = rec
            block = Exon(int(f[3]), int(f[4]))
            if feature == "exon":
                rec.exons.append(block)
            else:
                rec.cds.append(block)
    for rec in tx.values():
        rec.sort_features()
    if coding_only:
        tx = {k: v for k, v in tx.items() if v.has_cds and v.utr3_blocks()}
    return tx


def build_utr3_index(
    transcripts: Dict[str, Transcript]
) -> Dict[str, List[Tuple[int, int, str]]]:
    """Per-chromosome sorted list of ``(start, end, transcript_id)`` 3'UTR blocks."""
    idx: Dict[str, List[Tuple[int, int, str]]] = defaultdict(list)
    for tid, t in transcripts.items():
        for b in t.utr3_blocks():
            idx[t.chrom].append((b.start, b.end, tid))
    for c in idx:
        idx[c].sort(key=lambda x: x[0])
    return dict(idx)


def transcripts_overlapping(
    index: Dict[str, List[Tuple[int, int, str]]], chrom: str, pos: int
) -> List[str]:
    """All transcript IDs whose 3'UTR covers a 1-based genomic position."""
    blocks = index.get(chrom)
    if not blocks:
        return []
    hits: List[str] = []
    for start, end, tid in blocks:
        if start > pos:
            break
        if start <= pos <= end:
            hits.append(tid)
    return hits
