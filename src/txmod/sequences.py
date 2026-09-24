"""
txmod.sequences
====================

Genome sequence access and construction of REF / MUT spliced 3'UTR sequences.

Given a transcript model (:mod:`txmod.annotation`) and a genome FASTA, this
module builds the spliced 3'UTR in transcript orientation and applies a variant
to produce the mutant sequence. These paired REF/MUT sequences are the input to
m6A prediction and RBP motif scanning.

Genome access uses ``pyfaidx`` when available (random access via .fai) and falls
back to a plain in-memory FASTA reader for small references such as the bundled
example.
"""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from typing import Dict, Optional

from .annotation import Transcript, reverse_complement


class GenomeReader:
    """Random-access reader for a genome FASTA.

    Prefers ``pyfaidx`` (memory-efficient, requires a .fai index which pyfaidx
    creates on demand). Falls back to loading sequences into memory, which is
    fine for small test references but not for a whole genome.
    """

    def __init__(self, fasta_path: str):
        self.path = str(fasta_path)
        self._fa = None
        self._mem: Optional[Dict[str, str]] = None
        try:
            import pyfaidx  # type: ignore

            self._fa = pyfaidx.Fasta(self.path, as_raw=True, sequence_always_upper=True)
        except Exception:
            self._mem = self._load_plain(self.path)

    @staticmethod
    def _load_plain(path: str) -> Dict[str, str]:
        opener = gzip.open if path.endswith(".gz") else open
        seqs: Dict[str, str] = {}
        name = None
        buf: list[str] = []
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    if name is not None:
                        seqs[name] = "".join(buf).upper()
                    name = line[1:].strip().split()[0]
                    buf = []
                else:
                    buf.append(line.strip())
        if name is not None:
            seqs[name] = "".join(buf).upper()
        return seqs

    def _candidates(self, chrom: str):
        """Try common chromosome-naming variants (``1`` vs ``chr1``)."""
        c = str(chrom)
        yield c
        if c.startswith("chr"):
            yield c[3:]
        else:
            yield "chr" + c

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """Fetch a 1-based inclusive genomic interval as uppercase sequence."""
        if end < start:
            return ""
        for c in self._candidates(chrom):
            if self._fa is not None:
                try:
                    seq = self._fa[c][start - 1 : end]
                    if seq:
                        return str(seq).upper()
                except Exception:
                    continue
            elif self._mem is not None and c in self._mem:
                return self._mem[c][start - 1 : end].upper()
        raise KeyError(f"contig {chrom!r} not found in {self.path}")

    @property
    def contigs(self):
        if self._fa is not None:
            return list(self._fa.keys())
        return list((self._mem or {}).keys())


def spliced_utr3(transcript: Transcript, genome: GenomeReader) -> str:
    """Return the spliced 3'UTR sequence in transcript orientation."""
    blocks = transcript.utr3_blocks()
    if not blocks:
        return ""
    pieces = [genome.fetch(transcript.chrom, b.start, b.end) for b in blocks]
    seq = "".join(pieces)  # ascending genomic order
    if transcript.strand == "-":
        seq = reverse_complement(seq)
    return seq.upper().replace("U", "T")


@dataclass
class VariantSequences:
    """Paired REF/MUT 3'UTR sequences for one variant-transcript pair."""

    transcript_id: str
    gene_name: str
    chrom: str
    pos: int
    strand: str
    ref_allele: str      # genomic REF allele
    alt_allele: str      # genomic ALT allele
    offset_1based: int   # position within the spliced 3'UTR
    utr_len: int
    ref_seq: str
    mut_seq: str
    ref_base_transcript: str  # base at offset in transcript orientation
    alt_base_transcript: str

    @property
    def record_key(self) -> str:
        """Stable identifier, matching the pipeline's historical convention."""
        return (
            f"{self.transcript_id}|{self.gene_name}|{self.chrom}:{self.pos}"
            f"|REF={self.ref_allele}|ALT={self.alt_allele}"
            f"|strand={self.strand}|offset={self.offset_1based}|len={self.utr_len}"
        )


def apply_variant(
    transcript: Transcript,
    genome: GenomeReader,
    pos: int,
    ref_allele: str,
    alt_allele: str,
    utr3_seq: Optional[str] = None,
    check_reference: bool = True,
) -> Optional[VariantSequences]:
    """Build REF/MUT spliced 3'UTR sequences for a single-nucleotide variant.

    Only substitutions are supported (indels change 3'UTR length and are out of
    scope for position-matched m6A/PWM comparison).

    Returns ``None`` when the variant does not map inside this transcript's
    3'UTR. Raises ``ValueError`` when ``check_reference`` is set and the
    annotated reference base disagrees with the genome.
    """
    if len(ref_allele) != 1 or len(alt_allele) != 1:
        return None
    offset = transcript.genomic_to_utr3_offset(pos)
    if offset is None:
        return None
    ref_seq = utr3_seq if utr3_seq is not None else spliced_utr3(transcript, genome)
    if not ref_seq or offset > len(ref_seq):
        return None

    ref_g = ref_allele.upper().replace("U", "T")
    alt_g = alt_allele.upper().replace("U", "T")
    if transcript.strand == "-":
        ref_t = reverse_complement(ref_g)
        alt_t = reverse_complement(alt_g)
    else:
        ref_t, alt_t = ref_g, alt_g

    observed = ref_seq[offset - 1]
    if check_reference and observed != ref_t:
        raise ValueError(
            f"reference mismatch for {transcript.transcript_id} at {transcript.chrom}:{pos}: "
            f"annotation says {ref_t}, 3'UTR sequence has {observed}"
        )

    mut_seq = ref_seq[: offset - 1] + alt_t + ref_seq[offset:]
    return VariantSequences(
        transcript_id=transcript.transcript_id,
        gene_name=transcript.gene_name,
        chrom=transcript.chrom,
        pos=pos,
        strand=transcript.strand,
        ref_allele=ref_g,
        alt_allele=alt_g,
        offset_1based=offset,
        utr_len=len(ref_seq),
        ref_seq=ref_seq,
        mut_seq=mut_seq,
        ref_base_transcript=ref_t,
        alt_base_transcript=alt_t,
    )
