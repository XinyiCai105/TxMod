"""Tests for 3'UTR derivation and coordinate mapping — the error-prone core."""
import os

import pytest

from txmod.annotation import (Exon, Transcript, load_transcripts,
                                   reverse_complement, transcripts_overlapping,
                                   build_utr3_index, parse_attributes)

EX = os.path.join(os.path.dirname(__file__), "..", "examples")
GTF = os.path.join(EX, "mini_annotation.gtf")


def test_parse_attributes_gtf_and_gff3():
    a = parse_attributes('gene_id "G1"; transcript_id "T1"; gene_name "NAME";')
    assert a["transcript_id"] == "T1" and a["gene_name"] == "NAME"
    b = parse_attributes("ID=t1;Parent=transcript:T2;Name=x")
    assert b["Parent"] == "transcript:T2"


def test_reverse_complement():
    assert reverse_complement("ACGTN") == "NACGT"


def test_plus_strand_utr3_excludes_cds():
    """On + strand the 3'UTR starts one base after the last CDS base."""
    t = Transcript("T", "G", "chr1", "+",
                   exons=[Exon(100, 200), Exon(300, 400)],
                   cds=[Exon(100, 200), Exon(300, 350)])
    blocks = t.utr3_blocks()
    assert blocks == [Exon(351, 400)]
    assert t.utr3_length() == 50


def test_minus_strand_utr3_is_below_cds():
    """On - strand the 3'UTR lies at LOWER genomic coordinates than the CDS."""
    t = Transcript("T", "G", "chr1", "-",
                   exons=[Exon(100, 200), Exon(300, 400)],
                   cds=[Exon(150, 200), Exon(300, 400)])
    blocks = t.utr3_blocks()
    assert blocks == [Exon(100, 149)]


def test_multi_exon_utr3_offset_plus_strand():
    """Offsets accumulate across exons in transcript order."""
    t = Transcript("T", "G", "chr1", "+",
                   exons=[Exon(100, 200), Exon(300, 400)],
                   cds=[Exon(100, 150)])
    # 3'UTR blocks: 151-200 (50 nt) then 300-400 (101 nt)
    assert t.genomic_to_utr3_offset(151) == 1
    assert t.genomic_to_utr3_offset(200) == 50
    assert t.genomic_to_utr3_offset(300) == 51     # jumps the intron
    assert t.genomic_to_utr3_offset(250) is None   # intronic
    assert t.utr3_length() == 151


def test_multi_exon_utr3_offset_minus_strand():
    """On - strand offset 1 is the HIGHEST genomic 3'UTR coordinate."""
    t = Transcript("T", "G", "chr1", "-",
                   exons=[Exon(100, 200), Exon(300, 400)],
                   cds=[Exon(350, 400)])
    # 3'UTR blocks genomic: 100-200 and 300-349; transcript order starts at 349
    assert t.genomic_to_utr3_offset(349) == 1
    assert t.genomic_to_utr3_offset(300) == 50
    assert t.genomic_to_utr3_offset(200) == 51     # continues into lower block
    assert t.genomic_to_utr3_offset(100) == 151


def test_noncoding_transcript_has_no_utr3():
    t = Transcript("T", "G", "chr1", "+", exons=[Exon(100, 200)], cds=[])
    assert t.utr3_blocks() == []
    assert t.genomic_to_utr3_offset(150) is None


def test_load_example_gtf_and_isoforms():
    tx = load_transcripts(GTF)
    assert {"TXA1", "TXA2", "TXB1"} <= set(tx)
    # TXA1 has a two-block 3'UTR (1651-1700 and 2050-2500)
    assert len(tx["TXA1"].utr3_blocks()) == 2
    # the two GENEA isoforms have different 3'UTR lengths
    assert tx["TXA1"].utr3_length() != tx["TXA2"].utr3_length()
    assert tx["TXB1"].strand == "-"


def test_overlap_index_finds_all_isoforms():
    """A position in a shared 3'UTR region returns every containing isoform."""
    tx = load_transcripts(GTF)
    idx = build_utr3_index(tx)
    hits = transcripts_overlapping(idx, "chr1", 1660)
    assert "TXA1" in hits and "TXA2" in hits   # shared UTR region
    assert transcripts_overlapping(idx, "chr1", 999999) == []
