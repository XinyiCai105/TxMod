"""End-to-end test on the bundled mini reference: VCF -> transcript-resolved pairs."""
import os
import subprocess
import sys

import pytest

from txmod import prepare_pairs
from txmod.annotation import load_transcripts
from txmod.sequences import GenomeReader, spliced_utr3

HERE = os.path.dirname(__file__)
EX = os.path.abspath(os.path.join(HERE, "..", "examples"))
GTF = os.path.join(EX, "mini_annotation.gtf")
FA = os.path.join(EX, "mini_genome.fa")


@pytest.fixture()
def vcf(tmp_path):
    """A variant inside the shared 3'UTR of both GENEA isoforms."""
    p = tmp_path / "v.vcf"
    tx = load_transcripts(GTF)
    genome = GenomeReader(FA)
    # pick a position in TXA1's first 3'UTR block and read its true REF base
    block = tx["TXA1"].utr3_blocks()[0]
    pos = block.start + 5
    ref = genome.fetch("chr1", pos, pos)
    alt = "A" if ref != "A" else "C"
    p.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"chr1\t{pos}\tv1\t{ref}\t{alt}\t.\t.\t.\n"
    )
    return str(p), pos, ref, alt


def test_spliced_utr3_length_matches_blocks():
    tx = load_transcripts(GTF)
    genome = GenomeReader(FA)
    for tid, t in tx.items():
        assert len(spliced_utr3(t, genome)) == t.utr3_length()


def test_variant_maps_to_multiple_isoforms(vcf):
    """The transcript-resolved core: one variant, several isoform contexts."""
    path, pos, ref, alt = vcf
    pairs = prepare_pairs(path, GTF, FA)
    assert len(pairs) >= 2, "variant in a shared 3'UTR should hit both isoforms"
    assert {p.transcript_id for p in pairs} >= {"TXA1", "TXA2"}


def test_ref_base_matches_genome_and_mut_differs(vcf):
    path, pos, ref, alt = vcf
    for p in prepare_pairs(path, GTF, FA):
        assert p.ref_seq[p.offset_1based - 1] == p.ref_base_transcript
        assert p.mut_seq[p.offset_1based - 1] == p.alt_base_transcript
        assert len(p.ref_seq) == len(p.mut_seq)          # substitution keeps length
        assert p.ref_seq != p.mut_seq


def test_reference_mismatch_is_reported(tmp_path):
    """A wrong REF allele must not be silently accepted."""
    tx = load_transcripts(GTF)
    pos = tx["TXA1"].utr3_blocks()[0].start + 5
    genome = GenomeReader(FA)
    true_ref = genome.fetch("chr1", pos, pos)
    wrong = "A" if true_ref != "A" else "T"
    p = tmp_path / "bad.vcf"
    p.write_text("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                 f"chr1\t{pos}\t.\t{wrong}\tG\t.\t.\t.\n")
    with pytest.warns(UserWarning, match="reference-allele mismatch"):
        pairs = prepare_pairs(str(p), GTF, FA, on_mismatch="skip")
    assert pairs == []
    with pytest.raises(ValueError):
        prepare_pairs(str(p), GTF, FA, on_mismatch="raise")


def test_cli_prepare_runs(tmp_path, vcf):
    """The installed console script works end to end."""
    path, *_ = vcf
    out = tmp_path / "pairs.tsv"
    fa_out = tmp_path / "pairs.fa"
    r = subprocess.run(
        [sys.executable, "-m", "txmod.cli", "prepare",
         "--variants", path, "--gtf", GTF, "--genome", FA,
         "--out", str(out), "--fasta-out", str(fa_out)],
        capture_output=True, text=True,
    )
    assert r.returncode == 0, r.stderr
    assert out.exists() and fa_out.exists()
    body = out.read_text().strip().splitlines()
    assert len(body) >= 3                     # header + >=2 isoform rows
    assert body[0].split("\t")[0] == "record_key"
    assert fa_out.read_text().count(">") >= 4  # REF+MUT per pair
