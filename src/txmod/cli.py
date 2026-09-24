"""
txmod.cli
==============

Command-line interface.

Subcommands
-----------
``prepare``  variants + GTF + genome -> REF/MUT 3'UTR sequences (TSV or FASTA)
``m6a``      predict m6A tracks and Delta-prob for prepared pairs
``rbp``      scan an RBP PWM panel for binding alterations
``run``      the whole pipeline in one call
``panel``    build a systematic PWM panel from a CISBP-RNA download
"""

from __future__ import annotations

import argparse
import sys
from typing import List, Optional


def _add_variant_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--variants", required=True,
                   help="VCF (.vcf/.vcf.gz) or delimited table of substitutions")
    p.add_argument("--gtf", required=True, help="GTF/GFF3 transcript annotation")
    p.add_argument("--genome", required=True, help="Genome FASTA (indexed if large)")
    p.add_argument("--chrom-col", default="chrom", help="table mode: chromosome column")
    p.add_argument("--pos-col", default="pos", help="table mode: position column")
    p.add_argument("--ref-col", default="ref", help="table mode: REF allele column")
    p.add_argument("--alt-col", default="alt", help="table mode: ALT allele column")
    p.add_argument("--id-col", default=None, help="table mode: variant ID column")
    p.add_argument("--on-mismatch", choices=["skip", "raise"], default="skip",
                   help="what to do when the REF allele disagrees with the genome")


def _variant_kwargs(a: argparse.Namespace) -> dict:
    kw = {"on_mismatch": a.on_mismatch}
    if not str(a.variants).lower().endswith((".vcf", ".vcf.gz", ".bcf")):
        kw.update(chrom_col=a.chrom_col, pos_col=a.pos_col,
                  ref_col=a.ref_col, alt_col=a.alt_col, id_col=a.id_col)
    return kw


def cmd_prepare(a: argparse.Namespace) -> int:
    import pandas as pd

    from .pipeline import prepare_pairs

    pairs = prepare_pairs(a.variants, a.gtf, a.genome, **_variant_kwargs(a))
    if not pairs:
        print("[txmod] no variant-transcript pairs found", file=sys.stderr)
        return 1
    if a.fasta_out:
        with open(a.fasta_out, "w") as fh:
            for p in pairs:
                fh.write(f">{p.record_key}|REF\n{p.ref_seq}\n")
                fh.write(f">{p.record_key}|MUT\n{p.mut_seq}\n")
    df = pd.DataFrame([{
        "record_key": p.record_key, "transcript_id": p.transcript_id,
        "gene_name": p.gene_name, "chrom": p.chrom, "pos": p.pos, "strand": p.strand,
        "ref_allele": p.ref_allele, "alt_allele": p.alt_allele,
        "offset_1based": p.offset_1based, "utr_len": p.utr_len,
    } for p in pairs])
    df.to_csv(a.out, sep="\t", index=False)
    print(f"[txmod] {len(pairs)} variant-transcript pairs -> {a.out}")
    return 0


def cmd_run(a: argparse.Namespace) -> int:
    from .pipeline import run_all

    results, events = run_all(
        a.variants, a.gtf, a.genome,
        model_dir=a.model_dir, pwm_table=a.pwm,
        context=a.context, batch_size=a.batch_size,
        **_variant_kwargs(a),
    )
    if results.empty:
        print("[txmod] no variant-transcript pairs found", file=sys.stderr)
        return 1
    results.to_csv(a.out, sep="\t", index=False)
    print(f"[txmod] {len(results)} pairs -> {a.out}")
    if a.events_out and not events.empty:
        events.to_csv(a.events_out, sep="\t", index=False)
        print(f"[txmod] {len(events)} RBP binding-alteration events -> {a.events_out}")
    return 0


def cmd_panel(a: argparse.Namespace) -> int:
    import pandas as pd

    from .rbp import systematic_panel_from_cisbp

    keep = None
    if a.rbp_list:
        with open(a.rbp_list) as fh:
            keep = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    motifs, skipped = systematic_panel_from_cisbp(a.rbp_information, a.pwm_dir, keep)
    rows = []
    for m in motifs:
        probs = 2.0 ** m.log_matrix * 0.25
        for i in range(m.length):
            rows.append({"RBP": m.rbp, "motif_id": m.motif_id, "pos": i + 1,
                         "A": probs[i, 0], "C": probs[i, 1],
                         "G": probs[i, 2], "U": probs[i, 3]})
    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
    n_rbp = len({m.rbp for m in motifs})
    print(f"[txmod] panel: {n_rbp} RBPs, {len(motifs)} motifs -> {a.out}")
    for reason, ids in skipped.items():
        if ids:
            print(f"[txmod] skipped ({reason}): {len(ids)}", file=sys.stderr)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="txmod",
        description="Transcript-resolved interpretation of 3'UTR variants via "
                    "isoform-specific m6A and RBP-binding alteration.",
    )
    sub = p.add_subparsers(dest="command", required=True)

    pp = sub.add_parser("prepare", help="build REF/MUT 3'UTR sequences")
    _add_variant_args(pp)
    pp.add_argument("--out", required=True, help="output pair table (TSV)")
    pp.add_argument("--fasta-out", default=None, help="also write paired REF/MUT FASTA")
    pp.set_defaults(func=cmd_prepare)

    pr = sub.add_parser("run", help="run the full pipeline")
    _add_variant_args(pr)
    pr.add_argument("--out", required=True, help="output result table (TSV)")
    pr.add_argument("--events-out", default=None, help="per-event RBP table (TSV)")
    pr.add_argument("--model-dir", default=None,
                    help="iM6A model directory; omit to skip m6A prediction")
    pr.add_argument("--pwm", default=None,
                    help="RBP PWM table; omit to skip motif scanning")
    pr.add_argument("--context", type=int, default=10000)
    pr.add_argument("--batch-size", type=int, default=32)
    pr.set_defaults(func=cmd_run)

    pa = sub.add_parser("panel", help="build a PWM panel from a CISBP-RNA download")
    pa.add_argument("--rbp-information", required=True,
                    help="RBP_Information_all_motifs.txt")
    pa.add_argument("--pwm-dir", required=True, help="pwms_all_motifs/ directory")
    pa.add_argument("--rbp-list", default=None,
                    help="file of RBP names to keep (one per line)")
    pa.add_argument("--out", required=True, help="output PWM table (TSV)")
    pa.set_defaults(func=cmd_panel)

    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args) or 0)


if __name__ == "__main__":
    raise SystemExit(main())
