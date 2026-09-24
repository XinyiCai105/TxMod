"""
txmod.pipeline
===================

End-to-end orchestration: variants + annotation + genome -> transcript-resolved
m6A probability change and RBP binding alteration.

The pipeline is deliberately split so each stage can be run alone (see the
``prepare`` / ``m6a`` / ``rbp`` CLI subcommands). ``run_all`` wires them together
for the common case.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Dict, Iterable, List, Optional, Sequence

from . import m6a as m6a_mod
from . import rbp as rbp_mod
from .annotation import Transcript, load_transcripts
from .sequences import GenomeReader, VariantSequences, apply_variant, spliced_utr3
from .variants import Variant, pair_variants_with_transcripts, read_variants


def prepare_pairs(
    variant_path: str,
    gtf_path: str,
    genome_path: str,
    check_reference: bool = True,
    on_mismatch: str = "skip",
    **variant_kwargs,
) -> List[VariantSequences]:
    """Build REF/MUT 3'UTR sequences for every variant-transcript pair.

    Parameters
    ----------
    on_mismatch
        ``"skip"`` (default) drops pairs whose annotated REF allele disagrees with
        the genome and reports the count; ``"raise"`` fails loudly. Reference
        mismatches usually mean the variant file and the annotation/genome build
        disagree, so they are surfaced rather than silently accepted.
    """
    variants = list(read_variants(variant_path, **variant_kwargs))
    transcripts = load_transcripts(gtf_path)
    pairs = pair_variants_with_transcripts(variants, transcripts)
    genome = GenomeReader(genome_path)

    utr_cache: Dict[str, str] = {}
    out: List[VariantSequences] = []
    n_mismatch = 0
    for v, tid in pairs:
        t = transcripts[tid]
        if tid not in utr_cache:
            utr_cache[tid] = spliced_utr3(t, genome)
        try:
            vs = apply_variant(
                t, genome, v.pos, v.ref, v.alt,
                utr3_seq=utr_cache[tid], check_reference=check_reference,
            )
        except ValueError:
            n_mismatch += 1
            if on_mismatch == "raise":
                raise
            continue
        if vs is not None:
            out.append(vs)
    if n_mismatch:
        import warnings

        warnings.warn(
            f"{n_mismatch} variant-transcript pair(s) dropped on reference-allele "
            f"mismatch; check that the variant file, GTF and genome use the same build.",
            stacklevel=2,
        )
    return out


@dataclass
class PairResult:
    """One transcript-resolved result row."""

    record_key: str
    transcript_id: str
    gene_name: str
    chrom: str
    pos: int
    strand: str
    ref_allele: str
    alt_allele: str
    offset_1based: int
    utr_len: int
    ref_prob_at_mut: Optional[float] = None
    mut_prob_at_mut: Optional[float] = None
    delta_prob: Optional[float] = None
    m6a_class: Optional[str] = None
    motif_disrupted: Optional[bool] = None
    motif_created: Optional[bool] = None
    #: Offset of the scored adenosine relative to the mutation, in -2..+2.
    scored_A_rel: Optional[int] = None
    scored_A_offset_1based: Optional[int] = None
    n_A_in_window: int = 0
    #: Delta-prob read at the mutation index alone (single-index read).
    delta_prob_shift0: Optional[float] = None
    track_shift: Optional[int] = None
    n_rbp_alterations: int = 0


def run_all(
    variant_path: str,
    gtf_path: str,
    genome_path: str,
    model_dir: Optional[str] = None,
    pwm_table: Optional[str] = None,
    context: int = 10000,
    batch_size: int = 32,
    **variant_kwargs,
):
    """Run the full pipeline and return ``(results_df, rbp_events_df)``.

    m6A prediction is skipped when ``model_dir`` is None, and RBP scanning when
    ``pwm_table`` is None, so the pipeline degrades gracefully rather than failing
    when an optional resource is absent.
    """
    import pandas as pd

    pairs = prepare_pairs(variant_path, gtf_path, genome_path, **variant_kwargs)
    results = [
        PairResult(
            record_key=p.record_key,
            transcript_id=p.transcript_id,
            gene_name=p.gene_name,
            chrom=p.chrom,
            pos=p.pos,
            strand=p.strand,
            ref_allele=p.ref_allele,
            alt_allele=p.alt_allele,
            offset_1based=p.offset_1based,
            utr_len=p.utr_len,
        )
        for p in pairs
    ]

    # m6A layer
    if model_dir and pairs:
        models = m6a_mod.load_models(model_dir)
        ref_tracks = m6a_mod.predict_track([p.ref_seq for p in pairs], models,
                                           context=context, batch_size=batch_size)
        mut_tracks = m6a_mod.predict_track([p.mut_seq for p in pairs], models,
                                           context=context, batch_size=batch_size)
        # Measure the track-to-sequence registration on this run's own
        # sequences rather than trusting the module default, and fall back to
        # that default only when there is too little signal to decide.
        try:
            shift = m6a_mod.calibrate_shift(ref_tracks, [p.ref_seq for p in pairs])
        except RuntimeError:
            shift = m6a_mod.REGISTRATION_OFFSET

        for r, p, rt, mt in zip(results, pairs, ref_tracks, mut_tracks):
            res = m6a_mod.probability_change(
                rt, mt, p.offset_1based,
                ref_seq=p.ref_seq, mut_seq=p.mut_seq, shift=shift,
            )
            r.ref_prob_at_mut = res.ref_prob_at_mut
            r.mut_prob_at_mut = res.mut_prob_at_mut
            r.delta_prob = res.delta_prob
            r.m6a_class = res.m6a_class
            r.scored_A_rel = res.scored_A_rel
            r.scored_A_offset_1based = res.scored_A_offset_1based
            r.n_A_in_window = res.n_A_in_window
            r.delta_prob_shift0 = res.delta_prob_shift0
            r.track_shift = res.track_shift
            mc = m6a_mod.motif_change(p.ref_seq, p.mut_seq, p.offset_1based)
            r.motif_disrupted = mc["motif_disrupted"]
            r.motif_created = mc["motif_created"]

    # RBP layer
    event_rows: List[dict] = []
    if pwm_table and pairs:
        motifs = rbp_mod.load_panel(pwm_table)
        for r, p in zip(results, pairs):
            alts = rbp_mod.scan_pair(p.ref_seq, p.mut_seq, motifs)
            r.n_rbp_alterations = len(alts)
            for a in alts:
                event_rows.append({
                    "record_key": p.record_key,
                    "gene_name": p.gene_name,
                    "transcript_id": p.transcript_id,
                    "rbp": a.rbp,
                    "motif_id": a.motif_id,
                    "ref_norm": a.ref_norm,
                    "mut_norm": a.mut_norm,
                    "delta_norm": a.delta_norm,
                    "direction": a.direction,
                })

    results_df = pd.DataFrame([asdict(r) for r in results])
    events_df = pd.DataFrame(event_rows)
    return results_df, events_df
