"""
TxMod
=====

Transcript-resolved interpretation of somatic 3'UTR variants through
isoform-specific epitranscriptomic alteration.

A single genomic variant falls in the 3'UTR of several transcript isoforms and can
have a different consequence in each. TxMod evaluates each variant-transcript
pair separately: it reconstructs the spliced 3'UTR, predicts the change in m6A
methylation and scans for altered RBP binding.

Quick start
-----------
>>> from txmod import prepare_pairs
>>> pairs = prepare_pairs("variants.vcf", "annotation.gtf", "genome.fa")
>>> pairs[0].offset_1based, pairs[0].utr_len      # doctest: +SKIP

Command line::

    txmod prepare --variants v.vcf --gtf a.gtf --genome g.fa --out pairs.tsv
    txmod run --variants v.vcf --gtf a.gtf --genome g.fa \\
        --model-dir models/ --pwm panel.tsv --out results.tsv

Validation status of each layer
-------------------------------
The two layers do not carry equal evidential weight, and the difference matters
when interpreting output:

===============  ==========================================================
Layer            External validation
===============  ==========================================================
m6A (Delta-prob) iM6A itself is published and was assessed by its authors
                 (Luo et al., Nat Commun 2022). What is tested here is that
                 the mutation-associated changes land on measured
                 methylation: predicted losses coincide with GLORI sites at
                 29 times, and with m6A-Atlas v2.0 sites at 17 times, the
                 rate of matched neutral mutations.
RBP binding      Not supported by ENCODE eCLIP: none of the 20 testable
                 panel RBPs reaches q < 0.05 after Benjamini-Hochberg
                 correction. Treat binding alterations as hypotheses.
===============  ==========================================================

Quote enrichment values through
:func:`txmod.rbp.enrichment_with_power` so the site count travels with the
ratio.
"""

from __future__ import annotations

__version__ = "0.1.0"

from .annotation import (
    Exon,
    Transcript,
    build_utr3_index,
    load_transcripts,
    reverse_complement,
    transcripts_overlapping,
)
from .m6a import (
    RRACH,
    WINDOW,
    calibrate_shift,
    REGISTRATION_OFFSET,
    M6AResult,
    motif_change,
    isoform_divergence,
    methylated_base_index,
    overlaps_motif,
    probability_change,
)
from .pipeline import PairResult, prepare_pairs, run_all
from .rbp import (
    MIN_SITES_FOR_ENRICHMENT,
    BindingAlteration,
    Motif,
    best_normalised_score,
    enrichment_with_power,
    load_panel,
    motif_from_probabilities,
    scan_pair,
    systematic_panel_from_cisbp,
)
from .sequences import GenomeReader, VariantSequences, apply_variant, spliced_utr3
from .background import (
    RateModel,
    canonical_context,
    trinucleotide_rate_model,
    depletion_observed_expected,
    size_matched_enrichment,
)
from .variants import Variant, pair_variants_with_transcripts, read_variants

__all__ = [
    "RateModel",
    "canonical_context",
    "trinucleotide_rate_model",
    "depletion_observed_expected",
    "size_matched_enrichment",
    "__version__",
    # annotation
    "Exon", "Transcript", "load_transcripts", "build_utr3_index",
    "transcripts_overlapping", "reverse_complement",
    # sequences
    "GenomeReader", "VariantSequences", "spliced_utr3", "apply_variant",
    # variants
    "Variant", "read_variants", "pair_variants_with_transcripts",
    # m6a
    "REGISTRATION_OFFSET", "RRACH", "WINDOW", "M6AResult", "probability_change",
    "calibrate_shift",
    "methylated_base_index", "overlaps_motif", "motif_change", "isoform_divergence",
    # rbp
    "Motif", "BindingAlteration", "motif_from_probabilities",
    "best_normalised_score", "scan_pair", "load_panel", "systematic_panel_from_cisbp",
    "enrichment_with_power", "MIN_SITES_FOR_ENRICHMENT",
    # pipeline
    "PairResult", "prepare_pairs", "run_all",
]
