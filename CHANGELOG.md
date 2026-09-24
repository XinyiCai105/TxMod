# Changelog

## 0.1.0

First public release, the version used in the accompanying manuscript.

- **Inputs.** Standard VCF/TSV variants, GTF/GFF3 annotation and genome FASTA.
  Multi-exon and minus-strand 3′ UTRs are handled, and each variant is paired
  with every isoform whose 3′ UTR contains it (`txmod.annotation`,
  `txmod.sequences`, `txmod.variants`).
- **m6A layer** (`txmod.m6a`). Wrapper for the iM6A ensemble, with the
  one-nucleotide registration of the probability track made explicit and
  checked at run time. Δprob is read over a ±2 nt window around the
  substitution, and the single-index value is reported alongside it. RRACH
  disruption and creation calls, and inter-isoform divergence.
- **RBP layer** (`txmod.rbp`). PWM log-odds scoring with normalised scores,
  high-confidence binding-alteration calling, panel construction from a
  CISBP-RNA download restricted to a user-supplied protein list (skipped motifs
  are reported), and `enrichment_with_power()` so that site counts travel with
  enrichment ratios.
- **Matched backgrounds** (`txmod.background`). Trinucleotide somatic
  mutation-rate model with observed-to-expected depletion testing, and gene-set
  enrichment corrected for size and opportunity confounders by logistic
  regression. Both return the uncorrected and the corrected estimate together.
- **Pipeline and CLI** (`txmod.pipeline`, `txmod.cli`). `txmod prepare`,
  `txmod run` and `txmod panel`. TensorFlow, pyfaidx and statsmodels are
  optional extras, so the core installs without a deep-learning stack.
- **Manuscript figures** (`paper/`). Scripts that render every figure of the
  accompanying manuscript from the processed tables deposited at Zenodo.

Checked against the pipeline that produced the manuscript results: RRACH
annotation 100 % concordant over 20,000 records, normalised PWM scores
identical to within 3 × 10⁻¹⁶, and the Cancer Gene Census enrichment reproduced
from the deposited gene-level table (odds ratio 1.955 uncorrected, 1.119
corrected). 46 tests.
