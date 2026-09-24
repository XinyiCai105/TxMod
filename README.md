# TxMod

**Transcript-resolved interpretation of 3′UTR variants through isoform-specific
epitranscriptomic alteration.**

A genomic variant in a 3′UTR usually falls inside the 3′UTR of *several* transcript
isoforms, and its consequence can differ in each — a variant that destroys an m6A
site in one isoform may be inert in another because the surrounding sequence, or
the 3′UTR itself, is different. Tools that annotate a variant once per gene cannot
express this. TxMod evaluates every **variant × transcript** pair separately.

For each pair it: reconstructs the spliced 3′UTR, predicts the change in m6A
methylation (Δprob) with the iM6A ensemble, scans an RBP motif panel for altered
binding potential.

## Install

TxMod is not yet on PyPI. Install it from this repository:

```bash
pip install "git+<repo-url>"              # core: annotation, sequences, motifs
pip install "txmod[m6a] @ git+<repo-url>" # + TensorFlow, for iM6A prediction
pip install "txmod[all] @ git+<repo-url>" # + pyfaidx and all optional extras
```

From a clone:

```bash
git clone <repo-url> && cd txmod
pip install -e ".[dev]"
pytest
```

## Quick start

Everything below runs on the bundled miniature reference in `examples/` in a few
seconds — no genome download, no GPU.

```bash
txmod prepare \
    --variants examples/demo.vcf \
    --gtf      examples/mini_annotation.gtf \
    --genome   examples/mini_genome.fa \
    --out      pairs.tsv \
    --fasta-out pairs.fa
```

```
record_key                                                   transcript_id  gene_name  offset_1based  utr_len
TXA1|GENEA|chr1:1654|REF=C|ALT=A|strand=+|offset=4|len=501    TXA1           GENEA      4              501
TXA2|GENEA|chr1:1654|REF=C|ALT=A|strand=+|offset=4|len=250    TXA2           GENEA      4              250
```

One variant, two isoform contexts, **different 3′UTR lengths** — that is the point
of the tool.

Full pipeline (m6A + RBP layers):

```bash
txmod run \
    --variants variants.vcf --gtf gencode.gtf --genome GRCh38.fa \
    --model-dir /path/to/im6a_models --pwm panel.tsv \
    --out results.tsv --events-out rbp_events.tsv
```

`--model-dir` and `--pwm` are optional: omit either and that layer is skipped, so
you can run the annotation/sequence layer with no deep-learning stack installed.

Build a motif panel systematically from a CISBP-RNA download rather than by hand:

```bash
txmod panel \
    --rbp-information RBP_Information_all_motifs.txt \
    --pwm-dir pwms_all_motifs/ \
    --rbp-list stability_rbps.txt \
    --out panel.tsv
```

Motifs whose PWM file is missing or empty are **reported, not silently dropped**.

## Python API

```python
from txmod import prepare_pairs, load_panel, scan_pair

pairs = prepare_pairs("variants.vcf", "gencode.gtf", "GRCh38.fa")
motifs = load_panel("panel.tsv")

for p in pairs:
    events = scan_pair(p.ref_seq, p.mut_seq, motifs)
    if events:
        s = score_pair([{"rbp": e.rbp, "delta_norm": e.delta_norm} for e in events])
        print(p.transcript_id, p.offset_1based, s.tendency)
```

## Inputs

| Input | Format | Notes |
|---|---|---|
| Variants | VCF (`.vcf`/`.vcf.gz`) or TSV/CSV | Substitutions only; column names configurable for tables |
| Annotation | GTF or GFF3 (optionally gzipped) | Needs `exon` **and** `CDS` features to locate the stop codon |
| Genome | FASTA | `pyfaidx` used when available; `chr1`/`1` naming both accepted |
| m6A models | iM6A `.h5` ensemble | Optional; `humanRRACH10000_c{1..5}.h5` |
| RBP panel | PWM table | Optional; from `txmod panel` or your own |

Indels are skipped by design: they change 3′UTR length, so position-matched
REF/MUT comparison of a per-nucleotide track is not meaningful.

## Validation status differs by layer

The two layers do not carry equal evidential weight. Read this before quoting
any output.

| Layer | External validation |
|---|---|
| **m6A (Δprob)** | **The model is published and was assessed by its own authors; this package does not re-benchmark it.** iM6A's evaluation against m6A-CLIP across whole transcriptomes is reported in the publication that released it (Luo et al., *Nature Communications* 2022, 13: 2720). What the accompanying study establishes is a different claim — that the *mutation-associated changes* called here land on measured methylation: predicted losses coincide with a GLORI-measured site in **41.9 %** of cases against **1.4 %** of matched neutral mutations (**29-fold**, n = 7,645), and **17-fold** against m6A-Atlas v2.0. Gains behave the same way where the scored adenosine pre-exists (10.9-fold); where the mutation creates the adenosine, untreated-cell references cannot report it by construction. |
| **RBP binding** | **Not supported by the orthogonal test.** Of the 20 panel RBPs with released GRCh38 ENCODE eCLIP data, **none reaches q < 0.05** after Benjamini-Hochberg correction across the 20 tests. The strongest is FXR1 at **3.27-fold** (7/25 events in peaks vs 8.6 % background, p = 0.005, q = 0.100); IGF2BP1 is **1.29-fold** (14/61, p = 0.30, q = 0.87). Treat predicted binding alterations as hypotheses to test, not as findings. |

Quote enrichment values through `enrichment_with_power()` so the site count
travels with the ratio — a 2-site observation can read as "2.1× enriched" while
the Fisher test returns p = 0.23. Enrichments resting on fewer than 20 sites are
flagged `underpowered` and `headline_safe=False`.

## Two things to know before interpreting output

**1. Δprob is read over a window, not at one index.** A substitution can alter
methylation of any adenosine whose motif context it touches. An RRACH is five
nucleotides wide, so the affected methylated adenosine lies within ±2 nt of the
mutation; `probability_change()` scans that window and reports the adenosine that moves
most, along with `scored_A_rel` (−2…+2).

Reading a single index instead covers one motif position and silently drops the
other four. The difference is large: on the COSMIC 3′UTR set a single-index read
gives **2,299 events at a 4.19:1 loss:gain ratio**, the window read gives
**14,274 at 1.15:1** — so that imbalance is an artefact of *where the value is
read*, not a property of the mutations. Every result carries `delta_prob_shift0`
(the single-index value) so the two reads can be compared on the same run.

The track-to-sequence offset itself is measured, not assumed: `calibrate_shift()`
asks which shift places high-probability calls on the adenosine of an RRACH, and
raises rather than guessing when the evidence is ambiguous. `run_all()` calibrates
on your own sequences and falls back to `REGISTRATION_OFFSET` only if there is
too little signal to decide.

**2. `RRACH` is the motif.** It matches the training definition of the
`humanRRACH10000` ensemble. A different consensus can be supplied through
`motif=`, but pentamers outside the training definition are scored highly only
once a mutation converts them to RRACH, so they return gains and no losses. On
the reference dataset the six `D=U` pentamers of the broader DRACH consensus
behaved exactly that way — **181 events, all gains** — which is why RRACH is
what the package ships and what the study reports.



Reference-allele mismatches are surfaced too: if the variant file disagrees with
the genome, affected pairs are dropped with a warning (`--on-mismatch raise` to
fail hard). That almost always means mismatched genome builds.

## Modules

| Module | Role |
|---|---|
| `annotation` | GTF/GFF3 parsing, 3′UTR blocks, genomic ↔ spliced-offset mapping |
| `sequences` | Genome access, spliced 3′UTR, REF/MUT construction |
| `variants` | VCF/table readers, variant × transcript pairing |
| `m6a` | iM6A ensemble, windowed Δprob, shift calibration, RRACH motif calls, isoform divergence |
| `rbp` | PWM log-odds scoring, binding-alteration calling, panel building |
| `pipeline` | `prepare_pairs`, `run_all` orchestration |
| `cli` | `prepare` / `run` / `panel` subcommands |

## Reproducibility

The scoring logic reproduces the published pipeline that produced the original
results, so numbers computed with this package are directly comparable:

- RRACH annotation (overlap, disruption, creation): **100 % concordant** on 20,000 records
- Windowed Δprob: reproduces the published v2 event set (14,274 events); `delta_prob_shift0` reproduces the v1 set (2,299) to floating-point precision
- PWM normalised scores: identical to **3×10⁻¹⁶** (floating-point precision)

`pytest` (46 tests) covers minus-strand and multi-exon offset mapping, the
measured registration offset, the ±2 nt window read (including a case the
single-index read misses), RRACH motif calls, PWM normalisation bounds,
event thresholds, coefficient-scheme disagreement, and the
underpowered-enrichment guard.

For the biological validation of each layer, see *Validation status differs by
layer* above.

## Reproducing the manuscript figures

`paper/` holds the scripts that render every figure in the accompanying
manuscript from the processed tables deposited at Zenodo. See
[`paper/README.md`](paper/README.md).

## Citation

If you use TxMod, please cite the accompanying manuscript

## Licence

MIT — see `LICENSE`.

## Matched backgrounds (`txmod.background`)

Three results in the accompanying study changed once the comparison was drawn
against a *matched* background instead of a raw one, and two reversed. Those
comparisons are provided as reusable functions, because the failure mode is
generic: an apparent signal that is really a property of what the cases are made of.

| Function | Guards against | Study finding |
|---|---|---|
| `trinucleotide_rate_model` + `depletion_observed_expected` | sequence composition read as selection | apparent mutation depletion around m6A sites vanished (corrected ratio 1.00-1.01 across thresholds 0.3-0.8, all CIs spanning 1.0) |
| `size_matched_enrichment` | gene size read as biology | apparent cancer-gene enrichment vanished (OR 1.96 to 1.12, 95% CI 0.93-1.35, once adjusted) |

Each returns the uncorrected and corrected estimate together, so neither can be
reported without the other.

```python
from txmod.background import size_matched_enrichment

res = size_matched_enrichment(
    is_member=gene_table.is_CGC,
    has_event=gene_table.has_m6A_event,
    covariates={"log_utr": np.log10(gene_table.total_utr_length),
                "log_mut": np.log10(1 + gene_table.n_mutations),
                "n_tx": gene_table.n_transcripts},
)
res["uncorrected"]["odds_ratio"], res["corrected"]["logistic"]["odds_ratio"], res["verdict"]
```

Reproduction check against the study's own gene-level table: uncorrected odds ratio
**1.96** and logistic **1.12** (95% CI 0.93-1.35) for Cancer Gene Census membership,
matching the published values, with identical 2x2 counts and the same verdict
(`attenuated`). Only the logistic adjustment is offered, so one correction of the confounder is
reported rather than several for the reader to choose between.
