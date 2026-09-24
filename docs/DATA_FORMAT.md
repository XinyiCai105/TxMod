# Data formats

## `record_key`

A stable identifier for one variant × transcript pair:

```
<transcript_id>|<gene_name>|<chrom>:<pos>|REF=<ref>|ALT=<alt>|strand=<±>|offset=<n>|len=<n>
```

`offset` is 1-based within the **spliced** 3′UTR in transcript orientation; `len`
is the spliced 3′UTR length. On the minus strand, offset 1 is the *highest* genomic
3′UTR coordinate, and the REF/ALT alleles in the key are genomic (the
transcript-orientation bases are the reverse complement, available as
`ref_base_transcript` / `alt_base_transcript`).

## `prepare --out` (pair table)

| Column | Meaning |
|---|---|
| `record_key` | identifier above |
| `transcript_id`, `gene_name` | from the GTF |
| `chrom`, `pos` | genomic, 1-based |
| `strand` | transcript strand |
| `ref_allele`, `alt_allele` | genomic alleles |
| `offset_1based` | position within the spliced 3′UTR |
| `utr_len` | spliced 3′UTR length |

## `prepare --fasta-out`

Two records per pair, headers `<record_key>|REF` and `<record_key>|MUT`, sequences
in transcript orientation with `U`→`T`. This is the input format for the m6A
prediction step when run standalone.

## `run --out` (result table)

The pair-table columns plus:

| Column | Meaning |
|---|---|
| `ref_prob_at_mut`, `mut_prob_at_mut` | iM6A probability at `offset-1` of the track |
| `delta_prob` | `mut − ref`; negative = predicted m6A loss |
| `m6a_class` | `m6A_gain` / `m6A_loss` / `neutral` |
| `motif_disrupted`, `motif_created` | RRACH change at the mutated position (motif is a parameter; RRACH by default) |
| `n_rbp_alterations` | high-confidence binding-alteration events |

## `run --events-out` (per-event RBP table)

One row per (pair × RBP motif) passing threshold: `record_key`, `gene_name`,
`transcript_id`, `rbp`, `motif_id`, `ref_norm`, `mut_norm`, `delta_norm`,
`direction`.

## PWM table (`panel --out`, `--pwm` input)

Long format, one row per motif position:

| Column | Meaning |
|---|---|
| `RBP` | RBP gene symbol |
| `motif_id` | motif identifier |
| `pos` | 1-based position within the motif |
| `A`, `C`, `G`, `U` | base probabilities (rows renormalised on load) |

## Thresholds

Defaults, all overridable:

- binding alteration: `max(ref_norm, mut_norm) ≥ 0.80` and `|delta_norm| ≥ 0.10`
- m6A high-confidence event (analysis convention, not enforced by the API):
  `|delta_prob| ≥ 0.10`
