"""
txmod.background
=====================

Matched-background comparisons.

Three analyses in the TxMod study reached a different conclusion once the
comparison was made against a *matched* background rather than a raw one, and two
of them reversed. This module provides those comparisons as reusable functions,
because the failure mode they guard against is generic: an apparent signal that is
really a property of what the cases are made of.

The three settings
------------------
1. **Sequence composition** (:func:`trinucleotide_rate_model`,
   :func:`depletion_observed_expected`). Somatic mutation rate depends strongly on
   trinucleotide context. Functional elements defined by a sequence motif therefore
   sit in a non-random composition, and comparing raw mutation density inside and
   outside them confounds selection with mutability. In this dataset the RAC
   contexts that dominate predicted m6A sites mutate about an order of magnitude
   less often than CpG contexts, which fully accounted for an apparent depletion.

2. **Gene size** (:func:`size_matched_enrichment`). Longer, more isoform-rich,
   more heavily mutated genes have more opportunity to acquire any particular
   variant class. Testing membership of a curated gene set against an unmatched
   background of all genes therefore manufactures enrichment.

Each function returns the uncorrected and corrected estimate together, so a caller
cannot report one without seeing the other.

Optional dependency
-------------------
:func:`size_matched_enrichment` uses ``statsmodels`` for the logistic model
(``pip install txmod[stats]``). Without it, the corrected estimate is reported as
unavailable and only the uncorrected Fisher test is returned. The absence is
signalled in the returned dict rather than silently substituted.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def canonical_context(trinuc: str) -> Optional[str]:
    """Fold a trinucleotide onto its pyrimidine-referenced form.

    A context and its reverse complement describe the same physical position seen
    from either strand, so they are collapsed onto the form carrying a pyrimidine
    (C or T) at the centre. This halves 64 trinucleotides to 32 categories.

    Returns the folded trinucleotide, or ``None`` if any base is not canonical.

    Examples
    --------
    >>> canonical_context("ACA")
    'ACA'
    >>> canonical_context("TGT")     # the same context on the minus strand
    'ACA'
    """
    if len(trinuc) != 3 or any(b not in "ACGT" for b in trinuc):
        return None
    return trinuc if trinuc[1] in "CT" else _revcomp(trinuc)


@dataclass
class RateModel:
    """A trinucleotide somatic mutation rate model.

    The quantity of interest is whether a position is mutated, not what it changes
    to, so rates are per position and carry no substitution-type dimension.

    Attributes
    ----------
    rates
        ``trinucleotide -> per-position mutation rate``.
    context_counts
        ``trinucleotide -> number of positions of that context`` in the universe.
    mutation_counts
        ``trinucleotide -> number of mutated positions of that context``.
    """

    rates: Dict[str, float] = field(default_factory=dict)
    context_counts: Dict[str, int] = field(default_factory=dict)
    mutation_counts: Dict[str, int] = field(default_factory=dict)

    def total_rate(self, trinuc: str) -> float:
        """Mutation rate for one trinucleotide, folded onto its canonical form."""
        key = canonical_context(trinuc)
        return float(self.rates.get(key, 0.0)) if key else 0.0

    def expected(self, trinuc: str) -> float:
        """Expected mutations at a single position of this context."""
        return self.total_rate(trinuc)

    def summary(self, top: int = 8) -> List[Tuple[str, float]]:
        """Trinucleotides ranked by mutation rate, highest first."""
        tris = sorted(self.context_counts, key=self.total_rate, reverse=True)
        return [(t, self.total_rate(t)) for t in tris[:top]]


def trinucleotide_rate_model(
    sequences: Mapping[str, str],
    mutations: Iterable[Tuple[str, int]],
    margin: int = 1,
) -> RateModel:
    """Fit a trinucleotide mutation-rate model over a set of sequences.

    Parameters
    ----------
    sequences
        ``id -> sequence`` for the universe of interrogable positions (e.g. every
        3'UTR in the dataset).
    mutations
        Iterable of ``(sequence_id, position0)``; ``position0`` is 0-based within
        that sequence. A position mutated in several samples counts once, since
        the model estimates whether a position is mutated at all.
    margin
        Positions closer than this to either end are excluded, since their
        trinucleotide context is incomplete.

    Returns
    -------
    RateModel

    Notes
    -----
    The denominator is every interior position of every supplied sequence, so the
    model is calibrated on the same universe the comparison will be drawn from.
    Contexts are folded onto their pyrimidine-referenced form, giving 32
    categories, so a rate does not depend on which strand a sequence was taken
    from.
    """
    context_counts: Dict[str, int] = defaultdict(int)
    for seq in sequences.values():
        s = seq.upper()
        for i in range(margin, len(s) - margin):
            key = canonical_context(s[i - 1 : i + 2])
            if key is not None:
                context_counts[key] += 1

    # A position mutated in several samples is one mutated position, so deduplicate
    # before counting.
    seen = set()
    mut_counts: Dict[str, int] = defaultdict(int)
    for sid, pos0 in mutations:
        seq = sequences.get(sid)
        if seq is None or pos0 < margin or pos0 >= len(seq) - margin:
            continue
        if (sid, pos0) in seen:
            continue
        seen.add((sid, pos0))
        ctx = canonical_context(seq[pos0 - 1 : pos0 + 2].upper())
        if ctx is not None:
            mut_counts[ctx] += 1

    rates = {
        tri: mut_counts.get(tri, 0) / context_counts[tri]
        for tri in context_counts
        if context_counts[tri] > 0
    }
    return RateModel(rates=dict(rates),
                     context_counts=dict(context_counts),
                     mutation_counts=dict(mut_counts))


def depletion_observed_expected(
    sites: Iterable[Tuple[str, int]],
    sequences: Mapping[str, str],
    mutation_positions: Mapping[str, Sequence[int]],
    model: RateModel,
    centre_halfwidth: int = 25,
    flank_inner: int = 25,
    flank_outer: int = 100,
) -> dict:
    """Mutation depletion around a set of sites, raw and composition-corrected.

    For each site, observed and expected mutation counts are accumulated over a
    central window and a flanking window; expected counts come from ``model``, so
    the corrected statistic asks whether the centre carries fewer mutations *than
    its own sequence composition predicts*.

    Parameters
    ----------
    sites
        ``(sequence_id, position0)`` of each site of interest.
    sequences
        ``id -> sequence``, as used to fit ``model``.
    mutation_positions
        ``sequence_id -> sorted 0-based mutated positions``.
    model
        A fitted :class:`RateModel`.
    centre_halfwidth
        Central window is ``pos +/- centre_halfwidth``.
    flank_inner, flank_outer
        Flanking window is ``flank_inner < |offset| <= flank_outer``.

    Returns
    -------
    dict with ``n_sites``, ``raw_centre_flank_ratio``, ``centre_oe``, ``flank_oe``,
    ``corrected_ratio_of_ratios``, ``ci95`` (Poisson error propagation) and
    ``p_centre_depletion`` (one-sided Poisson test). A corrected value of 1.0 means
    no depletion beyond what composition predicts.

    Notes
    -----
    The raw ratio and the corrected ratio are both returned deliberately: the gap
    between them *is* the composition artefact, and reporting the corrected value
    alone hides how large that artefact was.
    """
    from scipy import stats

    obs_c = obs_f = 0.0
    exp_c = exp_f = 0.0
    n = 0
    for sid, pos in sites:
        seq = sequences.get(sid)
        if seq is None:
            continue
        muts = set(mutation_positions.get(sid, ()))
        lo, hi = 1, len(seq) - 1
        n += 1
        for off in range(-flank_outer, flank_outer + 1):
            i = pos + off
            if not (lo <= i < hi):
                continue
            d = abs(off)
            if d <= centre_halfwidth:
                bucket = "c"
            elif flank_inner < d <= flank_outer:
                bucket = "f"
            else:
                continue
            tri = seq[i - 1 : i + 2].upper()
            if len(tri) != 3 or any(b not in "ACGT" for b in tri):
                continue
            key = tri if tri[1] in "CT" else _revcomp(tri)
            e = model.expected(key)
            o = 1.0 if i in muts else 0.0
            if bucket == "c":
                obs_c += o
                exp_c += e
            else:
                obs_f += o
                exp_f += e

    centre_oe = obs_c / exp_c if exp_c > 0 else float("nan")
    flank_oe = obs_f / exp_f if exp_f > 0 else float("nan")
    ror = centre_oe / flank_oe if flank_oe else float("nan")

    # Poisson error propagation on the log of a ratio of ratios
    if obs_c > 0 and obs_f > 0:
        se_log = np.sqrt(1.0 / obs_c + 1.0 / obs_f)
        ci = (float(ror * np.exp(-1.96 * se_log)), float(ror * np.exp(1.96 * se_log)))
    else:
        ci = (float("nan"), float("nan"))

    if exp_c > 0 and flank_oe == flank_oe and flank_oe > 0:
        mu = exp_c * flank_oe          # centre count expected if centre behaved as flank
        p = float(stats.poisson.cdf(obs_c, mu))
    else:
        p = float("nan")

    # raw (uncorrected) density ratio: mutations per interrogated position
    n_centre_pos = 2 * centre_halfwidth + 1
    n_flank_pos = 2 * (flank_outer - flank_inner)
    density_c = obs_c / n_centre_pos
    density_f = obs_f / n_flank_pos
    raw = density_c / density_f if density_f > 0 else float("nan")

    return {
        "n_sites": n,
        "observed_centre": obs_c,
        "observed_flank": obs_f,
        "expected_centre": exp_c,
        "expected_flank": exp_f,
        "raw_centre_flank_ratio": float(raw),
        "centre_oe": float(centre_oe),
        "flank_oe": float(flank_oe),
        "corrected_ratio_of_ratios": float(ror),
        "ci95": ci,
        "p_centre_depletion": p,
        "no_depletion": bool(ci[0] <= 1.0 <= ci[1]) if ci[0] == ci[0] else None,
    }


# --------------------------------------------------------------------------- #
# 2. Size-matched gene-set enrichment
# --------------------------------------------------------------------------- #

def size_matched_enrichment(
    is_member: Sequence[bool],
    has_event: Sequence[bool],
    covariates: Mapping[str, Sequence[float]],
) -> dict:
    """Gene-set enrichment, uncorrected and corrected for size confounders.

    Asks whether genes carrying some event are enriched for membership of a curated
    set (e.g. a cancer gene census), first against the raw background and then
    controlling for covariates that drive *opportunity* rather than biology --
    typically total regulatory-region length, variant count and isoform number.

    Parameters
    ----------
    is_member
        Per-gene membership of the curated set.
    has_event
        Per-gene event status (the exposure being tested).
    covariates
        ``name -> per-gene value``. Values spanning orders of magnitude (lengths,
        counts) should be supplied already log-transformed by the caller, or pass
        them raw and read ``covariate_transform`` in the result.

    Returns
    -------
    dict with an ``uncorrected`` block (Fisher) and a ``corrected`` block holding
    the ``logistic`` estimate with odds ratio, confidence interval, p value and
    covariate list, plus ``verdict`` in {``survives``, ``attenuated``,
    ``disappears``, ``underpowered``}.

    Notes
    -----
    Both estimates are returned together so neither can be quoted in isolation: the
    uncorrected odds ratio is what the correction exists to explain, and reporting
    only the corrected one hides the size of the artefact.
    """
    from scipy.stats import fisher_exact

    member = np.asarray(is_member, dtype=bool)
    event = np.asarray(has_event, dtype=bool)
    if member.shape != event.shape:
        raise ValueError("is_member and has_event must have the same length")
    n = member.size
    X = {k: np.asarray(v, dtype=float) for k, v in covariates.items()}
    for k, v in X.items():
        if v.size != n:
            raise ValueError(f"covariate {k!r} has length {v.size}, expected {n}")

    a = int((event & member).sum())
    b = int((event & ~member).sum())
    c = int((~event & member).sum())
    d = int((~event & ~member).sum())
    or_raw, p_raw = fisher_exact([[a, b], [c, d]])
    uncorrected = {
        "odds_ratio": float(or_raw),
        "p": float(p_raw),
        "counts": {"event_member": a, "event_nonmember": b,
                   "noevent_member": c, "noevent_nonmember": d},
        "member_rate_event": a / max(a + b, 1),
        "member_rate_noevent": c / max(c + d, 1),
    }

    corrected: Dict[str, dict] = {}

    # (a) logistic regression
    try:
        import statsmodels.api as sm

        design = np.column_stack([event.astype(float)] + [X[k] for k in sorted(X)])
        design = sm.add_constant(design, has_constant="add")
        fit = sm.Logit(member.astype(float), design).fit(disp=0)
        coef = float(fit.params[1])
        lo, hi = (float(x) for x in fit.conf_int()[1])
        corrected["logistic"] = {
            "odds_ratio": float(np.exp(coef)),
            "ci95": (float(np.exp(lo)), float(np.exp(hi))),
            "p": float(fit.pvalues[1]),
            "coefficient": coef,
            "covariates": sorted(X),
            "converged": bool(fit.mle_retvals.get("converged", True)),
        }
    except ImportError:
        corrected["logistic"] = {"error": "statsmodels not installed"}

    est = corrected.get("logistic", {})
    or_adj = est.get("odds_ratio")
    p_adj = est.get("p")
    if a + c < 10 or or_adj is None:
        verdict = "underpowered"
    elif or_adj < 1.2 and (p_adj is None or p_adj > 0.05):
        verdict = "disappears"
    elif or_adj < 0.75 * float(or_raw):
        verdict = "attenuated"
    else:
        verdict = "survives"

    return {"uncorrected": uncorrected, "corrected": corrected, "verdict": verdict,
            "n_genes": n, "n_members": int(member.sum()), "n_events": int(event.sum())}
