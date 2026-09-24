"""
txmod.rbp
==============

RBP binding-site scanning with position weight matrices, and mutation-induced
binding alteration.

Scoring
-------
A PWM is converted to log2-odds against a uniform background (0.25 per base).
A sequence is scanned with every window; the best window's score is normalised
to ``[0, 1]`` against the motif's own achievable range::

    norm = (best - min_possible) / (max_possible - min_possible)

REF and MUT sequences are scored independently and compared. A high-confidence
binding alteration requires the site to be a credible motif match in at least one
allele and the change to be non-trivial::

    max(ref_norm, mut_norm) >= MIN_MAX_NORM      (default 0.80)
    abs(mut_norm - ref_norm) >= MIN_ABS_DELTA    (default 0.10)

Panel definition
----------------
:func:`load_panel` reads a PWM table; :func:`systematic_panel_from_cisbp` builds
one from a CISBP-RNA download restricted to a set of RBPs (e.g. those annotated
to mRNA-stability Gene Ontology terms), which is how the shipped panel is defined
rather than by hand-picking.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

#: Minimum ``max(ref_norm, mut_norm)`` for a credible motif match.
MIN_MAX_NORM = 0.80
#: Minimum ``abs(delta_norm)`` for a reported binding alteration.
MIN_ABS_DELTA = 0.10

_BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3}


@dataclass
class Motif:
    """A log2-odds scoring matrix for one RBP motif."""

    rbp: str
    motif_id: str
    log_matrix: np.ndarray  # (L, 4)
    min_possible: float
    max_possible: float

    @property
    def length(self) -> int:
        return int(self.log_matrix.shape[0])


def motif_from_probabilities(rbp: str, motif_id: str, probs: np.ndarray) -> Motif:
    """Build a :class:`Motif` from an ``(L, 4)`` probability matrix (A, C, G, U)."""
    p = np.asarray(probs, dtype=float)
    if p.ndim != 2 or p.shape[1] != 4:
        raise ValueError(f"{motif_id}: expected an (L, 4) matrix, got {p.shape}")
    p = np.clip(p, 1e-6, None)
    p = p / p.sum(axis=1, keepdims=True)
    log_mat = np.log2(p / 0.25)
    return Motif(
        rbp=rbp,
        motif_id=motif_id,
        log_matrix=log_mat,
        min_possible=float(log_mat.min(axis=1).sum()),
        max_possible=float(log_mat.max(axis=1).sum()),
    )


def encode(seq: str) -> np.ndarray:
    """Encode a sequence to base indices; unknown bases become ``-1``."""
    s = seq.upper()
    out = np.full(len(s), -1, dtype=np.int8)
    for b, i in _BASE_IDX.items():
        if not b:
            continue
        arr = np.frombuffer(s.encode(), dtype=np.uint8)
        out[arr == ord(b)] = i
    return out


def best_normalised_score(seq_or_code, motif: Motif) -> float:
    """Best normalised motif score over all windows; NaN when unscorable."""
    code = seq_or_code if isinstance(seq_or_code, np.ndarray) else encode(seq_or_code)
    L = motif.length
    if len(code) < L:
        return float("nan")
    windows = np.lib.stride_tricks.sliding_window_view(code, L)
    valid = (windows >= 0).all(axis=1)
    if not valid.any():
        return float("nan")
    w = windows[valid]
    scores = motif.log_matrix[np.arange(L)[None, :], w].sum(axis=1)
    denom = motif.max_possible - motif.min_possible
    if denom <= 0:
        return float("nan")
    return float((scores.max() - motif.min_possible) / denom)


@dataclass
class BindingAlteration:
    """Predicted change in binding potential for one RBP motif."""

    rbp: str
    motif_id: str
    ref_norm: float
    mut_norm: float

    @property
    def delta_norm(self) -> float:
        return self.mut_norm - self.ref_norm

    @property
    def max_norm(self) -> float:
        return max(self.ref_norm, self.mut_norm)

    @property
    def direction(self) -> str:
        d = self.delta_norm
        if d > 0:
            return "increased_binding_potential"
        if d < 0:
            return "decreased_binding_potential"
        return "unchanged"

    def is_high_confidence(
        self, min_max_norm: float = MIN_MAX_NORM, min_abs_delta: float = MIN_ABS_DELTA
    ) -> bool:
        if np.isnan(self.ref_norm) or np.isnan(self.mut_norm):
            return False
        return self.max_norm >= min_max_norm and abs(self.delta_norm) >= min_abs_delta


def scan_pair(
    ref_seq: str,
    mut_seq: str,
    motifs: Iterable[Motif],
    high_confidence_only: bool = True,
    min_max_norm: float = MIN_MAX_NORM,
    min_abs_delta: float = MIN_ABS_DELTA,
) -> List[BindingAlteration]:
    """Score a REF/MUT sequence pair against a motif panel."""
    cr = encode(ref_seq)
    cm = encode(mut_seq)
    out: List[BindingAlteration] = []
    for m in motifs:
        ba = BindingAlteration(
            rbp=m.rbp,
            motif_id=m.motif_id,
            ref_norm=best_normalised_score(cr, m),
            mut_norm=best_normalised_score(cm, m),
        )
        if high_confidence_only and not ba.is_high_confidence(min_max_norm, min_abs_delta):
            continue
        out.append(ba)
    return out


def load_panel(pwm_table: str) -> List[Motif]:
    """Load motifs from a long-format PWM table.

    Expected columns: ``RBP``, ``motif_id``, ``pos``, ``A``, ``C``, ``G``, ``U``.
    """
    import pandas as pd

    df = pd.read_csv(pwm_table, sep="\t")
    required = {"RBP", "motif_id", "pos", "A", "C", "G", "U"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{pwm_table}: missing column(s) {sorted(missing)}")
    motifs: List[Motif] = []
    for (rbp, mid), g in df.groupby(["RBP", "motif_id"], sort=True):
        g = g.sort_values("pos")
        motifs.append(motif_from_probabilities(str(rbp), str(mid),
                                               g[["A", "C", "G", "U"]].to_numpy()))
    return motifs


def systematic_panel_from_cisbp(
    rbp_information: str,
    pwm_dir: str,
    keep_rbps: Optional[Iterable[str]] = None,
) -> Tuple[List[Motif], Dict[str, List[str]]]:
    """Build a motif panel from a CISBP-RNA download.

    Parameters
    ----------
    rbp_information
        ``RBP_Information_all_motifs.txt`` from the CISBP-RNA bulk download.
    pwm_dir
        Directory of ``<motif_id>.txt`` PWM files (``pwms_all_motifs/``).
    keep_rbps
        Restrict to these RBP names — e.g. genes annotated to mRNA-stability GO
        terms — which is how a systematic (non hand-picked) panel is defined.

    Returns
    -------
    ``(motifs, skipped)`` where ``skipped`` maps a reason to the motif IDs
    dropped, so an empty or missing PWM file is reported rather than hidden.
    """
    import os

    import pandas as pd

    keep = {str(x) for x in keep_rbps} if keep_rbps is not None else None
    info = pd.read_csv(rbp_information, sep="\t", dtype=str)
    if "RBP_Name" not in info.columns or "Motif_ID" not in info.columns:
        raise ValueError(f"{rbp_information}: expected RBP_Name and Motif_ID columns")

    motifs: List[Motif] = []
    skipped: Dict[str, List[str]] = {"missing_file": [], "empty_matrix": [], "unreadable": []}
    seen: set[tuple[str, str]] = set()
    for _, row in info.iterrows():
        rbp = row.get("RBP_Name")
        mid = row.get("Motif_ID")
        if not rbp or not mid or mid == ".":
            continue
        if keep is not None and rbp not in keep:
            continue
        if (rbp, mid) in seen:
            continue
        seen.add((rbp, mid))
        path = os.path.join(pwm_dir, f"{mid}.txt")
        if not os.path.exists(path):
            skipped["missing_file"].append(mid)
            continue
        try:
            mat = pd.read_csv(path, sep="\t")
        except Exception:
            skipped["unreadable"].append(mid)
            continue
        if mat.empty or not {"A", "C", "G", "U"} <= set(mat.columns):
            skipped["empty_matrix"].append(mid)
            continue
        if "Pos" in mat.columns:
            mat = mat.sort_values("Pos")
        motifs.append(motif_from_probabilities(str(rbp), str(mid),
                                               mat[["A", "C", "G", "U"]].to_numpy()))
    return motifs, skipped


#: Minimum number of sites below which a fold-enrichment is not interpretable.
MIN_SITES_FOR_ENRICHMENT = 20


def enrichment_with_power(
    n_hits: int,
    n_sites: int,
    background_rate: float,
    min_sites: int = MIN_SITES_FOR_ENRICHMENT,
) -> dict:
    """Fold enrichment over a background rate, with an explicit power flag.

    A fold-enrichment computed on a handful of sites is arithmetically valid and
    scientifically meaningless: a 2-site observation can read as "2.1x enriched"
    while the Fisher test returns p = 0.23. This helper always returns ``n`` and a
    ``underpowered`` flag next to the ratio so a caller cannot quote the ratio
    without also seeing what it rests on.

    Parameters
    ----------
    n_hits
        Sites overlapping the reference feature (e.g. a real m6A site or a CLIP peak).
    n_sites
        Sites tested.
    background_rate
        Expected overlap fraction for an unselected site from the same universe.
    min_sites
        Below this, ``underpowered`` is True and ``headline_safe`` is False.

    Returns
    -------
    dict with ``n_sites``, ``n_hits``, ``observed_rate``, ``expected_rate``,
    ``fold_enrichment``, ``fisher_p``, ``underpowered`` and ``headline_safe``.

    Examples
    --------
    >>> r = enrichment_with_power(2, 2, 0.476)
    >>> r["fold_enrichment"] > 2 and r["underpowered"]
    True
    >>> r["headline_safe"]
    False
    """
    from scipy.stats import fisher_exact

    if n_sites <= 0:
        raise ValueError("n_sites must be positive")
    if not 0.0 < background_rate < 1.0:
        raise ValueError("background_rate must lie strictly between 0 and 1")
    obs = n_hits / n_sites
    exp_hits = int(round(background_rate * n_sites))
    odds, p = fisher_exact(
        [[n_hits, n_sites - n_hits], [exp_hits, n_sites - exp_hits]],
        alternative="greater",
    )
    underpowered = n_sites < min_sites
    return {
        "n_sites": int(n_sites),
        "n_hits": int(n_hits),
        "observed_rate": float(obs),
        "expected_rate": float(background_rate),
        "fold_enrichment": float(obs / background_rate),
        "fisher_p": float(p),
        "underpowered": bool(underpowered),
        "headline_safe": bool(not underpowered and p < 0.05),
    }
