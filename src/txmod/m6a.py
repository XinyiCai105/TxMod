"""
txmod.m6a
==============

m6A methylation prediction and mutation-induced change in predicted
probability (Delta-prob).

The predictor is the iM6A ensemble (Keras/TensorFlow ``.h5`` models). TensorFlow
is imported lazily so the rest of the package -- annotation, sequence building,
motif scanning -- works without a deep-learning stack installed.

Registration
------------
The iM6A output track is offset by one nucleotide relative to the input
sequence: probability index ``i`` scores the methylation of the adenosine at
sequence index ``i + 1`` (0-based). :data:`REGISTRATION_OFFSET` records this and
:func:`methylated_base_index` applies it.

The offset is a property of the model rather than of any particular result
table, so it should be measured, not assumed. :func:`calibrate_shift` recovers
it from real sequences by asking which shift places high-probability calls on
the adenosine of an RRACH, and refuses to guess when the evidence is ambiguous.

Reading the probability change
------------------------
A substitution can alter methylation of any adenosine whose motif context it
touches, not only the base it replaces. An RRACH is five nucleotides wide, so
the methylated adenosine of an affected motif lies within +/- 2 nt of the
mutation. :func:`probability_change` scans that window and reports the adenosine whose
predicted probability changes most, with its offset relative to the mutation
(``scored_A_rel``, in -2..+2).

Reading one index instead covers a single motif position and silently drops
the other four. On the COSMIC 3'UTR set a single-index read gives 2,299 events at
a 4.19:1 loss:gain ratio; the window read gives 14,274 at 1.15:1, so that
imbalance is an artefact of the read position rather than a property of the
mutations. The single-index value is kept as ``delta_prob_shift0`` so the two
reads can be compared directly.

Motif
-----
``RRACH`` (``R=[AG] R=[AG] A C H=[ACT]``) is the motif used throughout: it is the
training definition of the ``humanRRACH10000`` ensemble. A different consensus can
be supplied through the ``motif=`` argument, but pentamers outside the training
definition are scored highly only once a mutation converts them to RRACH, which
produces gain-only calls.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

#: iM6A probability index ``i`` scores the adenosine at sequence index ``i + 1``.
REGISTRATION_OFFSET = 1

#: Model-matched consensus: the ``humanRRACH10000`` ensemble's own definition.
RRACH = re.compile(r"^[AG][AG]AC[ACT]$")

#: Half-width of the window scanned around a mutation, in nucleotides. An RRACH
#: is 5 nt wide, so a mutation anywhere in the motif is within 2 nt of its
#: methylated adenosine.
WINDOW = 2

_BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3, "U": 3}


def methylated_base_index(prob_index: int) -> int:
    """Sequence index of the adenosine scored by probability index ``prob_index``."""
    return prob_index + REGISTRATION_OFFSET


def one_hot_encode(seqs: Sequence[str]) -> np.ndarray:
    """One-hot encode equal-length sequences to ``(n, L, 4)``; unknown bases are zero."""
    n = len(seqs)
    L = len(seqs[0]) if n else 0
    x = np.zeros((n, L, 4), dtype=np.float32)
    for i, s in enumerate(seqs):
        for j, b in enumerate(s.upper()):
            k = _BASE_IDX.get(b)
            if k is not None:
                x[i, j, k] = 1.0
    return x


def load_models(model_dir: str, model_prefix: str = "humanRRACH10000", n_models: int = 5):
    """Load the iM6A ensemble. Imports TensorFlow lazily.

    Expects ``{model_dir}/{model_prefix}_c{1..n}.h5``.
    """
    import os

    from tensorflow.keras import backend as kb  # type: ignore
    from tensorflow.keras.models import load_model  # type: ignore

    def categorical_crossentropy_2d(y_true, y_pred):
        return -kb.mean(
            y_true[:, :, 0] * kb.log(y_pred[:, :, 0] + 1e-10)
            + y_true[:, :, 1] * kb.log(y_pred[:, :, 1] + 1e-10)
        )

    models = []
    for i in range(1, n_models + 1):
        path = os.path.join(model_dir, f"{model_prefix}_c{i}.h5")
        if not os.path.exists(path):
            raise FileNotFoundError(path)
        models.append(
            load_model(
                path,
                custom_objects={"categorical_crossentropy_2d": categorical_crossentropy_2d},
                compile=False,
            )
        )
    return models


def predict_track(
    seqs: Sequence[str],
    models,
    context: int = 10000,
    batch_size: int = 32,
) -> List[np.ndarray]:
    """Per-nucleotide m6A probability for each sequence.

    Sequences are N-padded by ``context // 2`` on both sides, predicted by the
    ensemble (mean), then cropped back so output length equals input length.
    """
    if not seqs:
        return []
    pad = context // 2
    max_len = max(len(s) for s in seqs)
    padded = [
        ("N" * pad) + s.upper().replace("U", "T") + ("N" * (max_len - len(s))) + ("N" * pad)
        for s in seqs
    ]
    x = one_hot_encode(padded)
    preds = [m.predict(x, batch_size=batch_size, verbose=0) for m in models]
    y = np.mean(preds, axis=0)
    prob = y[:, :, 1]
    if prob.shape[1] >= pad + max_len:
        prob = prob[:, pad : pad + max_len]
    return [prob[i, : len(s)] for i, s in enumerate(seqs)]


@dataclass
class M6AResult:
    """Predicted m6A probability change for one variant-transcript pair."""

    ref_prob_at_mut: float
    mut_prob_at_mut: float
    delta_prob: float
    max_ref_prob: float
    max_mut_prob: float
    #: Offset of the scored adenosine relative to the mutation, in -2..+2.
    #: ``None`` when the window contains no adenosine.
    scored_A_rel: Optional[int] = None
    #: 1-based sequence position of the scored adenosine.
    scored_A_offset_1based: Optional[int] = None
    #: Number of adenosines found in the window.
    n_A_in_window: int = 0
    #: Delta-prob read at the mutation index alone (single-index read).
    #: Kept so window and single-index reads can be compared on the same run.
    delta_prob_shift0: float = float("nan")
    #: Registration offset actually used.
    track_shift: int = REGISTRATION_OFFSET

    @property
    def m6a_class(self) -> str:
        if self.delta_prob > 0:
            return "m6A_gain"
        if self.delta_prob < 0:
            return "m6A_loss"
        return "neutral"


def probability_change(
    ref_track: np.ndarray,
    mut_track: np.ndarray,
    offset_1based: int,
    ref_seq: Optional[str] = None,
    mut_seq: Optional[str] = None,
    window: int = WINDOW,
    shift: int = REGISTRATION_OFFSET,
) -> M6AResult:
    """Largest Delta-prob among adenosines within ``window`` nt of the mutation.

    Every adenosine whose motif context the substitution can touch is scored, and
    the one whose probability moves most is reported. A base counts as an
    adenosine if it is ``A`` in either the reference or the mutant sequence, so
    substitutions that create a methylation site are captured alongside those
    that destroy one.

    ``ref_seq``/``mut_seq`` are the sequences the tracks were predicted from. If
    both are omitted the window cannot be restricted to adenosines and every
    position in it is considered, which is looser but never misses a site.

    Raises
    ------
    IndexError
        If ``offset_1based`` lies outside the tracks.
    """
    mut0 = offset_1based - 1
    n = min(len(ref_track), len(mut_track))
    if mut0 < 0 or mut0 >= n:
        raise IndexError(f"offset {offset_1based} outside track of length {n}")

    def prob_index(seq_pos: int) -> int:
        """Track index scoring the base at 0-based ``seq_pos``."""
        return seq_pos - shift

    best: Optional[Tuple[float, int, int, float, float]] = None
    n_a = 0
    for rel in range(-window, window + 1):
        pos = mut0 + rel                      # candidate methylated A, 0-based
        idx = prob_index(pos)
        if idx < 0 or idx >= n:
            continue
        if ref_seq is not None or mut_seq is not None:
            rb = ref_seq[pos].upper() if ref_seq and pos < len(ref_seq) else ""
            mb = mut_seq[pos].upper() if mut_seq and pos < len(mut_seq) else ""
            if "A" not in (rb, mb):
                continue
        n_a += 1
        rp, mp = float(ref_track[idx]), float(mut_track[idx])
        d = mp - rp
        if best is None or abs(d) > abs(best[0]):
            best = (d, rel, pos + 1, rp, mp)

    i0 = prob_index(mut0)
    d0 = (
        float(mut_track[i0]) - float(ref_track[i0])
        if 0 <= i0 < n
        else float("nan")
    )

    if best is None:
        return M6AResult(
            ref_prob_at_mut=float("nan"),
            mut_prob_at_mut=float("nan"),
            delta_prob=float("nan"),
            max_ref_prob=float(np.max(ref_track)) if len(ref_track) else float("nan"),
            max_mut_prob=float(np.max(mut_track)) if len(mut_track) else float("nan"),
            n_A_in_window=0,
            delta_prob_shift0=d0,
            track_shift=shift,
        )

    d, rel, pos1, rp, mp = best
    return M6AResult(
        ref_prob_at_mut=rp,
        mut_prob_at_mut=mp,
        delta_prob=d,
        max_ref_prob=float(np.max(ref_track)) if len(ref_track) else float("nan"),
        max_mut_prob=float(np.max(mut_track)) if len(mut_track) else float("nan"),
        scored_A_rel=rel,
        scored_A_offset_1based=pos1,
        n_A_in_window=n_a,
        delta_prob_shift0=d0,
        track_shift=shift,
    )


def calibrate_shift(
    tracks: Sequence[np.ndarray],
    seqs: Sequence[str],
    candidates: Sequence[int] = (0, 1, -1),
    min_prob: float = 0.2,
    min_calls: int = 8,
    motif: "re.Pattern[str]" = RRACH,
) -> int:
    """Measure the track-to-sequence offset from real sequences.

    For each candidate shift, counts how often a high-probability call lands on
    the adenosine of a ``motif`` pentamer. The winning shift must both take a
    strict majority of the calls and be supported by at least ``min_calls``
    high-probability positions.

    Raises
    ------
    RuntimeError
        If too few high-probability positions are available, or no candidate
        wins a majority. Guessing an offset silently mis-reads every result, so
        an ambiguous calibration is an error rather than a default.
    """
    votes = {s: 0 for s in candidates}
    total = 0
    for track, seq in zip(tracks, seqs):
        s = seq.upper().replace("U", "T")
        hits = np.flatnonzero(np.asarray(track) >= min_prob)
        for idx in hits:
            total += 1
            for shift in candidates:
                pos = int(idx) + shift
                if pos < 2 or pos + 3 > len(s):
                    continue
                if s[pos] == "A" and motif.match(s[pos - 2 : pos + 3]):
                    votes[shift] += 1

    if total < min_calls:
        raise RuntimeError(
            f"only {total} positions with prob >= {min_prob}; need >= {min_calls} "
            "to calibrate. Supply more sequences or lower min_prob."
        )
    best = max(votes, key=lambda k: votes[k])
    if votes[best] * 2 <= sum(votes.values()):
        raise RuntimeError(
            f"calibration ambiguous (votes {votes} over {total} calls); "
            "refusing to guess the registration offset."
        )
    return best


def overlaps_motif(
    seq: str,
    offset_1based: int,
    motif: "re.Pattern[str]" = RRACH,
) -> bool:
    """True when the mutated base falls inside any ``motif`` 5-mer.

    Defaults to RRACH, the ensemble's training definition. A different
    for the broader consensus.
    """
    if not seq or offset_1based is None:
        return False
    s = seq.upper().replace("U", "T").replace("-", "")
    mut0 = int(offset_1based) - 1
    n = len(s)
    for start in range(mut0 - 4, mut0 + 1):
        end = start + 5
        if start < 0 or end > n:
            continue
        if motif.match(s[start:end]):
            return True
    return False


def motif_change(
    ref_seq: str,
    mut_seq: str,
    offset_1based: int,
    motif: "re.Pattern[str]" = RRACH,
) -> Dict[str, bool]:
    """Whether the variant disrupts or creates ``motif`` at the mutated position.

    Defaults to RRACH, the ensemble's training definition. Pass ``motif`` to
    test a different pattern. The returned keys do not name the motif, so they
    stay correct whichever pattern is used.
    """
    r = overlaps_motif(ref_seq, offset_1based, motif=motif)
    m = overlaps_motif(mut_seq, offset_1based, motif=motif)
    return {
        "ref_in_motif": r,
        "mut_in_motif": m,
        "motif_disrupted": bool(r and not m),
        "motif_created": bool(m and not r),
    }


def isoform_divergence(delta_probs: Sequence[float]) -> Optional[float]:
    """Spread of Delta-prob across isoforms of one variant (max - min).

    The core transcript-resolved statistic: a large value means the same genomic
    variant has isoform-dependent consequences.
    """
    vals = [d for d in delta_probs if d is not None and not np.isnan(d)]
    if len(vals) < 2:
        return None
    return float(max(vals) - min(vals))
