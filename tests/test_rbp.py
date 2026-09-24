"""Tests for PWM scoring, event thresholds and enrichment power."""
import numpy as np
import pytest

from txmod.rbp import (BindingAlteration, best_normalised_score,
                       motif_from_probabilities, scan_pair)


def _perfect_A_motif(n=4):
    """A motif that only matches poly-A."""
    p = np.full((n, 4), 0.01)
    p[:, 0] = 0.97
    return motif_from_probabilities("TESTRBP", "M1", p)


def test_normalised_score_bounds_and_best_match():
    m = _perfect_A_motif()
    assert best_normalised_score("AAAA", m) == pytest.approx(1.0, abs=1e-6)
    # a poor match scores well below 1 but stays in [0,1]
    s = best_normalised_score("CCCC", m)
    assert 0.0 <= s < 0.5


def test_score_picks_best_window():
    m = _perfect_A_motif()
    assert best_normalised_score("CCCCAAAACCCC", m) == pytest.approx(1.0, abs=1e-6)


def test_too_short_sequence_is_nan():
    assert np.isnan(best_normalised_score("AA", _perfect_A_motif()))


def test_unknown_bases_do_not_crash():
    assert np.isnan(best_normalised_score("NNNN", _perfect_A_motif()))


def test_high_confidence_thresholds():
    ok = BindingAlteration("R", "M", ref_norm=0.90, mut_norm=0.70)
    assert ok.is_high_confidence()                # max 0.90, |delta| 0.20
    weak = BindingAlteration("R", "M", ref_norm=0.50, mut_norm=0.30)
    assert not weak.is_high_confidence()          # max below 0.80
    small = BindingAlteration("R", "M", ref_norm=0.90, mut_norm=0.88)
    assert not small.is_high_confidence()         # |delta| below 0.10
    assert ok.direction == "decreased_binding_potential"


def test_scan_pair_detects_loss_of_binding():
    m = _perfect_A_motif()
    events = scan_pair("CCAAAACC", "CCAACACC", [m])
    assert len(events) == 1 and events[0].delta_norm < 0







def test_enrichment_with_power_flags_small_n():
    """The n=2 / 2.1x case that must never be quoted as a headline."""
    from txmod.rbp import enrichment_with_power
    r = enrichment_with_power(n_hits=2, n_sites=2, background_rate=0.476)
    assert r["fold_enrichment"] > 2.0          # the ratio looks impressive
    assert r["fisher_p"] > 0.05                # but it is not significant
    assert r["underpowered"] is True
    assert r["headline_safe"] is False


def test_enrichment_with_power_accepts_well_powered_result():
    """The n=153 / 1.88x case that is quotable."""
    from txmod.rbp import enrichment_with_power
    r = enrichment_with_power(n_hits=137, n_sites=153, background_rate=0.476)
    assert r["fold_enrichment"] == pytest.approx(1.88, abs=0.02)
    assert r["underpowered"] is False
    assert r["headline_safe"] is True


def test_enrichment_rejects_impossible_inputs():
    from txmod.rbp import enrichment_with_power
    with pytest.raises(ValueError):
        enrichment_with_power(1, 0, 0.5)
    with pytest.raises(ValueError):
        enrichment_with_power(1, 10, 0.0)
    with pytest.raises(ValueError):
        enrichment_with_power(1, 10, 1.0)


