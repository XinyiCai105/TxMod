"""Tests for RRACH logic, +1 registration and Delta-prob extraction."""
import numpy as np
import pytest

from txmod.m6a import (REGISTRATION_OFFSET, isoform_divergence,
                            methylated_base_index, motif_change, overlaps_motif,
                            probability_change)


def test_registration_offset_is_one():
    """prob index i scores the adenosine at sequence index i+1."""
    assert REGISTRATION_OFFSET == 1
    assert methylated_base_index(10) == 11


def test_rrach_detection_positive_and_negative():
    # GGACT is a canonical RRACH; mutation at the central A (index 3, 1-based)
    assert overlaps_motif("GGACT", 3) is True
    # CCCCC contains no RRACH
    assert overlaps_motif("CCCCC", 3) is False


def test_rrach_window_spans_mutation():
    """A mutation anywhere within an RRACH 5-mer counts as overlapping."""
    seq = "TTGGACTTT"      # GGACT at 1-based 3..7
    for pos in range(3, 8):
        assert overlaps_motif(seq, pos) is True
    assert overlaps_motif(seq, 1) is False


def test_rrach_disruption_and_creation():
    ref, mut = "TTGGACTTT", "TTGGCCTTT"     # A->C destroys the RRACH
    d = motif_change(ref, mut, 5)
    assert d["motif_disrupted"] and not d["motif_created"]
    d2 = motif_change(mut, ref, 5)
    assert d2["motif_created"] and not d2["motif_disrupted"]


def test_probability_change_sign_and_class():
    ref = np.array([0.1, 0.8, 0.2])
    mut = np.array([0.1, 0.3, 0.2])
    r = probability_change(ref, mut, offset_1based=2)
    assert r.delta_prob == pytest.approx(-0.5)
    assert r.m6a_class == "m6A_loss"
    r2 = probability_change(mut, ref, offset_1based=2)
    assert r2.m6a_class == "m6A_gain"


def test_probability_change_out_of_range_raises():
    with pytest.raises(IndexError):
        probability_change(np.zeros(3), np.zeros(3), offset_1based=99)


def test_window_finds_site_the_single_index_read_misses():
    """The mutation sits on the D of an RRACH; the methylated A is 2 nt away.

    A single-index read at the mutation offset lands on a flat part of the
    track and reports nothing. The window read finds the destroyed site.
    """
    #        pos: 0123456789
    ref_seq = "TTGGACTTTT"      # GGACT at 2..6, methylated A at index 4
    mut_seq = "TTTGACTTTT"      # G->T at index 2 (the D position)
    ref = np.zeros(10)
    mut = np.zeros(10)
    # +1 registration: probability index i scores sequence index i + 1,
    # so the A at sequence index 4 is scored by track index 3.
    ref[3] = 0.9
    mut[3] = 0.1

    r = probability_change(ref, mut, offset_1based=3, ref_seq=ref_seq, mut_seq=mut_seq)
    assert r.delta_prob == pytest.approx(-0.8)
    assert r.scored_A_rel == 2                     # A is 2 nt 3' of the mutation
    assert r.scored_A_offset_1based == 5           # 1-based position of that A
    assert r.m6a_class == "m6A_loss"
    # the old single-index read saw nothing at this position
    assert r.delta_prob_shift0 == pytest.approx(0.0)


def test_window_scores_created_site():
    """A substitution that creates an adenosine is scored at that new A."""
    ref_seq = "TTGGCCTTTT"
    mut_seq = "TTGGACTTTT"      # C->A at index 4 creates GGACT
    ref = np.zeros(10)
    mut = np.zeros(10)
    mut[3] = 0.7                # track index 3 scores sequence index 4

    r = probability_change(ref, mut, offset_1based=5, ref_seq=ref_seq, mut_seq=mut_seq)
    assert r.delta_prob == pytest.approx(0.7)
    assert r.scored_A_rel == 0                     # the mutated base itself
    assert r.m6a_class == "m6A_gain"


def test_window_reports_no_adenosine():
    ref_seq = mut_seq = "TTTTTTTTTT"
    r = probability_change(np.zeros(10), np.zeros(10), offset_1based=5,
                     ref_seq=ref_seq, mut_seq=mut_seq)
    assert r.n_A_in_window == 0
    assert np.isnan(r.delta_prob)
    assert r.scored_A_rel is None


def test_calibrate_shift_recovers_plus_one():
    from txmod.m6a import calibrate_shift

    seq = ("GGACT" + "TTTTT") * 6                  # RRACH every 10 nt
    track = np.zeros(len(seq))
    for start in range(0, len(seq) - 4, 10):
        track[start + 2 - 1] = 0.9                 # A at start+2, shift +1
    shift = calibrate_shift([track], [seq], min_calls=4)
    assert shift == 1


def test_calibrate_shift_refuses_when_uninformative():
    from txmod.m6a import calibrate_shift

    with pytest.raises(RuntimeError):
        calibrate_shift([np.zeros(50)], ["T" * 50])



def test_isoform_divergence():
    assert isoform_divergence([0.5, -0.2, 0.1]) == pytest.approx(0.7)
    assert isoform_divergence([0.4]) is None      # needs >= 2 isoforms
