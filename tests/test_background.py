"""Tests for matched-background comparisons.

These target the failure modes the module exists to prevent: a composition
artefact read as selection, an opportunity artefact read as biology, and a
group-identity artefact read as discrimination. Each test constructs data where
the correct answer is known by construction.
"""
import numpy as np
import pytest

from txmod.background import (
    canonical_context,
    trinucleotide_rate_model,
    depletion_observed_expected,
    size_matched_enrichment,
)


# --- context folding ------------------------------------------------------- #

def test_canonical_context_pyrimidine_reference():
    assert canonical_context("ACA") == "ACA"
    # the same physical context read from the minus strand
    assert canonical_context("TGT") == "ACA"


def test_canonical_context_rejects_noncanonical():
    assert canonical_context("ANA") is None
    assert canonical_context("AC") is None


# --- rate model ------------------------------------------------------------ #

def test_rate_model_recovers_known_rates():
    # one sequence of repeated ACA, every third C mutated to T
    seq = "ACA" * 100
    seqs = {"t1": seq}
    muts = [("t1", i) for i in range(1, len(seq) - 1) if seq[i] == "C"][:10]
    m = trinucleotide_rate_model(seqs, muts)
    assert m.mutation_counts["ACA"] == 10
    assert m.context_counts["ACA"] > 0
    assert m.rates["ACA"] == pytest.approx(10 / m.context_counts["ACA"])


def test_rate_model_counts_a_position_once():
    """A position mutated in several samples is one mutated position."""
    seq = "ACA" * 100
    seqs = {"t1": seq}
    once = [("t1", i) for i in range(1, len(seq) - 1) if seq[i] == "C"][:10]
    m_once = trinucleotide_rate_model(seqs, once)
    m_dup = trinucleotide_rate_model(seqs, once + once + once)
    assert m_dup.mutation_counts["ACA"] == m_once.mutation_counts["ACA"] == 10


def test_rate_model_reproduces_low_rac_high_cpg_ordering():
    """The study's mechanism: RAC contexts mutate far less than CpG contexts."""
    gac = "GAC" * 200      # an m6A-type context
    acg = "ACG" * 200      # a CpG-containing context
    seqs = {"rac": gac, "cpg": acg}
    muts = ([("rac", i) for i in range(1, 599) if gac[i] == "A"][:2] +
            [("cpg", i) for i in range(1, 599) if acg[i] == "C"][:40])
    m = trinucleotide_rate_model(seqs, muts)
    rac_key = "GAC" if "GAC"[1] in "CT" else "GTC"
    assert m.total_rate("ACG") > m.total_rate(rac_key)


# --- depletion ------------------------------------------------------------- #

def _uniform_universe(n_seq=6, length=400, seed=0):
    rng = np.random.default_rng(seed)
    return {f"t{i}": "".join(rng.choice(list("ACGT"), size=length)) for i in range(n_seq)}


def test_depletion_no_signal_when_mutations_are_uniform():
    """Uniformly scattered mutations must give a corrected ratio spanning 1.0."""
    seqs = _uniform_universe()
    rng = np.random.default_rng(1)
    mp = {k: sorted(rng.choice(np.arange(1, len(v) - 1), size=120, replace=False))
          for k, v in seqs.items()}
    muts = [(k, p) for k, ps in mp.items() for p in ps]
    model = trinucleotide_rate_model(seqs, muts)
    sites = [(k, 200) for k in seqs]
    r = depletion_observed_expected(sites, seqs, mp, model)
    assert r["n_sites"] == len(seqs)
    lo, hi = r["ci95"]
    assert lo <= 1.0 <= hi, f"expected no depletion, got {r['corrected_ratio_of_ratios']:.3f}"


def test_depletion_detects_a_real_hole():
    """When the centre is genuinely emptied of mutations, the test must fire."""
    seqs = _uniform_universe()
    rng = np.random.default_rng(2)
    mp = {}
    for k, v in seqs.items():
        cand = [p for p in range(1, len(v) - 1) if abs(p - 200) > 25]
        mp[k] = sorted(rng.choice(cand, size=120, replace=False))
    muts = [(k, p) for k, ps in mp.items() for p in ps]
    model = trinucleotide_rate_model(seqs, muts)
    r = depletion_observed_expected([(k, 200) for k in seqs], seqs, mp, model)
    assert r["observed_centre"] == 0
    assert r["corrected_ratio_of_ratios"] < 0.5


def test_depletion_returns_raw_and_corrected_together():
    seqs = _uniform_universe(n_seq=3)
    mp = {k: list(range(10, 390, 7)) for k in seqs}
    muts = [(k, p) for k, ps in mp.items() for p in ps]
    model = trinucleotide_rate_model(seqs, muts)
    r = depletion_observed_expected([(k, 200) for k in seqs], seqs, mp, model)
    for key in ("raw_centre_flank_ratio", "corrected_ratio_of_ratios",
                "centre_oe", "flank_oe", "ci95", "p_centre_depletion"):
        assert key in r, f"{key} must always be reported"


# --- size-matched enrichment ---------------------------------------------- #

def test_size_matched_enrichment_removes_pure_opportunity_artefact():
    """Membership and event both driven by size, with no direct association."""
    rng = np.random.default_rng(3)
    n = 4000
    size = rng.lognormal(mean=1.0, sigma=1.0, size=n)
    p_member = 1 / (1 + np.exp(-(np.log(size) - 1.0)))
    p_event = 1 / (1 + np.exp(-(np.log(size) - 1.5)))
    member = rng.random(n) < p_member
    event = rng.random(n) < p_event
    res = size_matched_enrichment(member, event, {"log_size": np.log(size)})
    assert res["uncorrected"]["odds_ratio"] > 1.2, "confounded design should look enriched"
    log = res["corrected"]["logistic"]
    if "odds_ratio" in log:
        assert log["p"] > 0.01 or abs(log["odds_ratio"] - 1) < 0.25
    assert res["verdict"] in {"disappears", "attenuated", "underpowered"}


def test_size_matched_enrichment_keeps_a_real_association():
    rng = np.random.default_rng(4)
    n = 4000
    size = rng.lognormal(mean=1.0, sigma=0.5, size=n)
    event = rng.random(n) < 0.3
    # membership depends on the event itself, not on size
    p_member = np.where(event, 0.5, 0.1)
    member = rng.random(n) < p_member
    res = size_matched_enrichment(member, event, {"log_size": np.log(size)})
    log = res["corrected"]["logistic"]
    if "odds_ratio" in log:
        assert log["odds_ratio"] > 2.0
        assert log["p"] < 1e-6
    assert res["verdict"] == "survives"


def test_size_matched_enrichment_reports_both_blocks():
    rng = np.random.default_rng(5)
    n = 500
    res = size_matched_enrichment(rng.random(n) < 0.2, rng.random(n) < 0.3,
                                  {"x": rng.normal(size=n)})
    assert "uncorrected" in res and "corrected" in res
    assert "odds_ratio" in res["uncorrected"]
    assert "logistic" in res["corrected"], "the corrected estimate is never omitted"
