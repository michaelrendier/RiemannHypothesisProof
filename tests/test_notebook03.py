# SPDX-License-Identifier: GPL-3.0-only
"""Tests for the θ-clock, Gram crossings, echo spectrum, harness, margins, balance and PRIMEGAME (notebook 03)."""
import numpy as np
import pytest
from rhp import bigpath, clock, gram, echo, harness, lehmer, balance, primegame

ZEROS = np.loadtxt(__import__("pathlib").Path(__file__).resolve().parents[1] / "data" / "zeros_1_to_3600.csv")


def test_carrier_within_one_spacing():
    sh = clock.shift_in_spacings(ZEROS)
    assert np.abs(sh).max() < 1.0 and 0.1 < np.abs(sh).mean() < 0.4


def test_jitter_matches_shift():
    S = clock.jitter(ZEROS)
    assert abs(S.mean()) < 1e-2 and np.abs(S).max() < 1.0


def test_negative_crossings():
    r = gram.negative_crossings(3600.0)
    assert len(r["g"]) == 3066 and len(r["neg"]) == 197
    assert [int(r["n"][i]) for i in r["neg"][:3]] == [126, 134, 195]


def test_gram_intervals_adjacent_to_failures():
    r = gram.negative_crossings(3600.0)
    cnt = gram.gram_interval_counts(r["n"], r["g"], ZEROS)
    bad = set(np.nonzero(cnt != 1)[0].tolist())
    adj = set()
    for i in r["neg"]: adj |= {int(i) - 1, int(i)}
    assert bad <= adj and len(bad) == 393


def test_echo_is_the_explicit_formula():
    t = np.linspace(0.3, 4.1, 600)
    F = echo.zero_sum(t, ZEROS, 300.0)
    P, _, _ = echo.explicit_prediction(t, 300.0)
    assert np.abs(F - P).max() < 1e-8


def test_harness_controls():
    c = harness.candidates(ZEROS)
    s = {k: harness.score(v, ZEROS) for k, v in c.items()}
    assert s["carrier  θ(E)=(k−3/2)π"]["weyl"] < 1e-6 and s["carrier  θ(E)=(k−3/2)π"]["spacing_var"] < 1e-6
    lad = s["ladder (n+½)ħ  [H_NN, xp in a periodic u-box]"]
    assert lad["weyl"] > 50 and lad["locate"] > 20          # H_NN as coded fails
    assert 0.1 < s["GUE (unfolded)"]["spacing_var"] < 0.25 and s["Poisson"]["spacing_var"] > 0.8
    assert s["zeros (reference)"]["locate"] == 0.0


def test_lehmer_pair_near_1977():
    p = bigpath.pathway_big(2100.0, 0.02)
    zz = ZEROS[(ZEROS > 20) & (ZEROS < 2100)]
    m = lehmer.gap_margins(p, zz)
    j = int(np.argmin(m["rel"]))
    assert abs(m["gamma"][j] - 1977.1739) < 1e-2 and m["rel"][j] < 0.02
    t, z = lehmer.refine_peak(1977.173944, 1977.271446)
    assert abs(z - 0.0221) < 1e-3


def test_balance_helpers():
    from mpmath import mp, mpc, gamma, sin, pi as mpi, log
    s = mpc(0.3, 100)
    ref = float(log(abs(2**s * mpi ** (s - 1) * sin(mpi * s / 2) * gamma(1 - s))))
    assert abs(balance.ln_abs_chi(0.3, 100.0) - ref) < 1e-9
    D, A, B = balance.D_field(0.5, np.array([50.0, 120.0]))
    assert np.abs(D).max() < 1e-9


def test_primegame():
    g = primegame.primegame()
    assert [next(g) for _ in range(12)] == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    ev = primegame.events(7)
    assert [(e[0], e[1]) for e in ev][:5] == [(2, 19), (3, 69), (5, 280), (7, 707), (11, 2363)]
    # the video's counts = steps + firings of 15/14
    assert [ev[i][1] + ev[i][2][11] for i in (2, 3, 4)] == [281, 710, 2375]
    tr = primegame.trace(1500)
    assert set(tr[:, 4:].ravel().tolist()) == {0, 1}
    assert primegame.closing_table(13) == {d: 13 for d in range(1, 13)} and primegame.orbit_period(12, 8) == 3
