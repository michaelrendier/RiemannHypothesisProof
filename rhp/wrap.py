# SPDX-License-Identifier: GPL-3.0-only
"""Primes on the spring: two different springs, kept apart.

SPRING 1 — the sieve spring (integers).  Wrap the integer line around a circle of circumference m: n sits at angle 2πn/m.
  n is a multiple of m exactly when it lands on the circle's 0 spot.  Every prime p is the 0 spot of its own circle, and n is
  prime iff it lands on the 0 spot of no circle p ≤ √n.  ESTABLISHED (this is the sieve of Eratosthenes / wheel factorisation).
  OURS-CALCULATED here: every twin pair (p, p+2), p > 3, straddles the 0 spot of the circumference-6 circle, and the allowed
  centres on the primorial circle P#  are  Π_{5≤p≤P}(p − 2)  of  P#  — the 0 spot gets thinner at every level.

SPRING 2 — the wheel spring (Euler product).  Prime p is the rotation e^{−it ln p}, radius p^{−σ}, height t.
  Two wheels p < q turn at rates ln p, ln q; they drift apart in phase by t·ln(q/p).  For a twin pair that is t·ln(1 + 2/p) ≈ 2t/p:
  twins ride the same circle and part slowly; the GAP is carried by the beat, not by the values.
  The Hardy spring (the zeros) turns at θ′(t) = ½ ln(t/2π) rad per unit height and tightens as t grows; the prime spikes of the
  explicit formula sit at t = ln n and crowd like eᵗ/t.  Zeros tighten logarithmically; primes crowd exponentially: that is the tension.
"""
from __future__ import annotations

import math
import numpy as np

from . import uns


def sieve(n):
    """Boolean primality array for 0..n (numpy sieve)."""
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return s


def primes_to(n):
    return np.flatnonzero(sieve(n))


def twin_pairs(n):
    """Lesser members p of twin pairs (p, p+2) with p + 2 ≤ n."""
    s = sieve(n)
    return np.flatnonzero(s[:-2] & s[2:])


def straddle_check(n):
    """Every twin pair with p > 3 has centre p + 1 ≡ 0 (mod 6).  Returns (pairs tested, pairs that straddle the 0 spot)."""
    tw = twin_pairs(n)
    tw = tw[tw > 3]
    return int(len(tw)), int(np.sum((tw + 1) % 6 == 0))


def primorial(p):
    out = 1
    for q in primes_to(p):
        out *= int(q)
    return out


def allowed_centres(P):
    """Residues c mod P# for which c−1 and c+1 are both coprime to P#: where a twin pair can be centred on the P#-circle."""
    m = primorial(P)
    ok = [c for c in range(m) if math.gcd(c - 1, m) == 1 and math.gcd(c + 1, m) == 1]
    return m, ok


def allowed_centre_count(P):
    """Π_{5≤p≤P}(p − 2) — the count of allowed centres on the P# circle (p = 2 and p = 3 fix centres ≡ 0 mod 6)."""
    out = 1
    for q in primes_to(P):
        if q >= 5:
            out *= int(q) - 2
    return out


def wheel_phase_lag(p, q, t):
    """Phase drift between the wheels of the primes p < q at height t:  t·ln(q/p)  (radians)."""
    return t * math.log(q / p)


def antiphase_height(p, q):
    """The first height at which the two wheels are in antiphase:  π / ln(q/p)."""
    return math.pi / math.log(q / p)


def hardy_turn_rate(t):
    """Radians per unit height of the Hardy spring:  θ′(t) ≈ ½ ln(t/2π)."""
    return np.array([uns.theta_prime(x) for x in np.atleast_1d(t)])


def hardy_pitch(t):
    """Height of one full 2π turn of the Hardy spring:  2π/θ′(t) — shrinks as t grows (the spring tightens)."""
    return 2 * np.pi / hardy_turn_rate(t)


def prime_density_in_t(t):
    """Prime spikes of the explicit formula sit at t = ln n; their density in t is eᵗ/t (π′(x)·x at x = eᵗ)."""
    t = np.asarray(t, dtype=float)
    return np.exp(t) / t


def residue_circle(P):
    """The primes' addresses on the P# circle: for p ≤ P that circle's 0 spot; for p > P the coprime residues."""
    m = primorial(P)
    res = sorted(c for c in range(m) if math.gcd(c, m) == 1)
    return m, res


# ── how far the wheels' grip reaches: the truncated Euler product against latitude ─────────────
def log_abs_euler(sigma, t, Ks, pmax=300_000):
    """ln|P_K(σ + it)| for P_K = Π_{p ≤ p_K}(1 − p^{−s})⁻¹, for each K in Ks (the K-th prime wheels), as an array
    (len(Ks), len(t)).  Branch-free (real part only).  Harmonics k ≥ 2 included up to k = 40 for p < 50, k ≤ 4 beyond."""
    t = np.asarray(t, float)
    P = primes_to(pmax)
    lp = np.log(P)
    Kmax = max(Ks)
    out = np.zeros((len(Ks), len(t)))
    acc = np.zeros(len(t))
    ks = sorted(Ks); j = 0
    for i in range(Kmax):
        p, l = P[i], lp[i]
        kmax = 40 if p < 50 else 4
        for k in range(1, kmax + 1):
            acc += np.exp(-k * sigma * l) * np.cos(k * t * l) / k
        if i + 1 == ks[j]:
            out[Ks.index(ks[j])] = acc
            j += 1
    return out


def euler_grip(sigma, t, Ks):
    """Median |ln|ζ| − ln|P_K|| over the t grid for each K: how well K wheels reproduce |ζ| at latitude σ."""
    from mpmath import zeta, mpc
    true = np.array([float(np.log(abs(complex(zeta(mpc(sigma, x)))))) for x in t])
    L = log_abs_euler(sigma, t, list(Ks))
    return np.median(np.abs(true[None, :] - L), axis=1)
