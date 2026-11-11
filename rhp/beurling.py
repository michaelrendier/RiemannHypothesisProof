# SPDX-License-Identifier: GPL-3.0-only
"""Beurling generalized primes: what the prime number theorem needs, measured.   :CALCULATED (exact integer/float computation from this file).

A Beurling system: a multiset P of reals > 1 (the 'primes'); the generalized integers are all finite products; ζ_P(s) = ∏_{q∈P}(1 − q^{−s})^{−1} = Σ n^{−s}.
The multiplicative structure is that of ℕ; nothing else is imposed (no functional equation, no lattice).
Systems built here
  true      the primes (control)
  gated     each prime p taken m(p) = 1 + sgn cos(τ ln p) ∈ {0, 2} times.  Density 1 on average; the gate puts a pair of singularities on σ = 1 at s = 1 ± ikτ (k odd)
            ψ_P(x)/x = 1 + (4/π) Σ_{j≥0} (−1)^j/(2j+1) · Re[x^{i(2j+1)τ}/(1 + i(2j+1)τ)] + o(1)       (sgn cos expanded in harmonics, then PNT for the true primes)
  random    each integer n ≥ 2 is a prime with probability 1/ln n (Cramér); E ψ = x, Var ψ ≈ x ln x
  Z[√−5]    the non-UFD level: irreducible elements against prime ideals
"""
from __future__ import annotations
import numpy as np
from .sedenion_slots import primes_upto


# ---------------------------------------------------------------- the systems
def gated_multiplicity(primes, tau):
    return (1 + np.sign(np.cos(tau * np.log(primes)))).astype(np.int64)


def psi_beurling(primes, mult, xs):
    """ψ_P(x) = Σ_{k ≥ 1} Σ_{q ∈ P, q^k ≤ x} ln q (each copy counted).  primes sorted ascending."""
    xs = np.asarray(xs, float); w = np.cumsum(mult * np.log(primes)); out = np.zeros(len(xs))
    for k in range(1, int(np.log2(xs.max())) + 1):
        idx = np.searchsorted(primes, xs ** (1.0 / k), side="right")
        out += np.where(idx > 0, w[np.maximum(idx - 1, 0)], 0.0)
    return out


def gated_prediction(xs, tau, J=80):
    """ψ_P(x)/x predicted from the harmonics of the gate (no fitting)."""
    u = np.log(np.asarray(xs, float)); s = np.zeros(len(u))
    for j in range(J):
        k = 2 * j + 1; s += (4 / np.pi) * (-1) ** j / k * np.real(np.exp(1j * k * tau * u) / (1 + 1j * k * tau))
    return 1 + s


def generalized_integers(primes, mult, X):
    """a[n] = number of formal products equal to n (n ≤ X) in the system; N(x) = cumsum a."""
    a = np.zeros(X + 1, np.int64); a[1] = 1
    for q, m in zip(primes[primes <= X].tolist(), mult[primes <= X].tolist()):
        for _ in range(m):
            new = a.copy(); qk = q
            while qk <= X: new[qk::qk] += a[1:X // qk + 1]; qk *= q
            a = new
    return a


def random_system(X, seed, chunk=10_000_000):
    """Cramér primes up to X: sorted array of the selected integers."""
    rng = np.random.default_rng(seed); out = []
    for s in range(2, X + 1, chunk):
        n = np.arange(s, min(s + chunk, X + 1)); out.append(n[rng.random(len(n)) < 1.0 / np.log(n)])
    return np.concatenate(out)


def rms_over(psi, xs, lo, hi):
    """RMS of (ψ(x) − x)/√x over the sample points lo ≤ x ≤ hi."""
    m = (xs >= lo) & (xs <= hi); return float(np.sqrt(np.mean(((psi[m] - xs[m]) / np.sqrt(xs[m])) ** 2)))


# ---------------------------------------------------------------- the non-UFD level Z[√−5]  (class number 2, discriminant −20)
def z5_irreducibles(X):
    """(prime ideals of norm ≤ X, principal prime ideals, irreducible elements of norm ≤ X up to units ±1).
    Prime ideals: p ≡ 1, 3, 7, 9 (mod 20) split (two of norm p), principal iff p ≡ 1, 9; p ≡ 11, 13, 17, 19 inert (one, norm p², principal); 2 ramified (norm 2, not principal); 5 ramified (norm 5, principal).
    Class group ℤ/2: two non-principal ideals multiply to a principal one, so the irreducibles are exactly the principal prime ideals and the products of two non-principal prime ideals."""
    p = primes_upto(X); r = p % 20
    split_pr = p[(r == 1) | (r == 9)]; split_np = p[(r == 3) | (r == 7)]; inert = p[(r >= 11) & (p * p <= X)] if False else p[np.isin(r, (11, 13, 17, 19)) & (p.astype(np.int64) ** 2 <= X)]
    principal = np.concatenate([np.repeat(split_pr, 2), inert.astype(np.int64) ** 2, [5] if X >= 5 else []]).astype(np.int64)
    nonprin = np.sort(np.concatenate([np.repeat(split_np, 2), [2] if X >= 2 else []])).astype(np.int64)
    n_prin = int(len(principal)); n_all = n_prin + len(nonprin)
    O = int(np.searchsorted(nonprin, X // nonprin, side="right").sum()); D = int((nonprin * nonprin <= X).sum())
    pairs = (O + D) // 2
    return n_all, n_prin, n_prin + pairs


def z5_irreducibles_brute(X):
    """Direct count of irreducible elements of ℤ[√−5] of norm ≤ X up to ±1, by testing every proper divisor (small X only)."""
    import math
    els = {}
    r = int(math.isqrt(X)) + 1
    for a in range(-r, r + 1):
        for b in range(-r, r + 1):
            n = a * a + 5 * b * b
            if 1 < n <= X: els.setdefault(n, []).append((a, b))
    def divides(be, al):
        c, d = be; a, b = al; N = c * c + 5 * d * d
        return (a * c + 5 * b * d) % N == 0 and (b * c - a * d) % N == 0
    cnt = 0
    for n, lst in els.items():
        for al in lst:
            if (al[0] < 0) or (al[0] == 0 and al[1] < 0): continue
            red = False
            for d in range(2, n):
                if n % d == 0 and d in els and any(divides(be, al) for be in els[d]): red = True; break
            cnt += not red
    return cnt
