# SPDX-License-Identifier: GPL-3.0-only
"""The Fermat cubic x³ + y³ = z³ as an L-function: a_p from Fermat-side counts, its zeros from the functional equation, and the zero-phase scan that finds them from the a_p alone.   :CALCULATED.

E: x³ + y³ = z³ is the elliptic curve of conductor N = 27 (isogeny class 27a, CM by ℤ[ω]); L(E, s) = ∏_p (1 − a_p p^{−s} + p^{1−2s})^{−1} (a_3 = 0, bad prime), a_p = p + 1 − #E(F_p).
a_p = 0 for p ≡ 2 (mod 3); for p ≡ 1 (mod 3) write 4p = L² + 27M² with L ≡ 1 (mod 3): a_p = −L (Gauss; checked against brute-force counts).  Analytic centre s = 1.
Completed function Λ(s) = (√27/2π)^s Γ(s) L(E, s) = w Λ(2 − s), w = +1 (checked: Λ(1 + it) is real).  Approximate functional equation:
    Λ(s) = Σ_n a_n [ (√27/2πn)^s Γ(s, 2πn/√27) + w (√27/2πn)^{2−s} Γ(2 − s, 2πn/√27) ].
Zero-phase scan (the Landau–Gonek detector used for ζ in docs/zerophase), no pole so no main term to subtract:   A_E(t) = | Σ_{n ≤ x} Λ_E(n) n^{−1−it} | / ln x,   Λ_E(p^k) = (α^k + β^k) ln p.
"""
from __future__ import annotations
import math
import numpy as np
from .sedenion_slots import primes_upto

N_COND = 27


def ap_table(X):
    """array a[0..X]: a[p] for primes p ≤ X (0 elsewhere)."""
    a = np.zeros(X + 1, np.int64); isp = np.zeros(X + 1, bool); isp[primes_upto(X)] = True
    M = 1
    while 27 * M * M <= 4 * X:
        Lmax = int(math.isqrt(4 * X - 27 * M * M)); L = np.arange(-Lmax - 3, Lmax + 4); L = L[(L % 3) == 1]
        v = L * L + 27 * M * M; ok = (v % 4 == 0) & (v // 4 <= X); L, p = L[ok], v[ok] // 4
        good = isp[p]; a[p[good]] = -L[good]; M += 1
    return a


def an_table(N, a):
    """multiplicative coefficients a_n of L(E, s) for n ≤ N (uses a_p from `a`, which must reach N)."""
    spf = np.zeros(N + 1, np.int64)
    for p in range(2, N + 1):
        if spf[p] == 0: spf[p::p] = np.where(spf[p::p] == 0, p, spf[p::p])
    an = np.zeros(N + 1, float); an[1] = 1
    for n in range(2, N + 1):
        p = int(spf[n]); m = n; k = 0
        while m % p == 0: m //= p; k += 1
        ap = float(a[p]) if p != 3 else 0.0; c = [1.0, ap]
        for _ in range(2, k + 1): c.append(ap * c[-1] - p * c[-2] if p != 3 else 0.0)
        an[n] = c[k] * an[m]
    return an


def lam_von_mangoldt(P, a, x):
    """(n, Λ_E(n)) for all prime powers n = p^k ≤ x, p ≠ 3 (Λ_E(p^k) = (α^k + β^k) ln p; s_1 = a_p, s_2 = a_p² − 2p, s_k = a_p s_{k−1} − p s_{k−2})."""
    ns, ls = [], []
    for p in P:
        p = int(p)
        if p == 3: continue
        ap = float(a[p]); s0, s1 = 2.0, ap; q = p; lp = math.log(p)
        while q <= x:
            ns.append(q); ls.append(s1 * lp); q *= p; s0, s1 = s1, ap * s1 - p * s0
    return np.array(ns, float), np.array(ls, float)


def scan(ts, ns, ls, chunk=200):
    """A_E(t) = |Σ Λ_E(n) n^{−1−it}| / ln x on the grid ts."""
    lx = math.log(ns.max()); ln = np.log(ns); w = ls / ns; out = np.empty(len(ts))
    for i in range(0, len(ts), chunk): out[i:i + chunk] = np.abs((np.exp(-1j * np.outer(ts[i:i + chunk], ln)) * w).sum(axis=1))
    return out / lx


def completed_L(t, an, nterms=100, dps=50):
    """Λ(1 + it) by the approximate functional equation (mpmath), w = +1."""
    import mpmath as mp
    mp.mp.dps = dps; s = mp.mpc(1, t); c = mp.sqrt(N_COND) / (2 * mp.pi); tot = mp.mpc(0)
    for n in range(1, nterms + 1):
        if an[n] == 0: continue
        x = 2 * mp.pi * n / mp.sqrt(N_COND)
        tot += an[n] * ((c / n) ** s * mp.gammainc(s, x) + (c / n) ** (2 - s) * mp.gammainc(2 - s, x))
    return tot


def find_zeros(an, tmin, tmax, step, nterms=100, dps=50):
    """Zeros of Re Λ(1 + it) (Λ is real on the line): sign changes on a grid, refined by bisection."""
    import mpmath as mp
    ts = np.arange(tmin, tmax + 1e-12, step); vals = [float(mp.re(completed_L(t, an, nterms, dps)) * mp.e ** (mp.pi * t / 2)) for t in ts]
    out = []
    for i in range(len(ts) - 1):
        if vals[i] * vals[i + 1] < 0:
            lo, hi, flo = ts[i], ts[i + 1], vals[i]
            for _ in range(40):
                mid = 0.5 * (lo + hi); fm = float(mp.re(completed_L(mid, an, nterms, dps)) * mp.e ** (mp.pi * mid / 2))
                if flo * fm <= 0: hi = mid
                else: lo, flo = mid, fm
            out.append(0.5 * (lo + hi))
    return np.array(out), ts, np.array(vals)
