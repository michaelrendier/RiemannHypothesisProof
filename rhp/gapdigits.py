# SPDX-License-Identifier: GPL-3.0-only
"""The Fermat domain as the spacing of the primes: the composites between consecutive primes and their digit weight.   :CALCULATED (exact integer counts) unless marked.

Position x = 10^D (D = the digit coordinate).  Between consecutive primes p < q lie g − 1 = q − p − 1 composites (the Fermat domain: numbers with factors).  Their 'weight' in the sense of the
2026-04 conversation is the number of digits: each n has d(n) = ⌊log₁₀ n⌋ + 1 digits; the gap's digit weight is W = Σ_{p<n<q} d(n).  Predictions from the prime number theorem (smooth, ESTABLISHED):
  mean gap ≈ ln x,  composites per gap ≈ ln x − 1,  composite density 1 − 1/ln x,  digit weight per gap ≈ (ln x − 1)(D + 1).
A gap that contains a power of ten has composites of two different digit counts: a digit transition (Gemini's 'digit transition lag').
"""
from __future__ import annotations
import math
import numpy as np
from .sedenion_slots import primes_upto


def decade_stats(P, m):
    """For primes in [10^m, 10^(m+1)): count, mean gap, mean composites per gap, mean digit weight per gap, empirical gap entropy (nats), digit-transition gaps."""
    lo, hi = 10 ** m, 10 ** (m + 1); i0, i1 = np.searchsorted(P, lo), np.searchsorted(P, hi)
    p = P[i0:i1 + 1].astype(np.int64); g = np.diff(p); w = np.zeros(len(g))
    # digit weight of the composites in each gap: all n in (p, q) have d(n) = m+1 digits except in the gap that straddles 10^(m+1) (excluded by the window)
    w = (g - 1) * (m + 1)
    cnt = np.bincount(g); pr = cnt[cnt > 0] / cnt.sum(); H = float(-(pr * np.log(pr)).sum())
    x = math.sqrt(lo * hi)
    return dict(m=m, primes=int(len(g)), mean_gap=float(g.mean()), mean_composites=float((g - 1).mean()), mean_digit_weight=float(w.mean()), gap_entropy_nats=H,
                pnt_mean_gap=math.log(x), pnt_mean_composites=math.log(x) - 1, composite_density=1 - len(g) / (hi - lo))


def digit_transition_gaps(P, mmax):
    """For each power of ten 10^m (m = 1 … mmax): the gap (p, q) that contains it, its length and its composites split by digit count (m digits below 10^m, m+1 at and above)."""
    out = []
    for m in range(1, mmax + 1):
        t = 10 ** m; i = np.searchsorted(P, t, side="right"); p, q = int(P[i - 1]), int(P[i])
        below = t - p - 1 if p < t else 0; above = q - t if q > t else 0                # composites n with p < n < 10^m  and  10^m ≤ n < q
        out.append(dict(m=m, p=p, q=q, gap=q - p, composites_below=below, composites_at_or_above=above - (1 if t == q else 0)))
    return out
