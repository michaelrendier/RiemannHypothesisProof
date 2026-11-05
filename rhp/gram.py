# SPDX-License-Identifier: GPL-3.0-only
"""Gram points, the real-axis crossings of the pathway, and Gram's law.

The flat pathway meets the real axis only where Im ζ = 0, i.e. at Gram points θ(g_n) = nπ (or at a zero).
There ζ(½ + i g_n) = (−1)^n Z(g_n): POSITIVE if Gram's law holds, NEGATIVE if it fails.  (ESTABLISHED; :CALCULATED below)
"""
from __future__ import annotations
import numpy as np
from mpmath import mpf, siegelz
from . import bigpath, clock


def gram_points(T=3600.0, nmax=3200):
    n = np.arange(0, nmax)
    t = np.maximum(17.8456 + n * 2 * np.pi / np.log(np.maximum(n, 1) + 8), 17.9)
    for _ in range(100):
        t = t - (bigpath.theta_fast(t) - n * np.pi) / clock.theta_prime_num(t)
    ok = (np.abs(bigpath.theta_fast(t) - n * np.pi) < 1e-9) & (t < T)
    return n[ok], t[ok]


def real_axis_values(n, g, refine_below=0.02):
    """ζ(½ + i g_n) = (−1)^n Z(g_n). Z from Riemann–Siegel; any |Z| < refine_below is re-evaluated with mpmath."""
    Z = bigpath.Z_rs(g)
    for i in np.nonzero(np.abs(Z) < refine_below)[0]:
        Z[i] = float(siegelz(mpf(g[i])))
    return ((-1.0) ** n) * Z


def negative_crossings(T=3600.0):
    n, g = gram_points(T)
    v = real_axis_values(n, g)
    neg = np.nonzero(v < 0)[0]
    return dict(n=n, g=g, value=v, neg=neg)


def gram_interval_counts(n, g, zeros):
    """Zeros in each Gram interval [g_n, g_{n+1})."""
    return np.array([int(np.sum((zeros >= g[i]) & (zeros < g[i + 1]))) for i in range(len(g) - 1)])
