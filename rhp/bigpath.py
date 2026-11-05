# SPDX-License-Identifier: GPL-3.0-only
"""Fast pathway ζ(½+it) to heights with thousands of zeros: Riemann–Siegel (ESTABLISHED), numpy.

Z(t) = 2 Σ_{n ≤ N} n^{-1/2} cos(θ(t) − t ln n) + R(t),  N = ⌊√(t/2π)⌋,
R(t) ≈ (−1)^{N−1} (t/2π)^{−1/4} Ψ(p),  p = √(t/2π) − N,  Ψ(p) = cos(2π(p² − p − 1/16)) / cos(2πp)
(first correction only; error O(t^{−5/4}), checked against mpmath in tests).
θ by Stirling for t ≥ 20, mpmath below.  ζ(½+it) = e^{−iθ} Z.
"""
from __future__ import annotations
import numpy as np
from mpmath import mpf, siegeltheta, zetazero

from . import uns


def theta_fast(t):
    t = np.asarray(t, float)
    out = np.empty_like(t)
    big = t >= 20.0
    tb = t[big]
    out[big] = (tb / 2) * np.log(tb / (2 * np.pi)) - tb / 2 - np.pi / 8 + 1 / (48 * tb) + 7 / (5760 * tb**3) + 31 / (80640 * tb**5)
    out[~big] = [float(siegeltheta(mpf(x))) for x in t[~big]]
    return out


def Z_rs(t):
    """Hardy Z by Riemann–Siegel with the first correction term. Vectorised in t."""
    t = np.asarray(t, float)
    th = theta_fast(t)
    a = np.sqrt(t / (2 * np.pi))
    N = np.floor(a).astype(int)
    Z = np.zeros_like(t)
    for n in range(1, int(N.max()) + 1):
        m = N >= n
        Z[m] += np.cos(th[m] - t[m] * np.log(n)) / np.sqrt(n)
    Z *= 2
    p = a - N
    psi = np.cos(2 * np.pi * (p * p - p - 1 / 16)) / np.cos(2 * np.pi * p)
    Z += ((-1.0) ** (N - 1)) * (t / (2 * np.pi)) ** -0.25 * psi
    return Z


def pathway_big(T=3600.0, dt=0.02, t0=20.0):
    """Sample t ∈ [t0, T] (t0 ≥ 20 so Stirling and RS apply). Returns dict like rhp.uns.pathway."""
    t = np.arange(t0, T + dt / 2, dt)
    th = theta_fast(t)
    Z = Z_rs(t)
    z = np.exp(-1j * th) * Z
    with np.errstate(divide="ignore"):
        u = np.log10(np.abs(z))
    return dict(t=t, theta=th, Z=Z, zeta=z, u=u, phi_signed=-th)


def zeros_from_Z(t, Z, iters=45):
    """Zeros by sign change of Z on the grid, refined by bisection on the RS function."""
    i = np.nonzero(np.sign(Z[:-1]) * np.sign(Z[1:]) < 0)[0]
    lo, hi = t[i].copy(), t[i + 1].copy()
    flo = Z[i].copy()
    for _ in range(iters):
        mid = (lo + hi) / 2
        fm = Z_rs(mid)
        left = np.sign(fm) == np.sign(flo)
        lo = np.where(left, mid, lo)
        flo = np.where(left, fm, flo)
        hi = np.where(left, hi, mid)
    return (lo + hi) / 2


def N_smooth(T):
    """θ(T)/π + 1 (the smooth part of N(T)); N(T) = N_smooth + S(T)."""
    return theta_fast(np.atleast_1d(T))[0] / np.pi + 1


def refine_zeros_mp(z0, dps=15):
    """Polish approximate zeros with mpmath's Hardy Z (secant). Slow (~0.05 s each) — cached to data/ by the notebook."""
    from mpmath import mp, findroot, siegelz
    mp.dps = dps
    out = []
    for x in z0:
        r = findroot(lambda s: siegelz(s), (mpf(x) - 0.01, mpf(x) + 0.01), solver="secant", tol=1e-12, maxsteps=30)
        out.append(float(r))
    return np.array(out)
