# SPDX-License-Identifier: GPL-3.0-only
"""Does ζ contain both Trees?  The Dirichlet-series bookkeeping, and where each Tree shows in the line spectrum.

Weights on n ≥ 2 (Λ = von Mangoldt):   Telperion  Λ(n)      →  T(s) = Σ Λ(n) n⁻ˢ = −ζ′/ζ(s)
                                        Laurelin   1 − Λ(n)  →  L(s) = Σ_{n≥2} (1 − Λ(n)) n⁻ˢ = ζ(s) − 1 + ζ′/ζ(s)
so  T + L = ζ − 1  exactly (the mingled whole is ζ; the split needs ζ′/ζ).      ESTABLISHED (Dirichlet algebra).
The interaction of the two is the convolution  1 ∗ Λ = ln :   −ζ′(s) = ζ(s)·T(s)  (Σ ln n · n⁻ˢ).             ESTABLISHED.
Everything below measures these on the line σ = σ₀ > 1 where the series converge and the t-Fourier transform has lines at u = ln n.
Labels: identity checks are :CALCULATED (floating point); line heights are :CALCULATED.
"""
from __future__ import annotations
import numpy as np
from mpmath import mp, mpc, zeta, log as mlog

mp.dps = 15


def vonmangoldt_table(N):
    lam = np.zeros(N + 1)
    sieve = np.ones(N + 1, bool); sieve[:2] = False
    for p in range(2, N + 1):
        if sieve[p]:
            sieve[p * p::p] = False
            q = p
            while q <= N:
                lam[q] = np.log(p); q *= p
    return lam


def parts(s):
    """(ζ, T, L) at s from the closed forms: ζ, −ζ′/ζ, ζ − 1 + ζ′/ζ."""
    s = mpc(s)
    z = zeta(s); dz = zeta(s, derivative=1)
    T = -dz / z
    return complex(z), complex(T), complex(z - 1 - T)


def parts_series(s, N=200000):
    """The same three from the weighted sums over n ≤ N (σ must be well above 1)."""
    lam = vonmangoldt_table(N); n = np.arange(2, N + 1, dtype=float)
    w = n ** (-s)
    return 1 + np.sum(w), np.sum(lam[2:] * w), np.sum((1 - lam[2:]) * w)


def identity_residuals(s=complex(2.0, 3.0), N=200000):
    """max of |T+L−(ζ−1)|, |−ζ′ − ζ·T|, and series-vs-closed-form."""
    z, T, L = parts(s)
    zs, Ts, Ls = parts_series(s, N)
    dz = complex(zeta(mpc(s), derivative=1))
    return dict(split=abs(T + L - (z - 1)), convolution=abs(-dz - z * T),
                series_zeta=abs(z - zs), series_T=abs(T - Ts), series_L=abs(L - Ls))


def line_heights(sigma=1.2, T=150.0, dt=0.05, ns=range(2, 31)):
    """t-Fourier lines (Hann window) of ζ, ζ′/ζ-part T, Laurelin L, ln ζ and ζ² at u = ln n.
    Returns {name: {n: height}} normalised so that a line of Dirichlet coefficient c_n has height c_n·n^{-σ} exactly (measured); the Dirichlet coefficient c_n n^{-σ}
    each is predicted to show:  ζ: 1,  T: Λ(n),  L: 1−Λ(n),  lnζ: Λ(n)/ln n,  ζ²: d(n)."""
    ts = np.arange(0, T, dt)
    zz = np.empty(len(ts), complex); TT = np.empty(len(ts), complex); lz = np.empty(len(ts), complex)
    for i, t in enumerate(ts):
        s = mpc(sigma, t)
        z = zeta(s); dz = zeta(s, derivative=1)
        zz[i] = complex(z); TT[i] = complex(-dz / z); lz[i] = complex(mlog(z))
    LL = zz - 1 - TT
    win = np.hanning(len(ts))
    out = {}
    for name, f in (("zeta", zz), ("T", TT), ("L", LL), ("lnzeta", lz), ("zeta2", zz * zz)):
        h = {}
        for n in ns:
            ph = np.exp(1j * ts * np.log(n))          # f contains n^{-σ} e^{-it ln n}; multiply back
            h[n] = abs(np.sum(f * ph * win)) / np.sum(win)
        out[name] = h
    return out


def divisor_count(n):
    return sum(1 for d in range(1, n + 1) if n % d == 0)


def expected_coeff(name, n, lam):
    return dict(zeta=1.0, T=lam[n], L=1 - lam[n], lnzeta=lam[n] / np.log(n),
                zeta2=float(divisor_count(n)))[name]
