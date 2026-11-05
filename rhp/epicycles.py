# SPDX-License-Identifier: GPL-3.0-only
"""The machinery that draws the pathway: prime wheels (σ>1) and the explicit formula (primes ↔ zeros).

log ζ(s) = Σ_p Σ_k p^{−ks}/k   (ESTABLISHED, σ > 1): each (p, k) is a wheel of radius p^{−kσ}/k turning at k·ln p per unit t.
ψ(x) = x − Σ_ρ x^ρ/ρ − ln 2π − ½ ln(1 − x^{−2})   (ESTABLISHED; ρ = ½ + iγ used for the zeros computed here).
A Gaussian factor e^{−(γ w)²/2} on each zero gives a continuous (smoothed) sieve; w → 0 recovers the jumps.
"""
from __future__ import annotations
import numpy as np
from mpmath import mp, mpc, zeta, log as mlog

mp.dps = 15


def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]


def logzeta_true(sig, t):
    """Continuous log ζ(σ+it) for σ > 1 (branch fixed by unwrapping, log ζ(σ) real at t = 0)."""
    z = np.array([complex(zeta(mpc(sig, x))) for x in t])
    return np.log(np.abs(z)) + 1j * np.unwrap(np.angle(z))


def logzeta_wheels(sig, t, K, kmax=None):
    """Sum of the first K prime wheels (all harmonics up to kmax chosen so the tail is below 1e-12)."""
    P = primes_upto(max(100, int(K * 20)))[:K]
    t = np.asarray(t, float)
    out = np.zeros(len(t), complex)
    for p in P:
        km = kmax or max(1, int(np.ceil(12 * np.log(10) / (sig * np.log(p)))))
        for k in range(1, km + 1):
            out += np.exp(-k * (sig + 1j * t) * np.log(p)) / k
    return out


def wheel_table(sig, t, Ks):
    """Max |log ζ − Σ wheels| over t for each K."""
    true = logzeta_true(sig, t)
    return [(K, float(np.max(np.abs(true - logzeta_wheels(sig, t, K))))) for K in Ks]


def psi_exact(x):
    """Chebyshev ψ(x) = Σ_{n ≤ x} Λ(n), exact."""
    x = np.atleast_1d(x).astype(float)
    n = int(x.max())
    lam = np.zeros(n + 1)
    for p in primes_upto(n):
        q = int(p)
        while q <= n:
            lam[q] = np.log(p); q *= p
    cs = np.cumsum(lam)
    return cs[np.floor(x).astype(int)]


def psi_zeros(x, gam, w=0.0):
    """ψ(x) − (trivial terms) from the zeros γ with smoothing width w (in ln x)."""
    x = np.atleast_1d(x).astype(float)
    lx = np.log(x)
    g = np.asarray(gam, float)
    wt = np.exp(-(g * w) ** 2 / 2)
    rho = 0.5 + 1j * g
    out = np.empty_like(x)
    for i, (xi, l) in enumerate(zip(x, lx)):
        out[i] = xi - 2 * np.real(np.sum(wt * np.exp(rho * l) / rho)) - np.log(2 * np.pi) - 0.5 * np.log(1 - xi**-2)
    return out


def afe_parts(sig, t):
    """Approximate functional equation (Hardy–Littlewood, ESTABLISHED): ζ(s) ≈ A + χ(s)·B with
    A = Σ_{n≤N} n^{−s}, B = Σ_{n≤N} n^{s−1}, N = ⌊√(t/2π)⌋, χ(s) = 2^s π^{s−1} sin(πs/2) Γ(1−s).
    At σ = ½: |χ| = 1 and B = conj(A), so the two halves have EQUAL modulus: ζ = A + e^{−2iθ}·conj(A).
    Returns (A, chiB, chi, zeta_true) as complex arrays. The truncation error of A + χB is O(t^{−1/4}) (no correction terms)."""
    from mpmath import gamma, sin, pi as mpi, mpf
    t = np.asarray(t, float)
    A = np.zeros(len(t), complex); B = np.zeros(len(t), complex); chi = np.zeros(len(t), complex); zt = np.zeros(len(t), complex)
    for i, x in enumerate(t):
        s = mpc(sig, x)
        N = max(1, int(np.floor(np.sqrt(x / (2 * np.pi)))))
        n = np.arange(1, N + 1)
        A[i] = np.sum(n ** (-complex(s)))
        B[i] = np.sum(n ** (complex(s) - 1))
        chi[i] = complex(2 ** s * mpi ** (s - 1) * sin(mpi * s / 2) * gamma(1 - s))
        zt[i] = complex(zeta(s))
    return A, chi * B, chi, zt
