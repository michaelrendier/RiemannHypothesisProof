# SPDX-License-Identifier: GPL-3.0-only
"""Flow tests: the Cartesian Laplacian of log|ζ|, the islands (primes) against the rigid points (zeros), and the Fermat finite-tone wave.   :CALCULATED unless marked.

ESTABLISHED: on the (σ_ζ, t) plane Δ log|ζ(σ + it)| = 2π Σ_ρ δ(s − ρ) − 2π δ(s − 1): the zeros are point sources of the potential log|ζ|, the pole at s = 1 a sink; arg ζ is the stream function (the zeros are also vortices).
Islands: the depth field h(n) = d(n) − 2 ≥ 0 (d = number of divisors) is zero exactly at the primes: the points no factor stream-line passes through.  Spacing statistics separate two kinds of point set:
Poisson (independent, variance of the normalised spacing = 1) and GUE-like (level repulsion, variance ≈ 0.18).
Fermat side: for p ≡ 1 (mod n) the 2g = (n−1)(n−2) Frobenius eigenvalues are ω_j = −J_j = √p e^{iθ_j}; F(s) = p^{gs} ∏(1 − ω_j p^{−s}) is real on σ_ζ = ½ (up to a constant), with 2g zeros per period 2π/ln p.
"""
from __future__ import annotations
import math
import numpy as np


def flux_around(center, r, n=200, dps=20):
    """(1/2π) ∮ ∂_r log|ζ| r dθ around `center` (a complex number): +1 per enclosed zero, −1 for the pole at s = 1."""
    import mpmath as mp
    mp.mp.dps = dps; c = mp.mpmathify(center); tot = 0.0
    for k in range(n):
        th = 2 * math.pi * k / n; e = mp.e ** (1j * th); s = c + r * e
        tot += float(mp.re(mp.zeta(s, derivative=1) / mp.zeta(s) * e)) * r * (2 * math.pi / n)
    return tot / (2 * math.pi)


def unfold_zeros(g):
    """Smooth counting function N̄(t) = θ(t)/π + 1 (Riemann–Siegel θ, asymptotic form), so that the mean spacing of N̄(γ_k) is 1."""
    g = np.asarray(g, float); th = g / 2 * np.log(g / (2 * math.pi)) - g / 2 - math.pi / 8 + 1 / (48 * g); return th / math.pi + 1


def spacing_stats(x):
    """normalised nearest-neighbour spacings s = Δx/⟨Δx⟩: (variance, fraction below 0.1, mean)."""
    s = np.diff(np.sort(x)); s = s / s.mean(); return float(s.var()), float(np.mean(s < 0.1)), float(s.mean())


def number_variance(x, L, nwin=4000, rng=None):
    """variance of the number of points in windows of length L (unit mean density): Poisson → L; GUE → (ln(2πL) + γ_E + 1)/π² for L ≫ 1."""
    rng = rng or np.random.default_rng(0); x = np.sort(x); lo, hi = x[0], x[-1] - L; st = rng.uniform(lo, hi, nwin)
    c = np.searchsorted(x, st + L) - np.searchsorted(x, st); return float(c.var())


def frobenius_angles(p, n):
    """θ_j = arg(−J_j/√p) for the (n−1)(n−2) generic Jacobi sums of x^n + y^n = 1 over F_p (p ≡ 1 mod n)."""
    from .fermat_curve import jacobi_sums
    J = jacobi_sums(int(p), int(n)); return np.array([np.angle(-v / math.sqrt(p)) for (a, b), v in J.items() if (a + b) % n != 0])


def fermat_F(s, p, th):
    """F(s) = ∏_j (p^{s/2} − √p e^{iθ_j} p^{−s/2}) (the real-on-the-line normalisation of the L-polynomial), s a complex array."""
    L = math.log(p); out = np.ones_like(np.asarray(s, complex))
    for t in th: out = out * (np.exp(np.asarray(s, complex) * L / 2) - math.sqrt(p) * np.exp(1j * t) * np.exp(-np.asarray(s, complex) * L / 2))
    return out


def fermat_wave(u, p, th):
    """W_F(u) = −Σ_j cos(θ_j u/ln p): the real standing wave of the finite-tone prime wave at u = k ln p (equal amplitudes, 2g tones, RMS √g)."""
    L = math.log(p); u = np.asarray(u, float); return -np.cos(np.outer(u, th) / L).sum(axis=1)
