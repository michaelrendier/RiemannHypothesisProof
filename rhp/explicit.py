# SPDX-License-Identifier: GPL-3.0-only
"""The explicit formulas as two-sided identities: tones on one side, primes (or point counts) on the other.   :CALCULATED.

Guinand–Weil (ESTABLISHED; for h even, analytic in |Im r| < ½ + δ, g(u) = (1/2π) ∫ h(r) e^{−iru} dr):
    Σ_ρ h(γ_ρ) = h(i/2) + h(−i/2) − g(0) ln π + (1/2π) ∫ h(r) Re ψ(¼ + ir/2) dr − 2 Σ_{n ≥ 1} Λ(n) n^{−½} g(ln n)
with ρ = ½ + iγ over all non-trivial zeros (γ and −γ), ψ = Γ'/Γ.  We use h(r) = e^{−(r−r₀)²/2σ²} + e^{−(r+r₀)²/2σ²}, for which g(u) = (σ/√2π) e^{−σ²u²/2}·2cos(r₀u).
Fermat side (ESTABLISHED, Weil; checked against exact counts in rhp.gf): #{xⁿ + yⁿ = 1 over F_{p^k}} = p^k − Σ_{a,b ≥ 1} (−J_p(a,b))^k for every k ≥ 1.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.special import digamma


def gaussian_pair(r0, sigma):
    h = lambda r: np.exp(-(np.asarray(r, float) - r0) ** 2 / (2 * sigma ** 2)) + np.exp(-(np.asarray(r, float) + r0) ** 2 / (2 * sigma ** 2))
    g = lambda u: sigma / math.sqrt(2 * math.pi) * np.exp(-sigma ** 2 * np.asarray(u, float) ** 2 / 2) * 2 * np.cos(r0 * np.asarray(u, float))
    return h, g


def guinand_weil(gamma, lam, r0, sigma, K=None, gridstep=None):
    """Both sides of Guinand–Weil for the Gaussian pair at (r0, σ).  gamma: positive zero ordinates; lam: Λ(n) for n = 0 … N (N must cover ln n ≲ 6/σ).
    Returns dict(zero_sum, gamma_term, pole_term, lnpi_term, prime_sum, rhs, residual)."""
    h, g = gaussian_pair(r0, sigma)
    gam = np.asarray(gamma if K is None else gamma[:K], float); zero_sum = float(2 * np.sum(h(gam)))
    pole = float(2 * np.real(np.exp(-(0.5j - r0) ** 2 / (2 * sigma ** 2)) + np.exp(-(0.5j + r0) ** 2 / (2 * sigma ** 2))))           # h(i/2) + h(−i/2)
    lnpi = float(g(0.0) * math.log(math.pi))
    step = gridstep or sigma / 40; r = np.arange(-(r0 + 10 * sigma), r0 + 10 * sigma, step)
    gterm = float(np.sum(h(r) * np.real(digamma(0.25 + 0.5j * r))) * step / (2 * math.pi))
    n = np.arange(1, len(lam)); w = lam[1:] / np.sqrt(n) * g(np.log(n)); prime = float(2 * np.sum(w))
    rhs = pole - lnpi + gterm - prime
    return dict(zero_sum=zero_sum, pole_term=pole, lnpi_term=lnpi, gamma_term=gterm, prime_sum=prime, rhs=rhs, residual=abs(zero_sum - rhs))
