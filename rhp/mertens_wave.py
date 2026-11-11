# SPDX-License-Identifier: GPL-3.0-only
"""The Möbius / Mertens wave: the zeros seen from the 1/ζ side.   :CALCULATED (floating point) except where marked ESTABLISHED.

ψ-wave (the ζ'/ζ side, 'T'):  ψ(x) = x − Σ_ρ x^ρ/ρ − …              residue of −ζ'/ζ at ρ is 1
M-wave (the 1/ζ side):        M(x) = Σ_{n≤x} μ(n) = Σ_ρ x^ρ/(ρ ζ'(ρ)) − 2 + Σ_{n≥1} (−1)^{n−1}(2π/x)^{2n}/(n (2n)! ζ(2n+1))     (simple zeros assumed; ESTABLISHED, Titchmarsh §14.27)
Both waves run on the SAME carriers x^{½+iγ}; the M tone is the ψ tone rescaled by −1/ζ'(ρ), i.e. by the slope of ζ at the zero (amplitude 1/|ζ'|, phase −arg ζ').
ζ'(ρ) = −i e^{−iθ(γ)} Z'(γ), from Hardy's Z (data/zetaprime_25000.tsv, experiments/zetaprime_at_zeros.py).
Dirichlet-convolution identity (ESTABLISHED): Σ_{n≤x} M(x/n) = 1 for x ≥ 1, the counting function ⌊x⌋ and M being each other's inverse; each zero-tone is annihilated by convolution with 1
because Σ_n (x/n)^ρ = x^ρ ζ(ρ) = 0 — the tone is in the null space of the counting side.
"""
from __future__ import annotations
import math
import numpy as np
from .sedenion_slots import primes_upto


def mobius_table(X):
    """μ(n) for n = 0 … X (int8; μ(0) = 0): sign flip and rem division by every prime ≤ √X, squares zero, one leftover prime flips."""
    mu = np.ones(X + 1, np.int8); mu[0] = 0; rem = np.arange(X + 1, dtype=np.int32)
    for p in primes_upto(int(X ** .5)):
        p = int(p); mu[p::p] *= -1; rem[p::p] //= p
        mu[p * p::p * p] = 0
    mu[rem > 1] *= -1
    return mu


def load_zp(path="data/zetaprime_25000.tsv", K=None):
    """(γ, ζ'(ρ)) arrays."""
    g, z = [], []
    for line in open(path):
        t = line.split()
        if t and t[0].isdigit():
            g.append(float(t[1])); z.append(complex(float(t[2]), float(t[3])))
            if K and len(g) >= K: break
    return np.array(g), np.array(z)


def residues(g, zp):
    """c_k = 1/(ρ_k ζ'(ρ_k)), the M-tone coefficients (a conjugate pair gives 2 Re[c x^ρ])."""
    return 1.0 / ((0.5 + 1j * g) * zp)


def trivial_terms(x, nmax=12):
    """−2 + Σ_{n=1}^{nmax} (−1)^{n−1} (2π/x)^{2n}/(n (2n)! ζ(2n+1))  (the s = 0 and trivial-zero residues)."""
    from math import factorial
    zeta_odd = lambda m: sum(k ** -m for k in range(1, 200)) + 200 ** (1 - m) / (m - 1)
    x = np.asarray(x, float); out = -2.0 * np.ones_like(x)
    for n in range(1, nmax + 1):
        out += (-1) ** (n - 1) * (2 * math.pi / x) ** (2 * n) / (n * factorial(2 * n) * zeta_odd(2 * n + 1))
    return out


def m_tones(x, g, c, K=None, chunk=500):
    """Σ_{k≤K} 2 Re[c_k x^{ρ_k}]  for each x (x ≥ 1)."""
    K = K or len(g); g = g[:K]; c = c[:K]; x = np.asarray(x, float); lx = np.log(x); out = np.empty(len(x))
    for a in range(0, len(x), chunk):
        sub = lx[a:a + chunk]; out[a:a + chunk] = 2 * np.real((np.exp(np.outer(sub, 0.5 + 1j * g[:K])) * c[:K]).sum(axis=1))
    return out


def m_reconstruct(x, g, c, K=None):
    return m_tones(x, g, c, K) + trivial_terms(x)


def mu_ping(n, g, c, K):
    """Unit-window jump of the full explicit formula (zero-tones plus the s = 0 and trivial-zero terms) across [n−½, n+½]: ≈ μ(n) ∈ {−1, 0, +1}."""
    n = np.asarray(n, float); return m_reconstruct(n + 0.5, g, c, K) - m_reconstruct(np.maximum(n - 0.5, 1.0), g, c, K)


def convolution_residual(x, g, c, K):
    """Σ_{n ≤ x} M_K(x/n) − 1 for an integer x (should be 0)."""
    n = np.arange(1, int(x) + 1); y = x / n
    return float(m_reconstruct(y, g, c, K).sum() - 1.0)


def classify(jump):
    """nearest of {−1, 0, +1}."""
    return np.clip(np.rint(jump), -1, 1).astype(int)
