# SPDX-License-Identifier: GPL-3.0-only
"""The echo spectrum: F_Γ(t) = Σ_{γ>0} w(γ) cos(γ t),  w = e^{−(γ/Γ)²/2}  against the explicit formula (Weil).

Σ_γ h(γ) = h(i/2) + h(−i/2) − g(0) ln π + (1/2π) ∫ h(r) Re ψ(¼ + ir/2) dr − 2 Σ_n Λ(n) n^{−1/2} g(ln n),
g(u) = (1/2π) ∫ h(r) e^{−iru} dr.  (ESTABLISHED.)  With h_t(r) = cos(rt) e^{−(r/Γ)²/2}: g_t(u) = ½[G(u−t) + G(u+t)],
G(u) = Γ (2π)^{−1/2} e^{−Γ²u²/2}.  So F_Γ(t) = ½·(RHS).  The pole term cosh(t/2) is the 'path of least primes'.
The prime powers n = p^k appear as NEGATIVE Gaussian bumps at t = ln n of area ½ Λ(n)/√n.
"""
from __future__ import annotations
import numpy as np
from mpmath import mp, digamma, mpc

mp.dps = 15


def lam_table(nmax):
    lam = np.zeros(nmax + 1)
    s = np.ones(nmax + 1, bool); s[:2] = False
    for i in range(2, int(nmax**0.5) + 1):
        if s[i]: s[i * i::i] = False
    for p in np.nonzero(s)[0]:
        q = int(p)
        while q <= nmax:
            lam[q] = np.log(p); q *= int(p)
    return lam


def zero_sum(t, zeros, Gamma):
    t = np.asarray(t, float)
    w = np.exp(-(zeros / Gamma) ** 2 / 2)
    return np.array([np.sum(w * np.cos(zeros * x)) for x in t])


def _re_psi(r):
    """Re ψ(¼ + i r/2): mpmath for |r| ≤ 40, asymptotic series beyond."""
    r = np.asarray(r, float)
    out = np.empty_like(r)
    small = np.abs(r) <= 40
    out[small] = [float(digamma(mpc(0.25, x / 2)).real) for x in r[small]]
    z = 0.25 + 0.5j * r[~small]
    out[~small] = np.real(np.log(z) - 1 / (2 * z) - 1 / (12 * z**2) + 1 / (120 * z**4) - 1 / (252 * z**6))
    return out


def explicit_prediction(t, Gamma, dr=0.05, nmax=None):
    t = np.asarray(t, float)
    G = lambda u: Gamma / np.sqrt(2 * np.pi) * np.exp(-(Gamma**2) * u**2 / 2)
    gt = lambda u, tt: 0.5 * (G(u - tt) + G(u + tt))
    r = np.arange(0.0, 7 * Gamma, dr)
    rp = _re_psi(r)
    wr = np.exp(-(r / Gamma) ** 2 / 2) * rp
    wr[0] *= 0.5
    I = np.array([np.sum(np.cos(r * x) * wr) * dr / np.pi for x in t])        # (1/2π)∫_{−∞}^{∞} = (1/π)∫_0^∞
    pole = 2 * np.cosh(t / 2) * np.exp(1 / (8 * Gamma**2))
    nmax = nmax or int(np.exp(t.max() + 8 / Gamma * 6)) + 2
    lam = lam_table(nmax)
    n = np.arange(2, nmax + 1)
    prime_terms = np.array([np.sum(lam[2:] / np.sqrt(n) * gt(np.log(n), x)) for x in t])
    total = pole - gt(0.0, t) * np.log(np.pi) + I - 2 * prime_terms
    return 0.5 * total, 0.5 * (-2 * prime_terms), 0.5 * (pole + I - gt(0.0, t) * np.log(np.pi))
