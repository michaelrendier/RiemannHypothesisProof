# SPDX-License-Identifier: GPL-3.0-only
"""Chebyshev bias from the zeros of L(s,χ₋₄): the Fermat (two-squares) side.   :CALCULATED
D(x) = Σ_{n≤x} χ₋₄(n)Λ(n) = ψ(x;4,1) − ψ(x;4,3) = −Σ_ρ x^ρ/ρ − c(χ) + Σ_k x^{−(2k−1)}/(2k−1)   (χ odd, trivial zeros at −1,−3,…; ESTABLISHED explicit formula)
Primes ≡ 1 mod 4 are Fermat's sums of two squares; ≡ 3 mod 4 are the inert primes.  D<0 means the inert primes lead (the bias).
c(χ) = L'/L(0,χ) is computed, not fitted."""
from __future__ import annotations
import numpy as np, mpmath as mp
CHI = [0, 1, 0, -1]

def Lam(t, dps=20):
    """completed Λ(½+it) = (π/4)^{-(s+1)/2} Γ((s+1)/2) L(s,χ₋₄), real on the line (root number 1)."""
    mp.mp.dps = dps; s = mp.mpf(.5) + 1j * t
    return ((mp.pi / 4) ** (-(s + 1) / 2) * mp.gamma((s + 1) / 2) * mp.dirichlet(s, CHI)).real

def zeros(tmax, step=0.1):
    out = []; ts = np.arange(1.0, tmax, step); v = [Lam(t) for t in ts]
    for a, b, fa, fb in zip(ts[:-1], ts[1:], v[:-1], v[1:]):
        if fa * fb < 0: out.append(float(mp.findroot(lambda t: Lam(t), (a, b), solver='anderson', tol=1e-12)))
    return np.array(out)

def c_chi():
    mp.mp.dps = 25; return float(mp.diff(lambda s: mp.dirichlet(s, CHI), 0) / mp.dirichlet(0, CHI))

def D_exact(X):
    """D(x) for x = 1…X (array index x)."""
    lam = np.zeros(X + 1); sv = np.ones(X + 1, bool); sv[:2] = False
    for p in range(2, int(X ** .5) + 1):
        if sv[p]: sv[p * p::p] = False
    for p in np.nonzero(sv)[0]:
        q = int(p)
        while q <= X: lam[q] = np.log(p) * (1 if q % 4 == 1 else -1 if q % 4 == 3 else 0); q *= int(p)
    return np.cumsum(lam)

def D_zeros(x, g, K=None, c=None):
    """−2Re Σ_{γ>0} x^ρ/ρ − c + Σ x^{-(2k-1)}/(2k-1)  (zeros symmetric about the real axis since χ is real)."""
    g = g[:K] if K else g; x = np.asarray(x, float); c = c_chi() if c is None else c
    rho = 0.5 + 1j * g; tones = -2 * np.real((np.exp(np.outer(np.log(x), rho)) / rho).sum(1))
    triv = 0.5 * np.log((1 + 1 / x) / (1 - 1 / x + 1e-300)) if False else np.arctanh(1 / x)
    return tones - c + triv
