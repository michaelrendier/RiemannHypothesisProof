# SPDX-License-Identifier: GPL-3.0-only
"""Tests of the RHP machinery against the 0_RB (ValaQuenta) blocks that touch ζ.  0_RB is read as the stated authority; every
function here either REPRODUCES a 0_RB statement from independent code, or reports where a 0_RB statement and this code part ways.

Where ζ lives in 0_RB (located 2026-11-07; paths under ValaQuenta/modules/):
  h_rb_hat/maths.py        Σ_RB = Σ_p p^{−σ}[R̂_p⊗∂̂ + ∂̂†⊗B̂_p].  ζ is the COUPLING p^{−σ} (the Euler product); R̂†=B̂ is the
                           functional equation; Red = x·p (Berry–Keating, 'what is'), Blue = ½p²+℘ (Fermat–Weierstrass, 'what cannot be').
  archimedes_screw/maths.py  the working axis u = ln x; ψ(eᵘ) = eᵘ − Σ tones(γ) − ln 2π − ½ln(1−e^{−2u}); jump ln p at u = ln p.
  .clauderc_canonical_maths  @PRIME WAVELENGTH / FACTOR WAVELENGTH: λ_prime(n)=γ_{n+1}−γ_n (zero spacing), λ_factor(p)=p (the circle).
                           @RING_THEORY_SPINE: n falls ⇔ Z/(n) has zero divisors ⇔ gcd(a,n)>1;  order-2 'cepstrum' = von Mangoldt Λ.
"""
from __future__ import annotations

import importlib.util
import math
import os
import pathlib
import numpy as np

VQ = pathlib.Path(os.environ.get("VALAQUENTA", "/home/rendier/Projects/ThePlace/ValaQuenta"))


def _load(relpath, name):
    p = VQ / relpath
    if not p.exists():
        return None
    spec = importlib.util.spec_from_file_location(name, p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def load_screw():
    return _load("modules/archimedes_screw/maths.py", "vq_screw")


def load_hrb():
    return _load("modules/h_rb_hat/maths.py", "vq_hrb")


# T1 ── the Archimedes Screw's explicit formula reproduces ψ from the zeros (independent of this repo's code) ──
def screw_reconstruction(x_list=(30.5, 100.5, 300.5), Ks=(10, 50, 200)):
    A = load_screw()
    out = []
    for x in x_list:
        ex = A.chebyshev_psi_exact(x)
        for K in Ks:
            ap = A.chebyshev_psi_explicit(x, A.zeros_upto(K))
            out.append((x, K, ex, ap, abs(ex - ap)))
    return out


def jumps_at_ln_p(N=200):
    """The screw says ψ jumps by ln p at u = ln p.  Returns the maximum |jump − ln p| over prime powers ≤ N."""
    A = load_screw()
    worst = 0.0
    for n, jump in A.leaf_drops(float(N)):
        lam = A.von_mangoldt(n)
        worst = max(worst, abs(jump - lam))
    return worst


# T2 ── Σ_RB 'balance': does it select σ = ½?  (what the code actually computes) ──
def sigma_rb_balance_vs_sigma(x=1.0, p=2.0, sigmas=(0.0, 0.25, 0.5, 1.0, 2.0)):
    H = load_hrb()
    return [(s, H.sigma_rb_evaluate(s, x, p)["balance"], H.sigma_rb_evaluate(s, x, p)["on_critical_line"]) for s in sigmas]


def sigma_rb_balance_vs_xp(sigma=0.5, pairs=((1.0, 2.0), (1.0, 1.0), (1.0, 0.5), (0.5, 0.5))):
    H = load_hrb()
    return [(x, p, H.sigma_rb_evaluate(sigma, x, p)["balance"], H.sigma_rb_evaluate(sigma, x, p)["on_critical_line"]) for x, p in pairs]


# T3 ── Red/Blue dominance by hemisphere: forward (Red) vs mirror (Blue) wave, from the polarisation ratio r = |χ| n^{2σ−1} ──
def dominance(n, sigma, t):
    from .sphere import amplitude_ratio
    r = amplitude_ratio(n, sigma, t)
    return "mirror (Blue) dominates" if r > 1 else ("forward (Red) dominates" if r < 1 else "balanced")


def dominance_table(t=100.0, sigmas=(0.25, 0.75), ns=(1, 2, 3, 4)):
    cut = math.sqrt(t / (2 * math.pi))
    return cut, {(s, n): dominance(n, s, t) for s in sigmas for n in ns}


# T4 ── λ_prime and λ_factor meet at the resonance: at t = 2πn² the zero spacing is exactly half the n-th wheel's period ──
def wavelength_pair(n):
    t = 2 * math.pi * n * n
    lam_prime_smooth = 2 * math.pi / math.log(t / (2 * math.pi))      # Riemann–von Mangoldt mean spacing
    lam_wheel = 2 * math.pi / math.log(n)                              # period of the wheel n in t
    return t, lam_prime_smooth, lam_wheel, lam_prime_smooth / lam_wheel


def measured_spacing(n, half_width=12.0):
    """Mean spacing of the actual zeros within ±half_width of t = 2πn²."""
    from mpmath import zetazero, mp
    from .uns import theta
    t = 2 * math.pi * n * n
    k0 = int(theta(t) / math.pi)                                       # ≈ zero index at t (N(T) ≈ θ/π + 1)
    gs = []
    k = max(1, k0 - 25)
    while True:
        g = float(zetazero(k).imag)
        if g > t + half_width:
            break
        if g >= t - half_width:
            gs.append(g)
        k += 1
        if k > k0 + 60:
            break
    return float(np.mean(np.diff(gs))), len(gs)
