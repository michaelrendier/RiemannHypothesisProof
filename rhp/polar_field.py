# SPDX-License-Identifier: GPL-3.0-only
"""The April 2026 entropy/inertia polar field S(r, θ), I(r, θ) and its ceilings, ported from the claude.ai widget `entropy_inertia_polar_field` (2026-04-10) as recorded in
ContextPlease/claude/hist_prime/RiemannHypothesisProof/PRIMER_2026-11-12_CEILING_DESCENT_FROM_CLAUDE_AI.md §3, §5.2.   :CALCULATED (float64, vectorised).

    L(r) = log_b(r + 1)
    S(r, θ) = clamp_[0, c]( c (1 − e^{−α L}) + 0.08 sin(6.28 L + θ) e^{−0.3 r} )          entropy: rises from 0, bounded by the speed-of-causality fraction c
    I(r, θ) = clamp_[0, 1]( 1/(1 + e^{−(L − depth/2)}) + 0.12 cos(6.28 L − θ + π/2) e^{−0.2 r} )     inertia: logistic in L with a phase-conjugate ripple
The constants 0.08, 6.28, 0.3, 0.12, 0.2 are the widget's hand-chosen values.  d* in the widget is the radius of closest approach min |S − I| over 0.1 < r < 3, all θ; the label 'd* ≈ 0.246' was a fixed string.
"""
from __future__ import annotations
import math
import numpy as np


def _L(r, base):
    r = np.asarray(r, float); return np.log(r + 1.0) if base == math.e else np.log(r + 1.0) / math.log(base)


def entropy(r, th, alpha=0.137, base=math.e, cfrac=1.0):
    r = np.asarray(r, float); L = _L(r, base)
    v = cfrac * (1 - np.exp(-alpha * L)) + 0.08 * np.sin(6.28 * L + th) * np.exp(-0.3 * r)
    return np.where(r < 0.001, 0.0, np.clip(v, 0.0, cfrac))


def inertia(r, th, depth=1.57, base=math.e):
    r = np.asarray(r, float); L = _L(r, base)
    v = 1 / (1 + np.exp(-(L - depth * 0.5))) + 0.12 * np.cos(6.28 * L - th + math.pi / 2) * np.exp(-0.2 * r)
    return np.where(r < 0.001, 0.0, np.clip(v, 0.0, 1.0))


def closest_approach(alpha=0.137, depth=1.57, cfrac=1.0, base=math.e, r_lo=0.1, r_hi=3.0, n_r=3000, n_th=720, refine=2):
    """min |S − I| over r ∈ (r_lo, r_hi), θ ∈ [−π, π): coarse grid then `refine` nested zooms.  Returns (min, r, θ)."""
    rs = np.linspace(r_lo, r_hi, n_r)[1:-1]; ths = np.linspace(-math.pi, math.pi, n_th, endpoint=False)
    for _ in range(refine + 1):
        R, T = np.meshgrid(rs, ths); D = np.abs(entropy(R, T, alpha, base, cfrac) - inertia(R, T, depth, base)); i = np.unravel_index(D.argmin(), D.shape)
        best = (float(D[i]), float(R[i]), float(T[i]))
        dr = (rs[1] - rs[0]) * 4; dt = (ths[1] - ths[0]) * 4
        rs = np.linspace(max(r_lo, best[1] - dr), min(r_hi, best[1] + dr), 801); ths = np.linspace(best[2] - dt, best[2] + dt, 801)
    return best


def base_crossings(alpha=0.137, depth=1.57, base=math.e, rmax=50.0, n=200001):
    """Sign changes of S0 − I0 (no oscillation) on (0, rmax], and min |S0 − I0|."""
    r = np.linspace(1e-6, rmax, n); L = _L(r, base); d = (1 - np.exp(-alpha * L)) - 1 / (1 + np.exp(-(L - depth * 0.5)))
    idx = np.nonzero(d[:-1] * d[1:] < 0)[0]; return r[idx], float(np.abs(d).min())


def one_d_forms(d=0.246, alpha=1 / 137.036, tau=math.pi / 2):
    """The 1-D forms of RH proof v3 §3.6: S(d) = 1 − e^{−α ln(d+1)}, I(d) = 1/(1 + e^{−(ln(d+1) − τ/2)})."""
    L = math.log(d + 1); return 1 - math.exp(-alpha * L), 1 / (1 + math.exp(-(L - tau / 2)))


def gaussian_curvature(F, x, y, h=1e-4):
    """K = (F_xx F_yy − F_xy²)/(1 + F_x² + F_y²)² of the graph z = F(x, y), central differences (vectorised)."""
    f = F; fx = (f(x + h, y) - f(x - h, y)) / (2 * h); fy = (f(x, y + h) - f(x, y - h)) / (2 * h)
    fxx = (f(x + h, y) - 2 * f(x, y) + f(x - h, y)) / h ** 2; fyy = (f(x, y + h) - 2 * f(x, y) + f(x, y - h)) / h ** 2
    fxy = (f(x + h, y + h) - f(x + h, y - h) - f(x - h, y + h) + f(x - h, y - h)) / (4 * h * h)
    return (fxx * fyy - fxy ** 2) / (1 + fx ** 2 + fy ** 2) ** 2


def surfaces(alpha=0.137, depth=1.57):
    """The three candidate surfaces whose curvature 'K = 0 at d*' could mean (the record does not say): z = S − I, z = S, z = I (as functions of Cartesian x, y)."""
    def pol(x, y): return np.hypot(x, y), np.arctan2(y, x)
    def S(x, y): r, t = pol(x, y); return entropy(r, t, alpha)
    def I(x, y): r, t = pol(x, y); return inertia(r, t, depth)
    return {"S-I": lambda x, y: S(x, y) - I(x, y), "S": S, "I": I}


# ---------------------------------------------------------------- the ceilings
def planck_temperature():
    """T_Planck = √(ħc⁵/G)/k_B from ħ (exact), c (exact), G (CODATA 2018), k_B (exact)."""
    import mpmath as mp
    mp.mp.dps = 30; c = mp.mpf(299792458); hbar = mp.mpf("1.054571817e-34"); G = mp.mpf("6.67430e-11"); kB = mp.mpf("1.380649e-23")
    TP = mp.sqrt(hbar * c ** 5 / G) / kB; lP = mp.sqrt(hbar * G / c ** 3)
    return dict(T_P=TP, l_P=lP, bekenstein_radius=hbar * c / (2 * mp.pi * kB * TP), l_P_over_2pi=lP / (2 * mp.pi))


def hagedorn_ceiling_kelvin(quadrillion_F=140.0):
    """The 'thermal information ceiling' as the code states it: 140 quadrillion °F × 5/9 (conjecture_proof.py: 140e15 °F; event_horizon_sim.py: OMEGA_K = 77.8e15 K 'from 140Q °F × 5/9').  The 5/9 scale only (the offset is negligible at 10¹⁷)."""
    return quadrillion_F * 1e15 * 5.0 / 9.0
