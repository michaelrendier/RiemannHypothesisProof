# SPDX-License-Identifier: GPL-3.0-only
"""The two hemispheres as one oscillation, and the equator as its nodal line: parities, flux, the helix of nodal lines, the lobe.

Notation: ξ(s) entire, ξ(1−s) = ξ(s), ξ(s̄) = conj ξ(s).  R : (σ, t) ↦ (1−σ, t) is the north–south mirror at fixed height.
ESTABLISHED (exact, checked in tests):
  ξ∘R = conj ξ     so  Re ξ is EVEN across the equator and Im ξ is ODD: two polarisations with different boundary conditions at σ = ½
                    (Im ξ = 0: Dirichlet — a held point, a fret;  ∂σ Re ξ = 0: Neumann — a free point).
  Cauchy–Riemann:  Re ξ and Im ξ are not independent; ∂σ Re = ∂t Im, ∂t Re = −∂σ Im.
  'Motion' defined here (the flux):  J = Im( conj ξ · ∂ξ/∂t ) = |ξ|²·∂t arg ξ  = ½ ∂σ |ξ|²  (by CR).  J ≡ 0 on the equator (ξ is real there): nothing flows.
     J is odd across the equator: the two hemispheres carry opposite flux.
OURS-CALCULATED: the maps below.  RH ⇔ (equivalent standard reformulation) |ξ(σ+it)| increases with σ for σ > ½, i.e. J > 0 in the north:
  a negative-flux pocket is the footprint of an off-line zero (tested here by injecting a quadruple; NOT a search for one).
READING (not a theorem): in the chart  (latitude π(σ−½), longitude θ(t) mod 2π)  each rung (Re ξ = 0 through a zero, Im ξ = 0 through a saddle)
  is a MERIDIAN ARC, and consecutive rungs step in longitude by the lobe angle (mean π, distributed like the GUE zero spacing × π):
  a ladder wrapped around the polar axis.  The equator is the one great circle all rungs cross at right angles.
"""
from __future__ import annotations

import numpy as np
from mpmath import mp, mpc, mpf, siegelz, siegeltheta, quad

from . import nodal as N


# ── parities and flux ─────────────────────────────────────────────────────────────────────────
def parity_defect(sig_offsets, ts):
    """max over samples of  |Re ξ(½+x+it) − Re ξ(½−x+it)|  and  |Im ξ(½+x+it) + Im ξ(½−x+it)|  (relative): both ≈ 0 ⇒ even/odd."""
    worst_re = worst_im = 0.0
    for x in sig_offsets:
        for t in ts:
            a, b = N.xi(mpc(0.5 + x, t)), N.xi(mpc(0.5 - x, t))
            sc = abs(a) + 1e-300
            worst_re = max(worst_re, float(abs(a.real - b.real) / sc)); worst_im = max(worst_im, float(abs(a.imag + b.imag) / sc))
    return worst_re, worst_im


def flux_two_ways(sig, t, h=1e-5):
    """(½∂σ|ξ|², Im(conj ξ ∂tξ)) at s = σ + it by central differences: equal by Cauchy–Riemann."""
    f = lambda a, b: N.xi(mpc(a, b))
    dsig = (abs(f(sig + h, t)) ** 2 - abs(f(sig - h, t)) ** 2) / (2 * h) / 2
    dt = (f(sig, t + h) - f(sig, t - h)) / (2 * h)
    return float(dsig), float((f(sig, t).conjugate() * dt).imag)


def flux_grid(sig, t, xifun=None, h=1e-4):
    """J = ½ ∂σ|ξ|² on a σ × t grid (rows t)."""
    xf = xifun or N.xi
    out = np.empty((len(t), len(sig)))
    for i, tt in enumerate(t):
        for j, ss in enumerate(sig):
            out[i, j] = float((abs(xf(mpc(ss + h, tt))) ** 2 - abs(xf(mpc(ss - h, tt))) ** 2) / (2 * h) / 2)
    return out


def rel_flux_grid(sig, t, xifun=None, h=1e-4):
    """Scale-free flux ∂σ ln|ξ| = J/|ξ|² = Re(ξ′/ξ) on a σ × t grid (same sign as J; finite where |ξ| is tiny)."""
    xf = xifun or N.xi
    out = np.empty((len(t), len(sig)))
    for i, tt in enumerate(t):
        for j, ss in enumerate(sig):
            a, b = abs(xf(mpc(ss + h, tt))), abs(xf(mpc(ss - h, tt)))
            out[i, j] = float((mp.log(a) - mp.log(b)) / (2 * h))
    return out


def north_violations(J, sig):
    """Count of north-hemisphere (σ > ½) grid points with J ≤ 0 (a flux reversal)."""
    cols = np.where(np.asarray(sig) > 0.5 + 1e-9)[0]
    return int(np.sum(J[:, cols] <= 0))


def with_offline_quadruple(delta, gamma0):
    """ξ with the four extra zeros of an off-line pair at ½ ± δ ± iγ₀ multiplied in (symmetric: keeps ξ(1−s) = ξ(s) and conjugation)."""
    r0 = mpc(0.5 + delta, gamma0)
    roots = [r0, r0.conjugate(), 1 - r0, 1 - r0.conjugate()]

    def f(s):
        m = mpc(1)
        for r in roots:
            m *= (s - r)
        return N.xi(s) * m
    return f


# ── the helix of nodal lines on the sphere ────────────────────────────────────────────────────
def chart_xyz(sig, t):
    """Sphere chart: latitude π(σ − ½), longitude θ(t) mod 2π (the Hardy phase): returns an (n,3) array."""
    th = np.array([float(siegeltheta(mpf(x))) for x in np.atleast_1d(t)])
    lat = np.pi * (np.asarray(sig, float) - 0.5)
    return np.column_stack([np.cos(lat) * np.cos(th), np.cos(lat) * np.sin(th), np.sin(lat)])


def nodal_segments(sig, t, Z):
    """Polylines (σ, t) of {Re ξ = 0} and {Im ξ = 0} from a grid Z = N.grid(sig, t)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    out = {}
    for name, f in (("re", Z.real), ("im", Z.imag)):
        cs = ax.contour(sig, t, f, levels=[0])
        out[name] = [seg for seg in cs.allsegs[0] if len(seg) > 1]
    plt.close(fig)
    return out


# ── the lobe: what one zero-to-zero turn sweeps ───────────────────────────────────────────────
def lobe_angles(zeros):
    """θ(γ_{k+1}) − θ(γ_k) for consecutive zeros: the angle swept by the Hardy phase over one lobe.  Mean π (Riemann–von Mangoldt)."""
    th = np.array([float(siegeltheta(mpf(g))) for g in zeros])
    return np.diff(th)


def lobe_area(g0, g1):
    """Area swept by the rose r = Z(t), φ = −θ(t) between consecutive zeros: ½ ∫ Z² θ′ dt (the lobe's area)."""
    from .uns import theta_prime
    f = lambda t: float(siegelz(mpf(t))) ** 2 * theta_prime(t) / 2
    return float(quad(f, [g0, g1]))
