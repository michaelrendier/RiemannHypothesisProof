# SPDX-License-Identifier: GPL-3.0-only
"""The ξ′/ξ flow on the strip, and the argument-principle accounting (Gauss's law for it).

Complex potential Φ = log ξ(s).  L(s) = ξ′/ξ(s) = 1/s + 1/(s−1) − ½ ln π + ½ ψ(s/2) + ζ′/ζ(s).
Velocity (u, v) = (Re L, −Im L), u along σ, v along t.  Zeros of ξ are sources + vortices.
Symmetry (ξ(s)=ξ(1−s), ξ(s̄)=conj ξ(s)):  L(1−s̄) = −conj L(s)  ⇒  u odd, v even under σ ↦ 1−σ.
Labels: ESTABLISHED (definitions); the grid checks are :CALCULATED (floating point, not a proof).
"""
from __future__ import annotations
import numpy as np
from mpmath import mp, mpc, mpf, zeta, digamma, log, pi, arg as marg

from . import bigpath

mp.dps = 15


def L(s):
    s = mpc(s)
    return complex(1 / s + 1 / (s - 1) - log(pi) / 2 + digamma(s / 2) / 2 + zeta(s, derivative=1) / zeta(s))


def flow_grid(sig, ts):
    """u = Re L, v = −Im L on a (σ × t) grid. Returns arrays shaped (len(ts), len(sig))."""
    U = np.zeros((len(ts), len(sig)))
    V = np.zeros_like(U)
    for i, t in enumerate(ts):
        for j, s in enumerate(sig):
            w = L(mpc(s, t))
            U[i, j], V[i, j] = w.real, -w.imag
    return U, V


def symmetry_residuals(pts):
    """max |u(s*) + u(s)| and |v(s*) − v(s)| over sample points, s* = 1 − s̄ (should be ~1e-13)."""
    du = dv = 0.0
    for s in pts:
        a, b = L(s), L(1 - np.conj(s))
        du = max(du, abs(b.real + a.real))
        dv = max(dv, abs(-b.imag - (-a.imag)))
    return du, dv


def S_of_T(T, step=0.02):
    """S(T) = (1/π) arg ζ(½+iT) by continuous variation along 2 → 2+iT → ½+iT.

    On the vertical leg arg stays in (−π/2, π/2) (|ζ(2+it) − 1| < 1), so the principal value starts the track;
    the horizontal leg is followed by unwrapping small steps. Returns (S, max |Δarg| per step).
    """
    a0 = float(marg(zeta(mpc(2, T))))
    sig = np.arange(2.0, 0.5 - step / 2, -step)
    ang = np.array([float(marg(zeta(mpc(s, T)))) for s in sig])
    ang[0] = a0
    d = np.diff(ang)
    d = (d + np.pi) % (2 * np.pi) - np.pi
    total = a0 + d.sum()
    return total / np.pi, float(np.max(np.abs(d)))


def argument_accounting(Ts, zeros_all):
    """For each height T: N from the argument principle (θ/π + 1 + S) against the number of zeros on the line below T."""
    rows = []
    for T in Ts:
        S, dmax = S_of_T(T)
        n_arg = bigpath.theta_fast([T])[0] / np.pi + 1 + S
        n_line = int(np.sum(zeros_all < T))
        rows.append(dict(T=float(T), S=S, N_arg=float(n_arg), N_line=n_line, gap=float(n_arg - n_line), step_max=dmax))
    return rows
