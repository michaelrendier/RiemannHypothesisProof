# SPDX-License-Identifier: GPL-3.0-only
"""The nodal map: where Re ξ = 0 and Im ξ = 0 on the strip, both hemispheres.

ξ(s) = ½ s(s−1) π^{−s/2} Γ(s/2) ζ(s) is entire (the poles of Λ are cancelled), ξ(s) = ξ(1−s), ξ(s̄) = conj ξ(s).
ESTABLISHED: on the critical line ξ(½+it) = F(t) is real, so Im ξ = 0 on the whole line; ∂ξ/∂σ = −i·F′(t) there (Cauchy–Riemann),
  so Im ξ is odd and Re ξ even across the equator, and a simple zero has Re ξ = 0 along the HORIZONTAL line t = γ locally.
OURS-CALCULATED (this file): the grid, the zero-crossing of Re ξ on the line (the 'rungs') and the saddles of the landscape.
"""
from __future__ import annotations

import numpy as np
from mpmath import mp, mpf, mpc, gamma, zeta, pi, diff

mp.dps = 15


def xi(s):
    with mp.workdps(20):
        s = mpc(s)
        if abs(s) < 1e-12 or abs(s - 1) < 1e-12:
            return mpc(0.5)
        return 0.5 * s * (s - 1) * pi ** (-s / 2) * gamma(s / 2) * zeta(s)


def grid(sig, t):
    """ξ on a σ × t grid (rows t, columns σ) as a complex numpy array."""
    out = np.empty((len(t), len(sig)), complex)
    for i, tt in enumerate(t):
        for j, ss in enumerate(sig):
            out[i, j] = complex(xi(mpc(ss, tt)))
    return out


def F(t):
    """The heartbeat F(t) = ξ(½ + it), real."""
    return float(xi(mpc(0.5, t)).real)


def rungs(tmax, n=4000):
    """Heights where Re ξ = 0 meets the equator = sign changes of F on (0, tmax]: the zeros seen as crossings of the 'rungs'."""
    t = np.linspace(0.5, tmax, n)
    f = np.array([F(x) for x in t])
    idx = np.flatnonzero(np.sign(f[:-1]) * np.sign(f[1:]) < 0)
    return [float(t[i] - f[i] * (t[i + 1] - t[i]) / (f[i + 1] - f[i])) for i in idx]


def heartbeat_extrema(tmax, n=4000):
    """Heights where F′(t) = 0: saddles of the analytic landscape ON the equator (ξ ≈ F − ½F″z² there)."""
    t = np.linspace(0.5, tmax, n)
    f = np.array([F(x) for x in t])
    d = np.diff(f)
    idx = np.flatnonzero(np.sign(d[:-1]) * np.sign(d[1:]) < 0)
    return [float(t[i + 1]) for i in idx]


def saddle_hessian(t0, h=1e-3):
    """At a heartbeat extremum t0: the Hessian of Re ξ in (σ, t) = [[−F″·(−1)…]] computed by finite differences. Returns eigenvalues."""
    s0 = mpc(0.5, t0)
    f = lambda a, b: float(xi(mpc(0.5 + a, t0 + b)).real)
    fss = (f(h, 0) - 2 * f(0, 0) + f(-h, 0)) / h ** 2
    ftt = (f(0, h) - 2 * f(0, 0) + f(0, -h)) / h ** 2
    fst = (f(h, h) - f(h, -h) - f(-h, h) + f(-h, -h)) / (4 * h * h)
    return np.linalg.eigvalsh(np.array([[fss, fst], [fst, ftt]]))


def offline_crossings(sig, t, Z):
    """Grid cells OFF the equator in which both Re ξ and Im ξ change sign (≈ points where Re ξ = Im ξ = 0, i.e. a zero off the line).
    Z is N.grid(sig, t).  Cells within one column of σ = ½ are skipped (Im ξ ≡ 0 there).  Returns a list of (σ, t) cell corners."""
    out = []
    re, im = Z.real, Z.imag
    for i in range(len(t) - 1):
        for j in range(len(sig) - 1):
            if abs(sig[j] - 0.5) < 1e-9 or abs(sig[j + 1] - 0.5) < 1e-9 or (sig[j] < 0.5 < sig[j + 1]):
                continue
            cr = re[i:i + 2, j:j + 2]; ci = im[i:i + 2, j:j + 2]
            if cr.min() < 0 < cr.max() and ci.min() < 0 < ci.max():
                out.append((float(sig[j]), float(t[i])))
    return out
