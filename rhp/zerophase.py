# SPDX-License-Identifier: GPL-3.0-only
"""Zero-phase windings against Archimedean (Ulam / Sacks / polygonal) windings: which arrangement do the primes follow?   :CALCULATED unless marked ESTABLISHED.

Notation (the code names the variation, never leaves it implicit):
  winding   how the integer n is turned into an angle
    'zero'   θ_k(n) = γ_k ln n           a zero-tone x^ρ = e^{(½+iγ)u} IS the logarithmic spiral r = √x, θ = γ ln x (ESTABLISHED, an identity); pitch angle arctan(1/2γ)
    'arch'   θ_h(n) = 2π h √n            Archimedean in √n.  h = 1: Sacks / Ulam (perfect squares on one ray).  The m-gon spiral is the same family with h = √(2/(m−2)):
                                         polygonal numbers P_m(k) = ((m−2)k² − (m−4)k)/2 lie on the ray θ = 2πk at n = (m−2)k²/2.  So 'which polygon' is one axis, h.
  lift      how (|z|, θ) becomes a point of the sphere
    'stereo' inverse stereographic projection with the equator at |z| = 1:  Z = (|z|² − 1)/(|z|² + 1) = tanh(u/2) for |z| = √x = e^{u/2}·…   conformal; the log spiral e^{ρu} lifts to a LOXODROME
    'linear' colatitude linear in ln n (the lift actually coded in AbrikosovTree/render/lattice_planes.py `_lift`)           not conformal; a log spiral is no longer a rhumb line
The critical line is the equator: for z = x^{ρ−½} = e^{(σ−½+iγ)u}, |z| = 1 exactly when σ = ½.
Spectral scans use ONE statistic for both windings:  A(θ) = | Σ_{n≤x} (Λ(n) − 1) n^{−½} e^{iθ(n)} | / ln x   (Landau–Gonek: the −1 removes the smooth part; at a zero ordinate the primes align and A → ~1)
"""
from __future__ import annotations
import math
import numpy as np


# ---------------------------------------------------------------- lifts
def lift_stereo(z):
    """Unit-sphere point (x, y, z) of the complex number(s) z, equator at |z| = 1, z → ∞ at the north pole."""
    z = np.asarray(z, complex); r2 = np.abs(z) ** 2
    return np.stack([2 * z.real / (r2 + 1), 2 * z.imag / (r2 + 1), (r2 - 1) / (r2 + 1)], axis=-1)


def lift_linear(u, theta, umax, margin=0.06):
    """colatitude linear in u = ln n (π(1 − margin) at u = 0, π·margin at u = umax), longitude θ."""
    u = np.asarray(u, float); colat = math.pi * ((1 - margin) - (1 - 2 * margin) * u / umax)
    return np.stack([np.sin(colat) * np.cos(theta), np.sin(colat) * np.sin(theta), np.cos(colat)], axis=-1)


def bearing(P):
    """Angle of the path's direction from the local meridian (north), radians, along a path P (m × 3 on the unit sphere): constant on a rhumb line."""
    X = 0.5 * (P[1:] + P[:-1]); T = P[1:] - P[:-1]; ez = np.array([0, 0, 1.0])
    N = ez - (X @ ez)[:, None] * X; N /= np.linalg.norm(N, axis=1)[:, None]; E = np.cross(ez, X); E /= np.linalg.norm(E, axis=1)[:, None]
    return np.arctan2((T * E).sum(1), (T * N).sum(1))


def pitch_deg(gamma):
    """pitch angle of the zero-tone's log spiral: arctan(1/2γ), degrees."""
    return math.degrees(math.atan(1.0 / (2 * gamma)))


def zero_spiral(u, gamma, sigma=0.5):
    """z(u) = e^{(σ + iγ)u}  (σ = ½: the zero-tone x^ρ itself, |z| = √x)."""
    return np.exp((sigma + 1j * gamma) * np.asarray(u, float))


# ---------------------------------------------------------------- the common statistic
def weights(X):
    """w[n] = (Λ(n) − 1)/√n for n = 1 … X (w[0] = 0), and ln X."""
    lam = np.zeros(X + 1); sv = np.ones(X + 1, bool); sv[:2] = False
    for p in range(2, int(X ** .5) + 1):
        if sv[p]: sv[p * p::p] = False
    for p in np.nonzero(sv)[0]:
        q = int(p)
        while q <= X: lam[q] = math.log(p); q *= int(p)
    n = np.arange(X + 1.0); n[0] = 1; w = (lam - 1.0) / np.sqrt(n); w[0] = 0.0
    return w, math.log(X)


def scan_zero(ts, w, lx, chunk=40):
    """A(t) with θ = t ln n."""
    n = np.arange(len(w), dtype=float); n[0] = 1; ln = np.log(n); out = np.empty(len(ts))
    for a in range(0, len(ts), chunk):
        out[a:a + chunk] = np.abs((np.exp(-1j * np.outer(ts[a:a + chunk], ln)) * w).sum(axis=1))
    return out / lx


def scan_arch(hs, w, lx, chunk=40):
    """A(h) with θ = 2π h √n."""
    n = np.arange(len(w), dtype=float); sq = np.sqrt(n); out = np.empty(len(hs))
    for a in range(0, len(hs), chunk):
        out[a:a + chunk] = np.abs((np.exp(2j * np.pi * np.outer(hs[a:a + chunk], sq)) * w).sum(axis=1))
    return out / lx


def polygon_h(m):
    """the Archimedean winding of the m-gon spiral: h = √(2/(m−2))."""
    return math.sqrt(2.0 / (m - 2))


def peaks(a, x, thr):
    """local maxima of a above thr: (positions, heights)."""
    i = np.nonzero((a[1:-1] > a[:-2]) & (a[1:-1] >= a[2:]) & (a[1:-1] > thr))[0] + 1
    return x[i], a[i]
