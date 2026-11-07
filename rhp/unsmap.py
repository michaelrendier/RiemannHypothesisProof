# SPDX-License-Identifier: GPL-3.0-only
"""Flat Cartesian → UNS: number-free mapping and artefact detection.  Nothing here mentions primes, factors or ζ.

For GenerationalLineage ('artefact detection and identification'): given a flat trace, say what shape made it.
Orientation: complex plane s = σ + it  →  Riemann sphere by the Cayley map w = (s−1)/s then inverse stereographic projection
(north pole w = 0, south pole w = ∞).  Aperture: the whole plane, with the point at infinity included.
ESTABLISHED (Möbius geometry):
  • every line and circle of the plane is a circle on the sphere (a straight line is a circle through the pole of the chart);
  • the map is conformal: angles survive, so a genuine cusp is still a cusp and a spiral of constant angle is still one (a loxodrome);
  • the line σ = ½ is the equator (|w| = 1): a GREAT circle; every other line σ = c is a small circle.
OURS-CALCULATED (this file): the detectors below, checked in tests/test_unsmap.py on synthetic shapes with known answers.
  1. cusp or reversal?  a GENUINE cusp stays a cusp after the lift (speed → 0 at the same parameter).  A reversal that comes from
     viewing a UNIFORM ROTATION edge-on (a 1-D back-and-forth) shows a constant analytic-signal amplitude and a linear phase:
     the Hilbert companion is the missing flat coordinate.
  2. runaway?   a plane curve that flies off to infinity and returns is a smooth passage over the pole: finite sphere curvature.
  3. spiral kind?  the radius–angle law fixes it: Archimedean (r linear in angle: constant pitch, the pancake coil), logarithmic
     (ln r linear in angle: a loxodrome on the sphere), Euler/Cornu (curvature linear in arc length).
"""
from __future__ import annotations

import numpy as np
from scipy.signal import hilbert


# ── the map ───────────────────────────────────────────────────────────────────────────────────
def cayley(s):
    s = np.asarray(s, dtype=complex)
    return (s - 1) / s


def stereo(w):
    """Inverse stereographic image of w on the unit sphere, north pole at w = 0.  w = ∞ → south pole (0,0,−1)."""
    w = np.asarray(w, dtype=complex)
    a2 = np.abs(w) ** 2
    out = np.empty(w.shape + (3,))
    inf = ~np.isfinite(a2)
    d = np.where(inf, 1.0, 1 + a2)
    out[..., 0] = np.where(inf, 0, 2 * w.real / d)
    out[..., 1] = np.where(inf, 0, 2 * w.imag / d)
    out[..., 2] = np.where(inf, -1, (1 - a2) / d)
    return out


def lift(s):
    """The plane curve s(t) as a curve on the sphere (Cayley then stereographic)."""
    return stereo(cayley(s))


def plane_fit(P):
    """Least-squares plane through 3-D points: (unit normal, offset d with n·x = d, rms residual).  A circle on the sphere lies in a plane;
    d = 0 means a GREAT circle."""
    c = P.mean(axis=0)
    u, s, vt = np.linalg.svd(P - c)
    n = vt[-1]
    d = float(n @ c)
    return n, d, float(np.sqrt(np.mean(((P - c) @ n) ** 2)))


# ── curvature: flat versus on the sphere ──────────────────────────────────────────────────────
def _derivs(P, t):
    d1 = np.gradient(P, t, axis=0, edge_order=2)
    d2 = np.gradient(d1, t, axis=0, edge_order=2)
    return d1, d2


def plane_curvature(z, t):
    """Signed curvature of the planar curve z(t)."""
    z1 = np.gradient(z, t, edge_order=2); z2 = np.gradient(z1, t, edge_order=2)
    return (z1.real * z2.imag - z1.imag * z2.real) / np.abs(z1) ** 3


def sphere_geodesic_curvature(P, t):
    """Geodesic curvature of a curve on the unit sphere: (γ′ × γ″)·γ / |γ′|³ (the sphere's normal at γ is γ itself)."""
    d1, d2 = _derivs(P, t)
    return np.einsum("ij,ij->i", np.cross(d1, d2), P) / np.linalg.norm(d1, axis=1) ** 3


def speed(P, t):
    return np.linalg.norm(np.gradient(P, t, axis=0, edge_order=2), axis=1)


# ── detector 1: cusp or reversal? ─────────────────────────────────────────────────────────────
def analytic_signal(x):
    """(amplitude, unwrapped phase) of a real trace x by the Hilbert transform: the hidden flat coordinate that completes a 1-D reversal
    into a rotation."""
    a = hilbert(np.asarray(x, float))
    return np.abs(a), np.unwrap(np.angle(a))


def edge_on_rotation_test(t, x, trim=0.1):
    """Is the back-and-forth trace x(t) a uniform rotation seen edge-on?  Returns dict(amp_cv, phase_linear_r2, reversals).
    Edge-on uniform rotation: amplitude constant (cv ≈ 0) and phase linear in t (r² ≈ 1).  A genuine cusp or a driven oscillator is not."""
    amp, ph = analytic_signal(x)
    k = int(len(t) * trim)
    sl = slice(k, len(t) - k)
    cv = float(np.std(amp[sl]) / np.mean(amp[sl]))
    p = np.polyfit(t[sl], ph[sl], 1)
    res = ph[sl] - np.polyval(p, t[sl])
    r2 = 1 - float(np.var(res) / np.var(ph[sl]))
    d = np.diff(np.sign(np.gradient(x, t)))
    return {"amp_cv": cv, "phase_linear_r2": r2, "reversals": int(np.sum(d != 0))}


def cusp_survives_lift(z, t, at):
    """For a planar curve with a suspected cusp at parameter `at`: (flat speed, sphere speed) there relative to their medians.
    Both → 0 together for a genuine cusp (conformal lift keeps it)."""
    P = lift(z)
    i = int(np.argmin(np.abs(t - at)))
    f = np.abs(np.gradient(z, t, edge_order=2)); g = speed(P, t)
    return float(f[i] / np.median(f)), float(g[i] / np.median(g))


# ── detector 2: runaway to infinity = a passage over the pole ────────────────────────────────
def pole_passage(w, t, big=50.0):
    """For a curve w(t) in the chart plane: (max |w|, min |w|, max |geodesic curvature| on the sphere) — a runaway with finite sphere
    curvature is smooth through the pole."""
    P = stereo(w)
    kg = sphere_geodesic_curvature(P, t)
    return float(np.max(np.abs(w))), float(np.min(np.abs(w))), float(np.max(np.abs(kg)))


# ── detector 3: what kind of spiral? ──────────────────────────────────────────────────────────
def spiral_laws(z, centre=0j):
    """Fit the three laws to a spiral z (points in order): returns r² of  r vs φ (Archimedean),  ln r vs φ (logarithmic),
    κ vs arc length (Euler/Cornu).  φ is the unwrapped polar angle about `centre`."""
    z = np.asarray(z) - centre
    z = z[:-1] if np.abs(z[-1]) < 1e-9 else z
    r = np.abs(z); phi = np.unwrap(np.angle(z))
    def r2(x, y):
        p = np.polyfit(x, y, 1)
        return 1 - float(np.var(y - np.polyval(p, x)) / np.var(y))
    s = np.concatenate([[0], np.cumsum(np.abs(np.diff(z)))])
    t = np.arange(len(z), dtype=float)
    kk = plane_curvature(z, t)                              # geometric curvature (parametrisation-free)
    return {"archimedean": r2(phi, r), "logarithmic": r2(phi, np.log(r)), "euler": r2(s[2:-2], kk[2:-2])}
