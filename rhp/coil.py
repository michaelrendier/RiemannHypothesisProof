# SPDX-License-Identifier: GPL-3.0-only
"""Spirals, the pancake coil, and the conical helices of the explicit formula.

POE pancake coil (ContextPlease/claude/hist_prime/POE/hardware_primer.md): 38 mm outer, 10 mm inner diameter, pitch 0.8 mm, 17 turns;
taps T1 / T3 / T5 / T_full = 1 / 3 / 5 / 17 inner turns.  It is an ARCHIMEDEAN spiral (constant radial pitch).

ESTABLISHED: θ(γ_k) = π(k − 1 − S) with S the bounded zero-count fluctuation (Riemann–von Mangoldt), so the Hardy phase advances π per zero;
  the explicit formula ψ(x) = x − Σ_ρ x^ρ/ρ − …  makes every zero a CONICAL HELIX in (u = ln x): radius ∝ e^{u/2}, pitch 2π/γ in u.
OURS-CALCULATED (this file, tests in tests/test_unsmap.py):
  • radius = zero index k (= N(T)) against the Hardy angle θ is an Archimedean spiral up to S: a pancake coil with ONE TURN PER TWO ZEROS;
    radius = height t is NOT (its pitch 4π/ln(t/2π) tightens);
  • the partial sums Σ_{n≤N} n^{−s} wind about the value ζ(s) along a LOGARITHMIC spiral (a loxodrome on the sphere): S(N) − ζ(s) ≈ N^{1−s}/(1−s)
    (Euler–Maclaurin), so d ln r/dφ = −(1−σ)/t.  σ > 1: inward (converges onto ζ); σ = 1: a CIRCLE (marginal); σ < 1: outward, r ∝ N^{1−σ}.
    The critical line is outward at 1/(2t); ζ(s) is the CENTRE of the spiral at every latitude, never a point it is pulled toward on the line;
  • in the explicit formula the primes are where all the conical helices come into phase (the jump ln p at u = ln p): the one spiral with no free choice.
"""
from __future__ import annotations

import math
import numpy as np
from mpmath import mpf, mpc, siegeltheta, zeta

PANCAKE = dict(d_out=38.0, d_in=10.0, pitch=0.8, turns=17)
TAPS = {"T1": 1, "T3": 3, "T5": 5, "T_full": 17}


def pancake(n_pts=4000, d_in=PANCAKE["d_in"], pitch=PANCAKE["pitch"], turns=PANCAKE["turns"]):
    """The coil as a curve: r = d_in/2 + pitch·φ/2π, φ ∈ [0, 2π·turns].  Returns (x, y, phi, r) in mm."""
    phi = np.linspace(0, 2 * np.pi * turns, n_pts)
    r = d_in / 2 + pitch * phi / (2 * np.pi)
    return r * np.cos(phi), r * np.sin(phi), phi, r


def zeros_on_coil(zeros, d_in=PANCAKE["d_in"], pitch=PANCAKE["pitch"]):
    """Place zero k (k = 1, 2, …) on a coil with ONE TURN PER TWO ZEROS: radius d_in/2 + pitch·(k−1)/2, angle θ(γ_k).  Returns (x, y, theta, r)."""
    k = np.arange(1, len(zeros) + 1)
    th = np.array([float(siegeltheta(mpf(g))) for g in zeros])
    r = d_in / 2 + pitch * (k - 1) / 2
    return r * np.cos(th), r * np.sin(th), th, r


def radius_law(zeros):
    """r² of a straight line for  (radius vs Hardy angle):  index k vs θ(γ_k)   — Archimedean (pancake) — and  height γ_k vs θ(γ_k)."""
    th = np.array([float(siegeltheta(mpf(g))) for g in zeros])
    k = np.arange(1, len(zeros) + 1, dtype=float)
    def r2(x, y):
        p = np.polyfit(x, y, 1)
        return 1 - float(np.var(y - np.polyval(p, x)) / np.var(y))
    return {"index_vs_angle": r2(th, k), "height_vs_angle": r2(th, np.asarray(zeros, float))}


def hardy_pitch_t(t):
    """Height per 2π turn of the Hardy angle: 2π/θ′(t) ≈ 4π/ln(t/2π) — what a coil wound on the height would need to be (a tightening pitch)."""
    return 4 * np.pi / np.log(np.asarray(t, float) / (2 * np.pi))


def tap_zeros(zeros):
    """Zeros at the end of each tap's turns (two zeros per turn): {tap: (zero index, γ)}."""
    return {name: (2 * n, float(zeros[2 * n - 1])) for name, n in TAPS.items()}


def partial_sum_path(t, N0, N1, sigma=0.5):
    """Partial sums S(N) = Σ_{n ≤ N} n^{−s} (s = σ + it) for N = N0..N1, relative to ζ(s): returns (N, S − ζ)."""
    s = complex(sigma, t)
    n = np.arange(1, N1 + 1)
    S = np.cumsum(n.astype(float) ** (-s))
    z = complex(zeta(mpc(sigma, t)))
    return n[N0 - 1:], S[N0 - 1:] - z


def tail_log_slope(t, N0=400, N1=6000, sigma=0.5):
    """Measured d ln r / dφ of the partial-sum spiral about ζ(s), and the prediction −(1−σ)/t."""
    n, d = partial_sum_path(t, N0, N1, sigma)
    phi = np.unwrap(np.angle(d)); r = np.abs(d)
    return float(np.polyfit(phi, np.log(r), 1)[0]), -(1 - sigma) / t


def explicit_path(u, zeros):
    """Q(u) = Σ_k e^{iγ_k u}/ρ_k  (ρ_k = ½ + iγ_k): the complex trace of the conical helices with the e^{u/2} cone divided out.
    E(u) = 2 Re Q(u) is the explicit formula's (x − ψ(x) + …)/√x; Im Q is the hidden quadrature of the 1-D trace."""
    u = np.asarray(u, float)
    rho = 0.5 + 1j * np.asarray(zeros, float)
    return np.exp(1j * np.outer(u, np.asarray(zeros, float))) @ (1 / rho)


def helix(gamma, u):
    """Conical helix of one zero:  (e^{u/2}cos γu, e^{u/2}sin γu, u) — a 3-D curve whose side view is the growing wave e^{u/2}cos γu."""
    u = np.asarray(u, float)
    r = np.exp(u / 2)
    return np.column_stack([r * np.cos(gamma * u), r * np.sin(gamma * u), u])


# ── not flat: the same objects as 3-D curves (UNS cylindrical chart) ─────────────────────────
def zeta_cylinder(ts, sigma=0.5):
    """ζ(σ + it) in the cylindrical UNS chart: (ρ cos φ, ρ sin φ, t) with φ = arg ζ and ρ = (2/π)·arctan|ζ| ∈ [0,1) (the value plane compactified).
    The t-axis is the set ζ = 0: at σ = ½ the pathway THREADS the axis at each zero (the Hardy rose seen with its height), off the line it circles the axis."""
    z = np.array([complex(zeta(mpc(sigma, float(t)))) for t in ts])
    rho = (2 / np.pi) * np.arctan(np.abs(z))
    return np.column_stack([rho * np.cos(np.angle(z)), rho * np.sin(np.angle(z)), np.asarray(ts, float)]), z


def axis_distance(ts, sigmas=(0.5, 0.6, 0.75, 1.0)):
    """min over t of |ζ(σ+it)| for each latitude: 0 exactly at the zeros on the line, and bounded away from 0 elsewhere in the window."""
    return {s: float(min(abs(complex(zeta(mpc(s, float(t))))) for t in ts)) for s in sigmas}


def partial_sum_cone(t, sigma=0.5, N1=4000, N0=2):
    """The partial-sum spiral about ζ(s) lifted by height u = ln N: (Re, Im of S(N) − ζ(s), ln N).  Radius ∝ N^{1−σ} = e^{(1−σ)u}:
    a CONE (flares) for σ < 1, a CYLINDER at σ = 1, a FUNNEL for σ > 1; the critical line flares as e^{u/2}, the same flare as the explicit formula."""
    n, d = partial_sum_path(t, N0, N1, sigma)
    return np.column_stack([d.real, d.imag, np.log(n)])


def explicit_cone(zeros, u):
    """The prime side as a conical spiral: (e^{u/2} Re Q, e^{u/2} Im Q, u) with Q = Σ_k e^{iγ_k u}/ρ_k.  Its radial excursions at u = ln p are the primes."""
    Q = explicit_path(u, zeros)
    r = np.exp(np.asarray(u) / 2)
    return np.column_stack([r * Q.real, r * Q.imag, np.asarray(u)])
