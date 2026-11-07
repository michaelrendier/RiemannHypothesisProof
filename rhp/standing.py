# SPDX-License-Identifier: GPL-3.0-only
"""Checks of an externally supplied 'standing wave' summation (Gemini, pasted 2026-11-07) against this repo's machinery.

Claims tested (the narrative is not tested; only statements with a mathematical content):
  C1  J(z) = i/z̄ acts as (r,θ) → (1/r, θ+π/2); four applications return the start.                        TRUE (J² = −z, J⁴ = id)
  C2  'the fixed set of J is r = 1'.                                                                       FALSE: J has no fixed point; r = 1 is INVARIANT
      (rotated a quarter turn).  J = (quarter turn) ∘ (w ↦ 1/w̄): the mirror of this repo (fixed pointwise on |w| = 1) composed with a rotation.
  C3  'J governs ζ / ξ'.                                                                                   NOT A SYMMETRY of ξ (tested): only the mirror part is.
  C4  Ψ = R(r)Θ(θ), R = A J₀(k ln r) + B Y₀(k ln r), 'node anchored at r = 1'.                             FALSE: J₀(0) = 1 (antinode), Y₀ diverges at r = 1.
  C5  ∫∫ |Ψ|² r dr dθ over the unit disc = π.                                                              FALSE: the integral given is ∫∫ r dr dθ = π (the disc's area,
      independent of Ψ); for Ψ = cos(k ln r)cos(mθ+π/4) the energy is (π/4)(1 + 1/(1+k²)).
  C6  Ψ is oscillatory in ln r (Helmholtz-like) with k = γ.                                                 NOT what ξ does: ξ is harmonic; across the equator the profile is
      cosh/sinh with transverse rate θ′(t) — oscillation along t, growth across σ (Cauchy–Riemann).   [fit below]
ESTABLISHED background: the real involution of the modular side is the Fricke / Atkin–Lehner map τ ↦ −1/(Nτ), which gives the functional equation of
L(E,s) with root number ε = ±1; the hyperbolic area of the PSL(2,ℤ) fundamental domain is π/3.
"""
from __future__ import annotations

import math
import numpy as np
from mpmath import mpc

from . import nodal as N
from .uns import theta_prime


def J(z):
    """J(z) = i / conj(z)."""
    z = np.asarray(z, dtype=complex)
    return 1j / np.conj(z)


def J_cycle(z):
    """(z, J z, J² z, J³ z, J⁴ z)."""
    out = [complex(z)]
    for _ in range(4):
        out.append(complex(J(out[-1])))
    return out


def J_fixed_point_defect(n=20000, seed=0):
    """min over random z of |J(z) − z|: bounded away from 0 (no fixed point); on r = 1 it is exactly |i/z̄ − z| = |z|·|i − 1|·… = √2."""
    rng = np.random.default_rng(seed)
    z = rng.normal(size=n) + 1j * rng.normal(size=n)
    th = rng.uniform(0, 2 * np.pi, n)
    on_circle = np.exp(1j * th)
    return float(np.min(np.abs(J(z) - z))), float(np.min(np.abs(J(on_circle) - on_circle))), float(np.max(np.abs(J(on_circle) - on_circle)))


def w_to_s(w):
    return 1 / (1 - np.asarray(w, dtype=complex))


def xi_ratio_under(op, s_points):
    """|ξ(s′)| / |ξ(s)| with s′ = Cayley⁻¹(op(Cayley s)): 1 everywhere iff op is a symmetry of |ξ|."""
    out = []
    for s in s_points:
        w = (s - 1) / s
        s2 = w_to_s(op(w))
        out.append(float(abs(N.xi(mpc(s2))) / abs(N.xi(mpc(s)))))
    return np.array(out)


def mirror(w):
    """w ↦ 1/w̄ : the north–south mirror of this repo (s ↦ 1 − s̄)."""
    return 1 / np.conj(np.asarray(w, dtype=complex))


def energy_gemini(k):
    """∫₀^{2π}∫₀¹ cos²(k ln r) cos²(mθ + π/4) r dr dθ  (m ≥ 1) = (π/4)(1 + 1/(1+k²)), closed form; and the area ∫∫ r dr dθ = π."""
    return math.pi / 4 * (1 + 1 / (1 + k * k)), math.pi


def energy_numeric(k, n=200001):
    """The radial factor by quadrature in u = ln r (∫_{−∞}^0 cos²(ku) e^{2u} du), times the angular π."""
    u = np.linspace(-14, 0, n)
    f = np.cos(k * u) ** 2 * np.exp(2 * u)
    return float(np.trapezoid(f, u)) * math.pi


def mean_transverse_growth(sigma, t0, t1, n=120):
    """Mean over t ∈ [t0,t1] of the growth rate across the equator ∂σ ln|ξ| at latitude σ, and the mean of θ′(t) = ½ln(t/2π)-ish over the same window.
    Across σ the profile is a GROWTH (cosh/sinh, Cauchy–Riemann) at the rate ≈ θ′(t) + Re ζ′/ζ-terms, not an oscillation cos(k ln r)."""
    from .chladni import rel_flux_grid
    t = np.linspace(t0, t1, n)
    J = rel_flux_grid(np.array([sigma]), t)[:, 0]
    return float(J.mean()), float(np.mean([theta_prime(x) for x in t]))


def flux_formula(s):
    """EXACT:  ∂σ ln|ξ(s)| = Re ξ′/ξ = Re[ 1/s + 1/(s−1) ] − ½ ln π + ½ Re ψ(s/2) + Re ζ′/ζ(s)   (from ξ = ½ s(s−1) π^{−s/2} Γ(s/2) ζ(s)).
    With θ′(t) = ½ Re ψ(¼ + it/2) − ½ ln π this is  θ′(t) + R(s),  R = Re[1/s + 1/(s−1)] + ½[Re ψ(s/2) − Re ψ(¼+it/2)] + Re ζ′/ζ(s)."""
    from mpmath import mp, digamma, zeta, log, pi
    s = mpc(s)
    with mp.workdps(30):
        zp = zeta(s, derivative=1) / zeta(s)
        return float((1 / s + 1 / (s - 1)).real - log(pi) / 2 + digamma(s / 2).real / 2 + zp.real)


def flux_identity_defect(sigma, t, h=1e-6):
    """|numeric ∂σ ln|ξ| − formula| (central difference): the identity is exact, so this is finite-difference error only."""
    from mpmath import mp, log
    a = abs(N.xi(mpc(sigma + h, t))); b = abs(N.xi(mpc(sigma - h, t)))
    num = float((log(a) - log(b)) / (2 * h))
    return abs(num - flux_formula(mpc(sigma, t)))


def remainder_terms(sigma, t):
    """The pieces of R(s) at one point: (Re[1/s+1/(s−1)], ½[Reψ(s/2) − Reψ(¼+it/2)], Re ζ′/ζ): only the last is not O(1/t) or O(σ−½)/t."""
    from mpmath import mp, digamma, zeta
    s = mpc(sigma, t)
    with mp.workdps(30):
        a = float((1 / s + 1 / (s - 1)).real)
        b = float((digamma(s / 2).real - digamma(mpc(0.25, t / 2)).real) / 2)
        c = float((zeta(s, derivative=1) / zeta(s)).real)
    return a, b, c


# ── the user's own J_N (ValaQuenta engine/units.py) against J above ─────────────────────────
def engine_J_agrees(n=200, seed=3):
    """Max difference between ValaQuenta's inversion_transform (r,θ) → (1/r, θ+π/2) and J = i/z̄ in polar form; plus whether the engine's own
    comment is right that only the RADIUS closes after two strokes.  Returns None if ValaQuenta is not on disk."""
    import importlib.util, pathlib, os, sys, types
    root = pathlib.Path(os.environ.get("VALAQUENTA", "/home/rendier/Projects/ThePlace/ValaQuenta"))
    f = root / "engine" / "units.py"
    if not f.exists():
        return None
    src = f.read_text().replace("from . import constants as C", "import math as _m\nclass C:\n    PI = _m.pi")
    mod = types.ModuleType("vq_units"); 
    try:
        exec(compile(src, str(f), "exec"), mod.__dict__)
    except Exception:
        return None
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(n):
        r = float(rng.uniform(0.2, 5)); th = float(rng.uniform(0, 2 * np.pi))
        rr, tt = mod.inversion_transform(r, th)
        z = J(r * np.exp(1j * th))
        worst = max(worst, abs(rr * np.exp(1j * tt) - z))
    (r2, t2), r_closes = mod.inversion_involution(1.7, 0.4)
    return {"max_diff": float(worst), "radius_closes_after_two": bool(r_closes), "angle_after_two": float(t2 - 0.4)}
