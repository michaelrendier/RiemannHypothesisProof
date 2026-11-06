# SPDX-License-Identifier: GPL-3.0-only
"""The equatorial construction: the critical line as the equator of the Riemann sphere.

Orientation : w = (s − 1)/s  (Cayley), s = σ + it.   Aperture: the whole sphere, resolved on σ ∈ (0,1).
ESTABLISHED (coordinate change, checked below by CALCULATED residuals):
  |w| = 1  ⇔  |s − 1| = |s|  ⇔  σ = ½          the critical line is the unit circle = the equator
  |w| < 1  ⇔  σ > ½  (north),  |w| > 1 ⇔ σ < ½  (south)
  s ↦ 1 − s   becomes   w ↦ 1/w ;  with conjugation  w ↦ 1/w̄  = reflection through the equator
  poles: w = 0 is s = 1, w = ∞ is s = 0 — the two poles of the completed Λ(s) = π^{−s/2} Γ(s/2) ζ(s)
  the line s = ½ + it has w = exp(i(π − 2 arctan 2t)): the whole infinite line wraps ONCE around the circle.
OURS (readings, not theorems): the Hardy phase θ(t) winds N(T) = θ/π + 1 + S times — the circle is its flat projection;
  the π-factor of Λ is the −(t/2)·ln π part of θ(t), a constant-rate rotation.
Nothing here is evidence about RH: RH is the statement that every zero lies on this circle.
"""
from __future__ import annotations

import numpy as np
from mpmath import mp, mpf, mpc, zeta, gamma, loggamma, log, pi, siegeltheta, sin, arg

from . import uns

mp.dps = 25


# ── the Cayley map and its two charts ────────────────────────────────────────────────────────
def cayley(s):
    """w = (s − 1)/s."""
    return (s - 1) / s


def cayley_inv(w):
    """s = 1/(1 − w)."""
    return 1 / (1 - w)


def to_xyz(w):
    """Stereographic image of w on the unit sphere, north pole at w = 0 (s = 1), south at ∞ (s = 0)."""
    w = np.asarray(w, dtype=complex)
    d = 1 + np.abs(w) ** 2
    return np.column_stack([2 * w.real / d, 2 * w.imag / d, (1 - np.abs(w) ** 2) / d])


def latitude_cayley(s):
    """Conformal latitude of s on the Cayley sphere (radians, 0 on the equator)."""
    w = np.abs(cayley(np.asarray(s, dtype=complex)))
    return np.arcsin((1 - w ** 2) / (1 + w ** 2))


def latitude_linear(sigma):
    """The chart that literally re-introduces π:  latitude = π(σ − ½).  Poles are the edges σ = 1, σ = 0."""
    return np.pi * (np.asarray(sigma, dtype=float) - 0.5)


def longitude_cayley(t):
    """Longitude of ½ + it on the Cayley equator, in (−π, π]:  π − 2·arctan(2t) (wrapped)."""
    a = np.pi - 2 * np.arctan(2 * np.asarray(t, dtype=float))
    return (a + np.pi) % (2 * np.pi) - np.pi


def longitude_hardy(t):
    """The dynamical longitude: the Hardy phase θ(t) mod 2π.  End-on, the spring is this circle.  Zeros scatter about its
    two antipodal carrier points (θ ≡ π/2 mod π) by the jitter S(γ): measured, frac(θ(γ)/π) has mean 0.50 and sd 0.25
    (a uniform scatter would be 0.29) — a weak bunching, NOT two spots.  The carrier statement is about INDEX matching."""
    return np.array([float(siegeltheta(mpf(x))) for x in np.atleast_1d(t)]) % (2 * np.pi)


def xyz_from(lat, lon):
    """Unit-sphere point for a latitude and longitude (radians)."""
    lat, lon = np.asarray(lat, float), np.asarray(lon, float)
    return np.column_stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])


def winding_hardy(t):
    """Turns of the Hardy phase, θ(t)/π (a zero every unit once S is added) — the winding the circle forgets."""
    return np.array([float(siegeltheta(mpf(x))) for x in np.atleast_1d(t)]) / np.pi


# ── completed zeta and the checks on the construction ───────────────────────────────────────
def Lam(s):
    """Λ(s) = π^{−s/2} Γ(s/2) ζ(s)."""
    with mp.workdps(30):                      # other rhp modules lower mp.dps; do not inherit their precision
        s = mpc(s)
        return pi ** (-s / 2) * gamma(s / 2) * zeta(s)


def mirror_residual(s):
    """|Λ(s) − Λ(1 − s)| / |Λ(s)|  — the functional equation, i.e. w ↦ 1/w."""
    with mp.workdps(30):
        a, b = Lam(s), Lam(1 - mpc(s))
        return float(abs(a - b) / abs(a))


def pole_residues(eps=1e-7):
    """(s·Λ(s) at s → 0,  (s − 1)·Λ(s) at s → 1): expected (−1, +1) — the poles of the sphere."""
    e = mpf(eps)
    return float((e * Lam(e)).real), float((e * Lam(1 + e)).real)


def equator_residual(gammas):
    """max ||w| − 1| over the zeros ½ + iγ (a coordinate identity: it measures only float error)."""
    w = cayley(0.5 + 1j * np.asarray(gammas, dtype=float))
    return float(np.max(np.abs(np.abs(w) - 1)))


def zeros_after_mirror(rho):
    """The four-fold orbit of a zero ρ under the sphere's symmetries: ρ, 1 − ρ (north↔south), ρ̄, 1 − ρ̄."""
    r = complex(rho)
    return [r, 1 - r, r.conjugate(), 1 - r.conjugate()]


# ── "π factored out": the phase offset ───────────────────────────────────────────────────────
def theta_split(t):
    """θ(t) = [Im lnΓ(¼ + it/2)] + [−(t/2)·ln π].  Returns (gamma_part, pi_part, total_from_mpmath).

    The π-part is a constant-rate rotation of the equator; the Γ-part carries all of the t·ln t growth."""
    t = np.atleast_1d(np.asarray(t, dtype=float))
    g = np.array([float(loggamma(mpc(0.25, x / 2)).imag) for x in t])
    p = -t / 2 * float(log(pi))
    tot = np.array([float(siegeltheta(mpf(x))) for x in t])
    return g, p, tot


def zeta_pi_factored(s):
    """ζ(s) = π^{s/2} · [Λ(s)/Γ(s/2)]: returns (π^{s/2}, the remainder).  Product reproduces ζ."""
    s = mpc(s)
    return pi ** (s / 2), Lam(s) / gamma(s / 2)


# ── the two counter-rotating waves of a term, and the polarisation they make ─────────────────
def abs_chi(sigma, t):
    """|χ(σ + it)| exactly, χ(s) = 2^s π^{s−1} sin(πs/2) Γ(1 − s)  (ζ(s) = χ(s) ζ(1 − s)).  ≈ (t/2π)^{½−σ}."""
    with mp.workdps(30):
        s = mpc(sigma, t)
        return float(abs(2 ** s * pi ** (s - 1) * sin(pi * s / 2) * gamma(1 - s)))


def amplitude_ratio(n, sigma, t):
    """r = |χ|·n^{2σ−1}: the amplitude of the mirror wave over the forward wave for the term n (σ + it).

    r = 1 on the equator for every n (equal counter-rotating waves = a standing wave), and off the equator only where
    n ≈ √(t/2π) — the Riemann–Siegel cut-off."""
    return abs_chi(sigma, t) * float(n) ** (2 * sigma - 1)


def polarisation(r):
    """Axial ratio b/a of the ellipse traced by two counter-rotating phasors of amplitudes 1 and r:  |1−r|/(1+r).
    0 = linear (a pure standing wave), 1 = circular (a pure travelling wave)."""
    return abs(1 - r) / (1 + r)


def ellipse_axes(a1, a2):
    """Semi-axes of the ellipse traced by a1·e^{iφ} + a2·e^{−iφ}:  (a1 + a2, |a1 − a2|)."""
    return a1 + a2, abs(a1 - a2)
