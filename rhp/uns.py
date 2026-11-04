# SPDX-License-Identifier: GPL-3.0-only
"""The flat pathway ζ(½+it) and its Universal-Native-Space (UNS) coordinates.

Orientation : log-polar in the ζ-plane (radius log10|ζ|, angle arg ζ) with t explicit;
              sphere chart polar = π·σ (equator σ = ½), azimuth = t mod 2π (winding kept apart).
Aperture    : the equator σ = ½, continuous in t.

Everything here is read off Hardy's identity  ζ(½+it) = e^{-iθ(t)} Z(t),  Z real:
  angle  of the flat point = -θ(t)         (+π where Z < 0)
  radius of the flat point = |Z(t)|        (signed radius ⇒ no corners at the origin)
  areal velocity           = Im(ζ̄ ζ′)     = -θ′(t) Z(t)²
Labels: ESTABLISHED (Hardy, Riemann–Siegel θ); the identities below are CALCULATED here.
"""
from __future__ import annotations

import numpy as np
from mpmath import mp, mpf, mpc, zeta, siegeltheta, digamma, log, pi, zetazero, exp

mp.dps = 20
LN10 = float(log(10))
HALF = mpf(1) / 2


def theta(t):
    """Riemann–Siegel θ(t)."""
    return float(siegeltheta(mpf(t)))


def theta_prime(t):
    """θ′(t) = ½·Re ψ(¼ + it/2) − ½·ln π   (exact)."""
    return float(HALF * digamma(mpc(0.25, mpf(t) / 2)).real - HALF * log(pi))


def zeta_line(t):
    """ζ(½ + it) as a Python complex."""
    return complex(zeta(mpc(HALF, mpf(t))))


def pathway(T=50.0, n=4001, t0=0.0):
    """Sample the pathway. Returns a dict of numpy arrays, one row per t.

    t, zeta (flat point), theta, thetap, Z (Hardy), u = log10|ζ|, phi_signed = −θ,
    L = −θ′Z² (areal velocity), sphere_polar = π/2, sphere_az = t mod 2π.
    """
    t = np.linspace(t0, T, n)
    z = np.array([zeta_line(x) for x in t])
    th = np.array([theta(x) for x in t])
    thp = np.array([theta_prime(x) for x in t])
    Z = np.real(np.exp(1j * th) * z)
    with np.errstate(divide="ignore"):
        u = np.log10(np.abs(z))
    return dict(t=t, zeta=z, theta=th, thetap=thp, Z=Z, u=u, phi_signed=-th,
                L=-thp * Z * Z, sphere_polar=np.full_like(t, np.pi / 2),
                sphere_az=np.mod(t, 2 * np.pi))


def zeros_below(T):
    """Imaginary parts γₙ < T of the nontrivial zeros (mpmath zetazero)."""
    g, n = [], 1
    while True:
        gam = float(zetazero(n).imag)
        if gam >= T:
            return np.array(g)
        g.append(gam)
        n += 1


def gram_points(T, n0=0):
    """Gram points gₙ (θ(gₙ) = nπ), n ≥ n0, below T. Newton on θ with exact θ′."""
    out, n = [], n0
    t = 17.8456 if n0 == 0 else None
    while True:
        # asymptotic start: θ(t) ≈ (t/2)ln(t/2π e) − π/8
        guess = t if t is not None else 2 * np.pi * (n + 0.125) / max(1e-9, np.log(max(n, 2)))
        x = float(guess)
        for _ in range(60):
            x_new = x - (theta(x) - n * np.pi) / theta_prime(x)
            if abs(x_new - x) < 1e-12:
                x = x_new
                break
            x = x_new
        if x >= T:
            return np.array(out)
        out.append(x)
        n += 1
        t = x + 2.0


def cylinder(p, floor=-3.0):
    """UNS cylinder coordinates (ρ, φ, z) of the spring: ρ = log10|ζ| − floor (clipped ≥ 0;
    the floor is a rendering choice — zeros sit on the axis), φ = arg ζ, z = t."""
    rho = np.clip(p["u"] - floor, 0.0, None)
    phi = np.angle(p["zeta"])
    return rho, phi, p["t"]
