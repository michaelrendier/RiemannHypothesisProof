# SPDX-License-Identifier: GPL-3.0-only
"""The Heartbeat: what the explicit formula's standing wave contains under its curve, and how ADD / SCALE / SIGN act on it.

AUTHORITATIVE SOURCE (author's own code): ValaQuenta/modules/bao_mass_gap/maths.py ('one standing wave per non-trivial zero') and
archimedes_screw/maths.py:  ψ(x) = x − Σ_ρ x^ρ/ρ − ln 2π − ½ ln(1 − x⁻²).   With u = ln x and E(u) = (x − ψ − ln 2π − …)/√x:
        E(u) = Σ_ρ e^{iγu}/ρ  (on the line)  =  Σ_{γ>0} (2/|ρ|) cos(γu − arg ρ).
Each zero pair is ONE standing wave = TWO COUNTER-ROTATING CIRCLES (ρ and ρ̄) of equal radius 1/|ρ| — the circle area is π/|ρ|².
ESTABLISHED (Hadamard product; Davenport ch. 12):  Σ_ρ 1/ρ = 1 + γ_E/2 − ½ ln 4π = 0.0230957…  ( = −B ),  hence  Σ_{γ>0} 1/|ρ|² = Σ_ρ 1/(ρ(1−ρ))/2 = same number
  when every zero is on the line (|ρ|² = ρ(1−ρ)).
OURS-CALCULATED:
  • area budget  π = [Σ_{γ>0} π/|ρ|²] + π·(½ ln 4π − γ_E/2) :  the zero circles hold 2.31 % of the unit circle's area; the rest is the Γ/π 'ground state' term;
  • mean power  ⟨E²⟩ = 2 Σ_{γ>0} 1/|ρ|² = 2 + γ_E − ln 4π = 0.046191…  (Parseval; measured on the exact ψ below);
  • ADD, SCALE, SIGN on one tone: a step Δu in ln x multiplies the tone by e^{(½+iγ)Δu}: SCALE e^{Δu/2} with ROTATION γΔu (a spiral similarity: ADD in u acts on the
    tone plane as SCALE∘rotation); conjugation ρ ↔ ρ̄ is SIGN (reflection); off the line the two circles have radii e^{±δu}: not a standing wave but an ellipse.
READING (not a theorem): the standing wave is the balanced state of two counter-rotating circles; the imbalance is the Telperion/Laurelin difference.
"""
from __future__ import annotations

import math
import numpy as np

EULER_GAMMA = 0.5772156649015329
SUM_INV_RHO = 1 + EULER_GAMMA / 2 - 0.5 * math.log(4 * math.pi)        # Σ_ρ 1/ρ = Σ_{γ>0} 1/(¼+γ²)
GROUND_SHARE = 0.5 * math.log(4 * math.pi) - EULER_GAMMA / 2            # 1 − SUM_INV_RHO
MEAN_POWER = 2 * SUM_INV_RHO                                            # 2 + γ_E − ln 4π


def tone_radius(gamma):
    return 1 / np.sqrt(0.25 + np.asarray(gamma, float) ** 2)


def tone_area(gamma):
    """Area of the circle of one zero's tone: π/|ρ|² = π/(¼+γ²)."""
    return math.pi * tone_radius(gamma) ** 2


def partial_share(zeros, K):
    """Σ_{k≤K} 1/|ρ_k|²: the share of the unit circle's area held by the first K zero circles (relative to π)."""
    return float(np.sum(tone_radius(zeros[:K]) ** 2))


def tail_estimate(T):
    """Σ_{γ>T} 1/(¼+γ²) ≈ (1/2π)∫_T^∞ ln(t/2π)/t² dt = (ln(T/2π) + 1)/(2πT)  (Riemann–von Mangoldt density)."""
    return (math.log(T / (2 * math.pi)) + 1) / (2 * math.pi * T)


def budget(zeros):
    """(partial sum over the given zeros, with tail estimate, exact value 1+γ_E/2−½ln4π)."""
    s = float(np.sum(tone_radius(zeros) ** 2))
    return s, s + tail_estimate(float(zeros[-1])), SUM_INV_RHO


def mean_power_measured(psi, u0=6.0, u1=16.0, n=200000):
    """⟨E²⟩ over u ∈ [u0,u1] on the exact ψ (long-run mean → 2 + γ_E − ln 4π when every zero is on the line)."""
    from .diverge import residual
    u = np.linspace(u0, u1, n)
    E = residual(u, psi)
    return float(np.mean(E ** 2)), float(np.mean(E))


def tone(u, gamma, scale_exp=0.5):
    """One zero's complex tone e^{(scale_exp + iγ)u}/ρ (its real part ×2 is the real standing wave when scale_exp = ½, divided by e^{u/2})."""
    rho = 0.5 + 1j * gamma
    return np.exp((scale_exp + 1j * gamma) * np.asarray(u, float)) / rho


def spiral_similarity(gamma, du):
    """The map tone(u) ↦ tone(u+du) is multiplication by e^{(½+iγ)du}: returns (scale factor, rotation angle)."""
    z = complex(np.exp((0.5 + 1j * gamma) * du))
    return abs(z), float(np.angle(z))


def counter_rotating_pair(u, gamma, delta=0.0, phase=0.0):
    """The two circles of one standing wave: f = e^{δu} e^{i(γu+φ)}/ρ (forward, +γ) and its conjugate m = conj(f) (the mirror ρ̄, −γ).
    They have EQUAL radius at every u, so f + m = 2 Re f is a REAL wave; on the line (δ = 0) the radius is constant 1/|ρ| (a standing wave in the normalised
    variable); off the line the common radius grows as e^{δu}: a dilation, not a different shape.  Returns (f, m, f + m)."""
    u = np.asarray(u, float)
    rho = 0.5 + delta + 1j * gamma
    f = np.exp(delta * u) * np.exp(1j * (gamma * u + phase)) / rho
    m = np.conj(f)
    return f, m, f + m
