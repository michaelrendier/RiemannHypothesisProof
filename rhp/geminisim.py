# SPDX-License-Identifier: GPL-3.0-only
"""The Gemini event-horizon sims, as plain functions (no GUI): what their curves ARE.

Source (read, not edited):  Ainulindale/AgeFirst/original_2026/ResearchPaper/FirstDraft/FirstAge/Gemini/
  Gemini_event_horizon_sim.py  (oldest, 37 lines): z_{n+1} = 1 + 1/z_n  drawn in polar at θ = linspace(0, π), between r = 1 and r = φ.
  Gemini_Transition_EHS.py     (65 lines): an inertial path r = t from 1.5 to 0.5, θ = t·(2/π)·π; below r = 1 it is inverted r → (1+GAP)/r, θ → θ + step.
  Gemini_Revised_2_EHS.py      (135 lines): four phases of 600 points — Inertial, Metric Swap, Heartbeat, Star — as get_coords(phase, t).
  event_horizon_sim.py         (470 lines, Claude Sonnet, April 2026): 5 phases on top of the same constants; golden-angle prime spiral.
Labels: transcription of Gemini's formulas (OURS-CALCULATED here: the measurements on them); no claim about ζ.
"""
from __future__ import annotations
import numpy as np

PI_INV = 2 / np.pi
PHI = (1 + 5 ** 0.5) / 2
GAP = 0.00070
GOLDEN_ANGLE = np.pi * (3 - np.sqrt(5))


def phi_flow(n=12):
    z = [1.0]
    for _ in range(n):
        z.append(1 + 1 / z[-1])
    return np.array(z)


def revised2(phase, t, N=600, primes=None):
    """(θ, r) of Gemini_Revised_2_EHS.get_coords for n = 1..N."""
    n = np.arange(1, N + 1)
    if phase == "Inertial":
        return n * (np.pi / 2) * PI_INV, n * 0.02 * t
    if phase == "Metric Swap":
        base = n * 0.02 * t
        return n * PI_INV + np.pi / 2, 1.0 / (base + 0.1)
    if phase == "Heartbeat":
        pulse = (np.sin(t * 0.5) + 1) / 2
        r_target = 1.0 + (PHI - 1.0) * pulse
        return n * PI_INV * t, (np.sqrt(n) * 0.5) * (r_target + GAP)
    p = primes if primes is not None else _primes(N)
    return p * GOLDEN_ANGLE * t, np.sqrt(p) * 0.4


def _primes(N):
    s = np.ones(20 * N, bool); s[:2] = False
    for i in range(2, int(len(s) ** .5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.flatnonzero(s)[:N]


def transition(h=1.0, n=200):
    t = np.linspace(1.5, 0.5, n); th = t * PI_INV * np.pi
    out = t >= 1; inn = ~out
    return (th[out], t[out]), (th[inn] + (np.pi / 2) * h, (1 / t[inn]) * (1 + GAP))


def heartbeat_xy(t, N=600):
    th, r = revised2("Heartbeat", t, N)
    return r * np.cos(th), r * np.sin(th)


def turns_per_index(t):
    """Heartbeat angle step per index in turns: (2/π)·t / 2π = t/π².  The swirl has ≈ 1/ω arms (t = 1: 9.87 ≈ 10 arms, slowly drifting)."""
    return t / np.pi ** 2


def fib_ratio_error(n=12):
    """|z_n − φ| of the oldest sim's flow z ↦ 1 + 1/z (z_n = F_{n+2}/F_{n+1}); error shrinks by 1/φ² per step."""
    return abs(phi_flow(n)[-1] - PHI)
