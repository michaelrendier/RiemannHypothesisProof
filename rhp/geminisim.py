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


# ── the ORIGINAL sim: Ainulindale/AgeFirst/original_2026/ResearchPaper/event_horizon_sim.py (143 lines, Gemini) ──────────────
E_INV = 2 / np.e


def original_phase(phase, t, n=5000, pi_f=PI_INV, e_f=E_INV):
    """(θ, r) of the original file's phases 0-2 at cycle time t ∈ [0, 1] (slider values default to 2/π and 2/e).
    0: Orientation — θ = 10π·u·pi_f, r = linspace(10, 0.1)·(1 − t/2).
    1: Metric swap (the WORM) — r = t·e_f/(r_b + 0.01), θ = θ_b + 5·sin(θ_b·t),  θ_b ∈ [0, 10π], r_b ∈ [10, 0.1].
    2: Fermat frequency (the ROSETTE) — θ ∈ [0, 100π], r = √θ·(1 + 0.2·sin(θ·pi_f·t))."""
    if phase == 0:
        return np.linspace(0, 10 * np.pi, n) * pi_f, np.linspace(10, 0.1, n) * (1 - t * .5)
    if phase == 1:
        tb = np.linspace(0, 10 * np.pi, n); rb = np.linspace(10, 0.1, n)
        return tb + np.sin(tb * t) * 5, (1 / (rb + 0.01)) * t * e_f
    th = np.linspace(0, 100 * np.pi, n)
    return th, np.sqrt(th) * (1 + 0.2 * np.sin(th * pi_f * t))


def worm_sidebands(t, beta=5.0, n=20000, thmax=10 * np.pi, nmax=8):
    """Phase 2's unit phasor e^{iθ}, θ = θ_b + β sin(ωθ_b), ω = t, is PHASE modulation: e^{iθ} = Σ_n J_n(β) e^{i(1+nω)θ_b}.
    Returns (n, |J_n(β)|, measured line amplitude at frequency 1 + nω) for n = −nmax…nmax (measured by direct projection)."""
    from scipy.special import jv
    tb = np.linspace(0, thmax, n, endpoint=False); z = np.exp(1j * (tb + beta * np.sin(t * tb)))
    out = []
    for k in range(-nmax, nmax + 1):
        out.append((k, abs(jv(k, beta)), abs(np.mean(z * np.exp(-1j * (1 + k * t) * tb)))))
    return out


def rosette_sidebands(k, n=200000, thmax=400 * np.pi):
    """Phase 3's z/√θ = e^{iθ}(1 + 0.2 sin kθ) is AMPLITUDE modulation: carrier at frequency 1 and sidebands at 1 ± k, each 0.1 the carrier.
    Returns measured amplitudes (carrier, 1+k, 1−k) by direct projection."""
    th = np.linspace(0, thmax, n, endpoint=False); z = np.exp(1j * th) * (1 + .2 * np.sin(k * th))
    f = lambda w: abs(np.mean(z * np.exp(-1j * w * th)))
    return f(1.0), f(1 + k), f(1 - k)


def convergents_ratio(p, q, depth=8):
    """Continued-fraction convergents (a/b) of ln p / ln q ≥ 1 — the near-resonances of two prime wheels.  ln p/ln q is irrational for
    distinct primes (unique factorisation), so no convergent is ever exact: the pathway never closes into a rosette."""
    x0 = x = np.log(p) / np.log(q); out = []; h0, h1, k0, k1 = 0, 1, 1, 0
    for _ in range(depth):
        a = int(np.floor(x)); h0, h1 = h1, a * h1 + h0; k0, k1 = k1, a * k1 + k0
        out.append((h1, k1, abs(x0 - h1 / k1)))
        if x - a < 1e-12: break
        x = 1 / (x - a)
    return out
