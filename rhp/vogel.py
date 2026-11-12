# SPDX-License-Identifier: GPL-3.0-only
"""Vogel / phyllotaxis spirals with the step angle given in DEGREES: decoding the 'Hamiltonian Convergence' widget (Gemini, 2026-04-11, screenshots Gem_Ham_1/2).   :CALCULATED (exact).

Point n sits at angle n·β degrees, radius √n.  The visible arms (parastichies) are the continued-fraction convergents p/q of β/360: q arms while the drift N·|β/360 − p/q| over the N points stays below ~1 turn.
Screenshots: α⁻¹ slider 130 → many arms (130/360 = 13/36 exactly: 36 arms); 143.494 → five straight arms (2/5 within 0.0014 turn).  The golden angle 360/φ² and α⁻¹ (read as degrees) share the first four partial quotients.
"""
from __future__ import annotations
import math

PHI = (1 + 5 ** .5) / 2
GOLDEN_DEG = 360 / PHI ** 2
ALPHA_INV = 137.035999084


def cf_convergents(x, n=10):
    """Continued fraction terms and convergents (p, q, |x − p/q|) of x."""
    a_s, conv = [], []; p0, q0, p1, q1 = 0, 1, 1, 0; y = x
    for _ in range(n):
        a = math.floor(y); a_s.append(a); p0, p1 = p1, a * p1 + p0; q0, q1 = q1, a * q1 + q0; conv.append((p1, q1, abs(x - p1 / q1)))
        if abs(y - a) < 1e-12: break
        y = 1 / (y - a)
    return a_s, conv


def arms(beta_deg, N=700, qmax=60):
    """[(q, drift in turns over N points)] for the convergent denominators q ≤ qmax of β/360."""
    _, c = cf_convergents(beta_deg / 360, 12); return [(q, N * e) for p, q, e in c if q <= qmax]


def drift_turns(beta1_deg, beta2_deg):
    """N at which two Vogel spirals with step angles β1, β2 have drifted apart by one full turn."""
    return 360 / abs(beta1_deg - beta2_deg)
