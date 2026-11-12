# SPDX-License-Identifier: GPL-3.0-only
"""Nodal structure of the two-sided standing wave: the prime wave ψ(x) − x and the Möbius wave M(x), and the strip field ξ(s) against 1/ξ(s).   :CALCULATED unless marked.

1-D (exact data): W_ψ(x) = (ψ(x) − x)/√x and W_M(x) = M(x)/√x are built on the same zero-tones (the M tone is the ψ tone rescaled by −1/ζ'(ρ)); correlation ≈ −0.82.  Nodes = sign changes.
  locking(A, B): for each node b of B the position inside the bracketing gap of A, d = min(b − a_left, a_right − b)/(a_right − a_left) ∈ [0, ½]; uniformly placed nodes give mean 0.25.
2-D: ESTABLISHED algebra: Re(1/ξ) = Re ξ/|ξ|² and Im(1/ξ) = −Im ξ/|ξ|², so the nodal lines of ξ and of 1/ξ are IDENTICAL; the zeros of ξ are the poles of 1/ξ and sit at the crossings.
Near a zero on the line, ξ(½ + δσ + i(γ + δt)) ≈ ξ_t(δt − iδσ) (ξ is real on the line, ξ_t = dξ/dt real): Re ξ = 0 is the horizontal rung δt = 0, Im ξ = 0 is the critical line δσ = 0: they cross at a right angle.
"""
from __future__ import annotations
import numpy as np


def sign_changes(a):
    """indices n ≥ 1 where the sign of a (zeros skipped) differs from the last nonzero sign."""
    s = np.sign(a); nz = np.nonzero(s)[0]; t = s[nz]; ch = nz[1:][t[1:] != t[:-1]]; return ch


def locking(A, B):
    """mean and distribution of the normalised position of each B node inside the bracketing gap of the A nodes; (mean, array)."""
    A = np.sort(A); i = np.searchsorted(A, B); ok = (i > 0) & (i < len(A)); i, b = i[ok], B[ok]; l, r = A[i - 1], A[i]
    d = np.minimum(b - l, r - b) / (r - l); return float(d.mean()), d


def null_locking(A, nB, lo, hi, rng, reps=20):
    """the same statistic for nB uniformly random nodes in [lo, hi): mean over reps (expected 0.25)."""
    return float(np.mean([locking(A, rng.uniform(lo, hi, nB))[0] for _ in range(reps)]))
