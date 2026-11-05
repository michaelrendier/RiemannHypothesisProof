# SPDX-License-Identifier: GPL-3.0-only
"""The θ-clock: the pathway in the phase variable θ(t), and the n = 1 'carrier' circle.

ζ ≈ A + χB.  The n = 1 pair gives 1 + e^{−2iθ}: a circle of radius 1 centred at 1 through the origin,
reached when θ ≡ π/2 (mod π).  Its zeros solve θ(E) = (k − 3/2)π  (k = 1, 2, …), the Weyl-law positions.
The wheels n ≥ 2 shift each zero by an amount set by S(γ):  θ(γ_k) = π(k − 3/2 − S(γ_k)).  (ESTABLISHED relation; :CALCULATED here)
"""
from __future__ import annotations
import numpy as np
from . import bigpath


def theta_prime_num(t, h=1e-5):
    t = np.asarray(t, float)
    return (bigpath.theta_fast(t + h) - bigpath.theta_fast(t - h)) / (2 * h)


def carrier_zeros(kmax, t_start=14.0):
    """E_k with θ(E_k) = (k − 3/2)π, k = 1..kmax, on the increasing branch (t > 2π), Newton."""
    k = np.arange(1, kmax + 1)
    target = (k - 1.5) * np.pi
    t = np.maximum(t_start, 2 * np.pi * np.exp(1 + np.lib.scimath.log(np.maximum(k, 2))) * 0 + 14.0 + 6.5 * k / np.log(np.maximum(k, 2) + 3))
    for _ in range(100):
        f = bigpath.theta_fast(t) - target
        t = t - f / theta_prime_num(t)
    assert np.max(np.abs(bigpath.theta_fast(t) - target)) < 1e-9
    return t


def jitter(zeros):
    """S(γ_k) = k − 3/2 − θ(γ_k)/π  (with k = 1, 2, … the zero's index): the zero-counting error at the zero (mean convention)."""
    k = np.arange(1, len(zeros) + 1)
    return k - 1.5 - bigpath.theta_fast(zeros) / np.pi


def shift_in_spacings(zeros):
    """(γ_k − E_k) / (mean local spacing): how far the other wheels move each zero from the carrier's prediction."""
    E = carrier_zeros(len(zeros))
    sp = np.pi / np.maximum(theta_prime_num(E), 1e-9)       # θ advances π per zero
    return (zeros - E) / sp
