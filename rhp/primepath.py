# SPDX-License-Identifier: GPL-3.0-only
"""The pathway of the first K primes: partial Euler product at s = ½ + it, flat, 3-axis and on the sphere; the torus of wheels; survivors; zero-height precision.

Wheel p is the complex circle  z_p = p^{−s} = p^{−σ} e^{−it ln p}  (radius p^{−σ}, rate ln p), one per prime, each in its OWN complex plane C_p.  K wheels = a point of a flat torus
T^K ⊂ C^K, and K ↦ the path is the SUM of the coordinates (a projection onto the diagonal).  log P_K(s) = Σ_{k≤K} −ln(1 − p_k^{−s}) is that path (continuous, no branch jumps since |p^{−s}| < 1);
P_K = exp(log P_K) is the partial Euler product.
ESTABLISHED: the Euler product, Σ_{p≤x} p^{−s} ≈ Li-term Ei((1−s)ln x) (prime number theorem), Mertens Π(1−1/p) ≈ e^{−γ}/ln y.
OURS-CALCULATED: everything measured below.  Nothing here is evidence about RH.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.special import expi

from .wrap import primes_to
from .unsmap import stereo

GAMMA1 = 14.134725141734693
EULER_GAMMA = 0.5772156649015329


def first_primes(K):
    """The first K primes (sieve to a bound from p_K < K(ln K + ln ln K) for K ≥ 6)."""
    lim = 20 if K < 6 else int(K * (math.log(K) + math.log(math.log(K)))) + 10
    return primes_to(lim)[:K].astype(float)


def log_product_path(t, P, sigma=0.5):
    """log P_K(σ+it) for K = 1..len(P): the cumulative Σ −ln(1 − p^{−s})."""
    s = sigma + 1j * t
    return np.cumsum(-np.log1p(-np.exp(-s * np.log(P))))


def prime_sum_path(t, P, sigma=0.5):
    """Σ_{p≤p_K} p^{−s} (prime terms only)."""
    return np.cumsum(np.exp(-(sigma + 1j * t) * np.log(P)))


def smooth_prime_sum(t, x, sigma=0.5):
    """PNT smooth part of Σ_{p≤x} p^{−s}:  Ei((1−s) ln x) − Ei((1−s) ln 2) (grows like x^{1−σ}/((1−s) ln x))."""
    s = sigma + 1j * t
    return expi((1 - s) * np.log(x)) - expi((1 - s) * math.log(2))


def wheel_sphere_radius(P, sigma=0.5):
    """R_K = √Σ p^{−2σ}: the radius of the sphere S^{2K−1} ⊂ C^K that holds the torus of the K wheels in the rms sense (σ = ½: √Σ 1/p ~ √ln ln p_K)."""
    return float(np.sqrt(np.sum(P ** (-2 * sigma))))


def survivor_density(K):
    """Π_{k≤K}(1 − 1/p_k): density of integers with no prime factor among the first K primes (not only primes!)."""
    return float(np.prod(1 - 1 / first_primes(K)))


def mertens_estimate(y):
    return math.exp(-EULER_GAMMA) / math.log(y)


def downsample(n, keep_first=65, m=6000):
    """Indices 0..n−1: all of the first `keep_first` plus about m log-spaced ones."""
    idx = set(range(min(keep_first, n)))
    idx.update(np.unique(np.logspace(math.log10(keep_first), math.log10(n - 1), m).astype(int)).tolist())
    idx.add(n - 1)
    return np.array(sorted(i for i in idx if i < n))


def sphere_path(logP):
    """Stereographic image of P_K = exp(log P_K): |P| = 1 is the equator, P = 0 the north pole, ∞ the south pole."""
    return stereo(np.exp(np.clip(logP.real, -300, 300)) * np.exp(1j * logP.imag))


def equator_crossings(logP):
    """Number of times |P_K| crosses 1 (sign changes of Re log P_K)."""
    sg = np.sign(logP.real)
    return int(np.sum(sg[1:] != sg[:-1]))


def round_zeros(zeros, decimals):
    return np.round(np.asarray(zeros, float), decimals)
