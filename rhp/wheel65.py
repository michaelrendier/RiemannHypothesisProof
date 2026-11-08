# SPDX-License-Identifier: GPL-3.0-only
"""The 65-wheel coil: the primes ≤ 313 as circles, used as an exact primality test below 317² = 100489.

Wheel p is the circle of circumference p; the integer n sits at angle 2πn/p and the wheel RINGS when n is on its 0 spot (p | n).
Ramanujan sum  c_p(n) = Σ_{a coprime p} e^{2πi a n/p} = p·[p|n] − 1  is that wheel's whole non-trivial spectrum, so the coil's
total response is  S(n) = Σ_p c_p(n) = Σ_{p|n} p − 65.   Silence (S = −65) means no wheel rings.
ESTABLISHED: this is trial division (Eratosthenes); the guarantee holds because 317 is the 66th prime and 317² > 99999 > 313².
OURS-CALCULATED: exhaustive check on every n < 317², and the cost/period/epicycle numbers below.
"""
from __future__ import annotations
import math
import numpy as np
from .wrap import primes_to

PRIMES = primes_to(313)            # 65 wheels
LIMIT = 317 * 317                  # the guarantee covers 2 ≤ n < LIMIT


def response(n):
    """Coil response S(n) = Σ_{p|n} p − 65 (Ramanujan sums of all 65 wheels)."""
    n = np.asarray(n)
    return (PRIMES[:, None] * (n[None, :] % PRIMES[:, None] == 0)).sum(0) - len(PRIMES)


def ringing(n):
    """Wheels (primes) that ring for n."""
    return [int(p) for p in PRIMES if n % p == 0]


def is_prime(n):
    """Exact for 2 ≤ n < 317²: prime iff no wheel rings, or only the wheel p = n itself."""
    if not 2 <= n < LIMIT:
        raise ValueError("outside the coil's guarantee (needs 2 ≤ n < 317²)")
    r = ringing(n)
    return not r or r == [n]


def wheels_read(n):
    """How many wheels are read before the first rings (early exit); a prime reads all 65."""
    for k, p in enumerate(PRIMES, 1):
        if n % p == 0:
            return k
    return len(PRIMES)


def log_period():
    """ln of the combined period of the 65 circles (CRT): ln(313#) = θ(313)."""
    return float(np.sum(np.log(PRIMES)))


def epicycle_rms(sigma=0.5):
    """rms radius of the Euler-wheel epicycle Σ_p p^{−σ−it} (frequencies ln p): √(Σ p^{−2σ})."""
    return float(np.sqrt(np.sum(PRIMES.astype(float) ** (-2 * sigma))))
