# SPDX-License-Identifier: GPL-3.0-only
"""Conway's PRIMEGAME (FRACTRAN, 1987): 14 fractions that emit the primes as exponents of 2.  ESTABLISHED.

State = the exponent vector over the registers 2,3,5,7 (counters) and 11,13,17,19,23,29 (0/1 flags).
Divisibility by d is the closing of a rotation by d on a circle of n positions: period n/gcd(n,d).
"""
from __future__ import annotations
from math import gcd
import numpy as np

PROG = [(17, 91), (78, 85), (19, 51), (23, 38), (29, 33), (77, 29), (95, 23), (77, 19), (1, 17), (11, 13), (13, 11), (15, 14), (15, 2), (55, 1)]
REGISTERS = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]


def primegame():
    """Generator of the primes (Listing). Pure-integer version of the Fraction loop."""
    n = 2
    while True:
        for a, b in PROG:
            if n % b == 0:
                n = n // b * a
                break
        k, m = 0, n
        while m % 2 == 0:
            m //= 2; k += 1
        if m == 1 and k > 1:
            yield k


def events(nprimes):
    """[(p, step)]: the step at which 2^p first appears, with per-fraction cumulative use counts at that step."""
    n, step, use, out = 2, 0, [0] * 14, []
    while len(out) < nprimes:
        for i, (a, b) in enumerate(PROG):
            if n % b == 0:
                n = n // b * a; use[i] += 1
                break
        step += 1
        k, m = 0, n
        while m % 2 == 0:
            m //= 2; k += 1
        if m == 1 and k > 1:
            out.append((k, step, list(use)))
    return out


def trace(steps):
    """Exponent vectors over REGISTERS for the first `steps` states."""
    n = 2; rows = []
    for _ in range(steps):
        for a, b in PROG:
            if n % b == 0:
                n = n // b * a
                break
        v, m = [], n
        for p in REGISTERS:
            k = 0
            while m % p == 0: m //= p; k += 1
            v.append(k)
        assert m == 1
        rows.append(v)
    return np.array(rows)


def fit_cost(ev):
    p = np.array([e[0] for e in ev], float); s = np.array([e[1] for e in ev], float)
    sel = p >= 11
    slope, icpt = np.polyfit(np.log(p[sel]), np.log(s[sel]), 1)
    return float(slope), float(np.exp(icpt))


def orbit_period(n, d):
    return n // gcd(n, d)


def closing_table(n):
    """For a circle of n positions: d -> period of rotation by d (a closed early orbit means period < n)."""
    return {d: orbit_period(n, d) for d in range(1, n)}
