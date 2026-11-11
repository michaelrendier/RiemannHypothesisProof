# SPDX-License-Identifier: GPL-3.0-only
"""The bit flip: what sets bit j of a prime's slot, computed two independent ways, and what disagreement says about composites.   :CALCULATED

Bit j of the `char` slot of p is 1 exactly when the Kronecker symbol (d_j/p) = −1 (p stays a hole in Q(√d_j)), d = −4, 8, −3, 5, −7, −11, 13.
Two routes to the same bit for ODD n (prime or not):
  exponent  e_j = d_j^((n−1)/2) mod n          (Euler's criterion: one modular exponentiation, O(log n) squarings)
  reciprocity  J_j = Jacobi(d_j / n)           (the Euclid-style algorithm with the reciprocity law and the supplements for −1 and 2: no exponent, no factoring)
For n prime the two routes agree on every base (Euler's criterion).  For composite n they agree on a base only if n is an Euler–Jacobi pseudoprime to it (Solovay–Strassen): the slot vector computed
two ways is a built-in complementarity check.  Vectorised int64 numpy (n < 3·10⁹ so products fit)."""
from __future__ import annotations
import numpy as np

DISC = (-4, 8, -3, 5, -7, -11, 13)


def modpow(a, e, n):
    """vectorised a^e mod n, element-wise (arrays of equal length, n < 3·10⁹)."""
    a = np.asarray(a, np.int64) % n; e = np.asarray(e, np.int64).copy(); n = np.asarray(n, np.int64)
    r = np.ones_like(n); b = a.copy()
    while e.any():
        odd = (e & 1) == 1
        r = np.where(odd, r * b % n, r); b = b * b % n; e >>= 1
    return r % n


def jacobi(a, n):
    """vectorised Jacobi symbol (a/n), n odd positive; returns int8 in {−1, 0, 1}.  Reciprocity algorithm: strip factors of 2 (using n mod 8), flip, use (a/n)(n/a) = (−1)^((a−1)(n−1)/4)."""
    a = np.asarray(a, np.int64) % np.asarray(n, np.int64); n = np.asarray(n, np.int64).copy(); a = a.copy()
    res = np.ones(len(n), np.int8)
    for _ in range(200):
        act = a != 0
        if not act.any(): break
        while True:
            ev = act & ((a & 1) == 0)
            if not ev.any(): break
            a = np.where(ev, a >> 1, a); n8 = n & 7
            res = np.where(ev & ((n8 == 3) | (n8 == 5)), -res, res)
        sw = act
        a2 = np.where(sw, n, a); n2 = np.where(sw, a, n)
        flip = sw & ((a2 & 3) == 3) & ((n2 & 3) == 3)
        res = np.where(flip, -res, res)
        a = np.where(sw, a2 % np.where(n2 == 0, 1, n2), a); n = n2
    return np.where(n == 1, res, 0).astype(np.int8)


def two_ways(ns, d):
    """(exponent residue as −1/0/1/other-flag, Jacobi symbol, agree) for base d at every n in ns (odd)."""
    ns = np.asarray(ns, np.int64); e = modpow(d, (ns - 1) // 2, ns); j = jacobi(d, ns).astype(np.int64)
    agree = (e == (j % ns))
    return e, j, agree


def bit_vectors(ns, bases=DISC):
    """per-base arrays: exp_bit (1 iff residue ≡ −1), jac_bit (1 iff symbol = −1), agree."""
    out = {}
    for d in bases:
        e, j, agree = two_ways(ns, d); out[d] = dict(exp_bit=(e == (np.asarray(ns) - 1)).astype(np.int8), jac_bit=(j == -1).astype(np.int8), agree=agree)
    return out
