# SPDX-License-Identifier: GPL-3.0-only
"""The Fermat side of the horizon, made exact: the Fermat curve xⁿ + yⁿ = zⁿ over F_p.   :CALCULATED (exact integer / 1e-9 float) unless marked ESTABLISHED.

Why here: over F_p the Fermat equation has a finite, exactly known spectrum.  ESTABLISHED (Weil 1949): for p ≡ 1 mod n the affine count is
    N_aff(p) = #{(x,y) ∈ F_p² : xⁿ + yⁿ = 1} = p + Σ_{a,b=1}^{n−1} J(χ^a, χ^b),   J(χ^a, χ^b) = Σ_u χ^a(u) χ^b(1−u),  χ of order n,
with |J| = √p for the (n−1)(n−2) = 2g terms where a + b ≢ 0 (n), and J = −χ^a(−1) (modulus 1) for the n − 1 terms with a + b ≡ 0.  So the 'prime wave' of the Fermat curve is a FINITE sum of 2g tones of fixed
modulus √p (the Riemann hypothesis for this curve is a theorem), the counterpart of ψ(x) = x − Σ x^ρ/ρ with infinitely many tones.  The Weil ceiling is the tone sum's maximum 2g√p.
The Fermat 'horizon' is where the main term p overtakes that ceiling: p + 1 − 2g√p > 3n guarantees a non-trivial point.  Below it a finite set of primes can be unreached: the exceptional primes P_n.
"""
from __future__ import annotations
import math
import numpy as np
from .sedenion_slots import primes_upto


def primitive_root(p):
    p = int(p)
    fs = [q for q in range(2, p) if (p - 1) % q == 0 and all(q % r for r in range(2, int(q ** .5) + 1))]
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in fs): return g
    return 1


def index_table(p):
    p = int(p)
    g = primitive_root(p); ind = np.zeros(p, np.int64); x = 1
    for j in range(p - 1): ind[x] = j; x = x * g % p
    return ind


def affine_count_bruteforce(p, n):
    p = int(p); n = int(n)
    """#{(x, y) ∈ F_p² : xⁿ + yⁿ = 1}."""
    pw = np.array([pow(x, n, p) for x in range(p)]); cnt = np.bincount(pw, minlength=p)
    return int(sum(cnt[a] * cnt[(1 - a) % p] for a in range(p)))


def jacobi_sums(p, n):
    p = int(p); n = int(n)
    """{(a, b): J(χ^a, χ^b)} for 1 ≤ a, b ≤ n−1, p ≡ 1 mod n, χ(g^j) = e^{2πi j/n} (χ(0) = 0)."""
    assert (p - 1) % n == 0
    ind = index_table(p); u = np.arange(2, p); out = {}
    for a in range(1, n):
        for b in range(1, n):
            out[(a, b)] = complex(np.exp(2j * np.pi * (a * ind[u] + b * ind[(1 - u) % p]) / n).sum())
    return out


def count_from_jacobi(p, n):
    """N_aff from the sum formula, rounded; also the (n−1)(n−2) tones' normalised angles and the max | |J| − √p | over them."""
    J = jacobi_sums(p, n); tot = p + sum(J.values()); gen = [v for (a, b), v in J.items() if (a + b) % n != 0]
    return tot, np.angle(np.array(gen) / math.sqrt(p)), max(abs(abs(v) - math.sqrt(p)) for v in gen)


def has_nontrivial(p, n):
    p = int(p); n = int(n)
    """Is there (x, y, z) mod p with xyz ≠ 0 and xⁿ + yⁿ = zⁿ?  (nonzero n-th powers R: ∃ u, v ∈ R with u + v ∈ R)."""
    R = np.unique(np.array([pow(x, n, p) for x in range(1, p)])); mask = np.zeros(p, bool); mask[R] = True
    for u in R:
        if mask[(R + u) % p].any(): return True
    return False


def exceptional_primes(n, pmax):
    """Primes p ≤ pmax, p ∤ n, with no non-trivial solution of xⁿ + yⁿ = zⁿ mod p (xyz ≠ 0): the places the Fermat equation cannot reach."""
    return [int(p) for p in primes_upto(pmax) if n % p and not has_nontrivial(int(p), n)]


def genus(n):
    return (n - 1) * (n - 2) // 2


def weil_threshold(n):
    """p > (g + √(g² + 3n − 1))² guarantees p + 1 − 2g√p > 3n, hence a non-trivial point (the 3n points with xyz = 0 are the trivial ones)."""
    g = genus(n); return (g + math.sqrt(g * g + 3 * n - 1)) ** 2
