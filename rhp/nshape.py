# SPDX-License-Identifier: GPL-3.0-only
"""The Shape: the Extinguishable Bulk (Laurelin) wrapped around the holes (Telperion), measured.   :CALCULATED / ESTABLISHED as marked.

Three exact descriptions of 'a bulk of everything a factor can reach, with the primes as the holes factoring cannot reach':
 1  the Ω-onion      layer k = {n : Ω(n) = k} (Ω = prime factors with multiplicity).  Primes = Ω 1 (the core), semiprimes wrap them (Ω 2) …  ESTABLISHED: Landau π_k(x) ~ x (ln ln x)^{k−1}/((k−1)! ln x)
 2  two factoring conics  hyperbola a² − b² = N (Fermat): a lattice point with a − b ≥ 3 exists  ⇔  N odd composite;  circle x² + y² = N (Euler): two primitive representations  ⇒  gcd(N, x₁y₂ − x₂y₁) is a factor
 3  the tower          the lattices of the division-algebra integers: ℤ (units 2), ℤ[i] = Z² (units 4), the Hurwitz quaternions = D4 (units 24), the integral octonions = E8 (units 240).  A rational prime p is a hole of a level when it is irreducible there.
    ESTABLISHED (theta series): norm-n counts  r₂ = 4 Σ_{d|n} χ₋₄(d),  Hurwitz 24 Σ_{d|n, d odd} d,  E8 240 σ₃(n).
"""
from __future__ import annotations
import math
import numpy as np


# ------------------------------------------------------------------ 1  the Ω-onion
def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for p in range(2, int(n ** .5) + 1):
        if s[p]: s[p * p::p] = False
    return np.nonzero(s)[0]


def big_omega(X):
    """Ω(n) for n = 0 … X (Ω(0) = Ω(1) = 0), uint8, by dividing out every prime ≤ √X and counting the leftover prime."""
    rem = np.arange(X + 1, dtype=np.int64); om = np.zeros(X + 1, np.uint8)
    for p in primes_upto(int(X ** .5)):
        q = int(p)
        while q <= X:
            rem[q::q] //= p; om[q::q] += 1; q *= int(p)
    om[rem > 1] += 1; om[:2] = 0
    return om


def layer_counts(om, xs, kmax=9):
    """counts[k] = #{2 ≤ n ≤ x : Ω(n) = k} for k = 1 … kmax (last = kmax and above), for each x in xs; also the count of n=1."""
    out = {}
    for x in xs:
        c = np.bincount(np.minimum(om[2:x + 1], kmax), minlength=kmax + 1)
        out[int(x)] = [int(v) for v in c[1:]]
    return out


def landau(x, k):
    """Landau's leading term for π_k(x) (n ≤ x with Ω(n) = k)."""
    L = math.log(math.log(x)); return x * L ** (k - 1) / (math.factorial(k - 1) * math.log(x))


# ------------------------------------------------------------------ 2  the two factoring conics
def hyperbola_points(X):
    """Odd N ≤ X with a lattice point (a, b), a − b = p ≥ 3, a + b = q ≥ p, on a² − b² = N, enumerated from the lattice (not from a factor list): N = p·q, p, q odd."""
    reach = np.zeros(X + 1, bool); p = 3
    while p * p <= X:
        qs = np.arange(p, X // p + 1, 2); reach[p * qs] = True; p += 2
    return reach


def hyperbola_search(Nmax):
    """For each odd N ≤ Nmax a direct search over a (a² − N a perfect square, a − b ≥ 3); no factoring.  Returns the boolean array 'has a non-trivial point'."""
    out = np.zeros(Nmax + 1, bool)
    for N in range(9, Nmax + 1, 2):
        a = np.arange(math.isqrt(N) + (0 if math.isqrt(N) ** 2 == N else 1), (N - 1) // 2 + 1)
        b2 = a * a - N; b = np.sqrt(b2).astype(np.int64)
        ok = (b * b == b2) & (a - b >= 3); out[N] = bool(ok.any())
    return out


def circle_primitive_reps(X):
    """number of primitive representations N = x² + y², x > y ≥ 1, gcd = 1, for every N ≤ X  (array of ints)."""
    cnt = np.zeros(X + 1, np.int32); r = math.isqrt(X)
    xs = np.arange(1, r + 1)
    for y in range(1, r + 1):
        x = xs[(xs > y) & (xs * xs + y * y <= X)]
        x = x[np.gcd(x, y) == 1]; np.add.at(cnt, x * x + y * y, 1)
    return cnt


def euler_factor(N):
    """Euler's method from two primitive representations: returns a factor of N or None."""
    reps = []; r = math.isqrt(N)
    for y in range(1, r + 1):
        x2 = N - y * y
        x = math.isqrt(x2)
        if x * x == x2 and x > y and math.gcd(x, y) == 1: reps.append((x, y))
    if len(reps) < 2: return None
    (a, b), (c, d) = reps[0], reps[1]
    g = math.gcd(N, a * d - b * c)
    return g if 1 < g < N else math.gcd(N, a * c + b * d)


# ------------------------------------------------------------------ 3  the tower
def sum_of_squares_exists(p, k):
    """does p have a representation as a sum of k squares of integers with a non-unit element: any with norm p > 1 (p = N(q), q a non-unit)."""
    if k == 1: return False                      # ℤ: N(x) = x², p is not a square: p is irreducible
    if k == 2:
        a = 0
        while a * a <= p:
            b2 = p - a * a; b = math.isqrt(b2)
            if b * b == b2 and a >= 0: return True
            a += 1
        return False
    return True                                  # k ≥ 4: Lagrange (checked by search in the experiment for small p)


def count_hurwitz_norm(n):
    """Number of Hurwitz quaternions of norm n: integer (a,b,c,d) with a²+b²+c²+d² = n plus half-odd quadruples with Σ(x²) = n."""
    r = math.isqrt(n) + 1; cnt = 0
    rng = range(-r, r + 1)
    for a in rng:
        for b in rng:
            s = a * a + b * b
            if s > n: continue
            for c in rng:
                s2 = s + c * c
                if s2 > n: continue
                d2 = n - s2; d = math.isqrt(d2)
                if d * d == d2: cnt += 1 if d == 0 else 2
    # half-integers: coordinates (2m+1)/2, i.e. (u/2) with u odd, Σ u² = 4n
    ru = int(2 * math.sqrt(n)) + 2; odds = [u for u in range(-ru, ru + 1) if u % 2]
    half = 0
    for u1 in odds:
        for u2 in odds:
            t = u1 * u1 + u2 * u2
            if t > 4 * n: continue
            for u3 in odds:
                t3 = t + u3 * u3
                if t3 > 4 * n: continue
                u4sq = 4 * n - t3; u4 = math.isqrt(u4sq)
                if u4 * u4 == u4sq and u4 % 2: half += 1 if u4 == 0 else 2
    return cnt + half


def e8_shell(m):
    """Vectors of E8 (all integer or all half-integer coordinates, coordinate sum even) with squared length m, brute force over the 8-dimensional lattice (m ≤ 8)."""
    def ints(depth, remaining, total):
        if depth == 8:
            return 1 if remaining == 0 and total % 2 == 0 else 0
        r = math.isqrt(remaining); s = 0
        for x in range(-r, r + 1): s += ints(depth + 1, remaining - x * x, total + x)
        return s
    def halves(depth, remaining4, total):                # coordinates u/2, u odd; remaining4 = 4m − Σu²
        if depth == 8:
            return 1 if remaining4 == 0 and (total // 2) % 2 == 0 and total % 2 == 0 else 0
        r = math.isqrt(remaining4); s = 0
        for u in range(-r - 1, r + 2):
            if u % 2 and u * u <= remaining4: s += halves(depth + 1, remaining4 - u * u, total + u)
        return s
    return ints(0, m, 0) + halves(0, 4 * m, 0)
