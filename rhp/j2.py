# SPDX-License-Identifier: GPL-3.0-only
"""The J₂ definition, tested: swap involution (wiki 51) against the two-stroke / four-cycle (spin) reading, as operators on T_N.   :CALCULATED / exact as marked.

Definition under test:  J₂ swaps Red and Blue; J₂² = identity (wiki 51), 'two strokes' (wiki 60), 'two rotations to return' (Cody, 2026-11-11).
Colour: blue = indices with bit 4 set (ptol.c shells: red {0–3, 8–11}, blue {4–7, 12–15}); red·red → red, red·blue → blue, blue·blue → red (by XOR of indices).
Candidates: L_k (left multiplication by the blue unit e_k, the grading's own shift), P_k (pure index swap i → i⊕4, all signs +), and the sign patterns between.
Three further objects:
  the Klein group of the functional equation   s → 1−s and s → s̄ ; their product s → 1−s̄ fixes the line; the orbit {ρ, 1−ρ, ρ̄, 1−ρ̄} is 4 off the line, 2 on it
  the arithmetic J₂                            the Galois involution of a quadratic field: a prime is FIXED (inert = irreducible = a hole of that level) or SWAPPED with its conjugate (split)
  the Frobenius reading of the slot rule       slot(p) = Frobenius of p in (ℤ/2)^m = Gal(ℚ(√d₁,…,√d_m)/ℚ); the trivial class = splits completely
"""
from __future__ import annotations
import math
import numpy as np
from .sedenion_slots import Alg, sign_table, kronecker, primes_upto, rule_slots, _DISC


# ------------------------------------------------------------ operators on T_N
def blue(N): return np.array([i for i in range(N) if i & 4]) if N >= 8 else np.array([], int)
def red(N): return np.array([i for i in range(N) if not i & 4])


def J_left(N, k=4): return Alg(N).Lmat(k)


def J_perm(N, k=4):
    P = np.zeros((N, N)); P[np.arange(N) ^ k, np.arange(N)] = 1.0; return P


def order(M, cap=16):
    I = np.eye(len(M)); A = M.copy()
    for n in range(1, cap + 1):
        if np.allclose(A, I): return n
        A = A @ M
    return None


def projective_order(M, cap=16):
    """order up to a global sign (the action on rays)."""
    I = np.eye(len(M)); A = M.copy()
    for n in range(1, cap + 1):
        if np.allclose(A, I) or np.allclose(A, -I): return n
        A = A @ M
    return None


def swaps_colour(M, N):
    r, b = red(N), blue(N)
    return bool(np.allclose(M[np.ix_(r, r)], 0) and np.allclose(M[np.ix_(b, b)], 0))


def mul_basis(N):
    return sign_table(N)


def product(x, y, N):
    S = sign_table(N); out = np.zeros(N); idx = np.arange(N)
    for i in range(N):
        np.add.at(out, idx ^ i, x[i] * S[i] * y)
    return out


def automorphism_defects(M, N, m=300, seed=0):
    """max |M(xy) − M(x)M(y)| and max |M(xy) − M(y)M(x)| over random x, y."""
    rng = np.random.default_rng(seed); a = b = 0.0
    for _ in range(m):
        x, y = rng.standard_normal(N), rng.standard_normal(N); p = M @ product(x, y, N)
        a = max(a, np.abs(p - product(M @ x, M @ y, N)).max()); b = max(b, np.abs(p - product(M @ y, M @ x, N)).max())
    return float(a), float(b)


def orbit_of_identity(N, k=4, steps=5):
    """J_left^n e₀ for n = 0 … steps, as (index, sign) pairs: the stroke cycle."""
    J = J_left(N, k); v = np.zeros(N); v[0] = 1; out = []
    for _ in range(steps + 1):
        i = int(np.nonzero(v)[0][0]); out.append((i, int(v[i]))); v = J @ v
    return out


# ------------------------------------------------------------ relation of every slot to the axis slot k
def slot_relation(N, j, k=4, alg=None):
    """'zero' j = 0 (identity), 'axis' j = k, else 'anti' (L_j L_k = −L_k L_j: J₂ reverses the circle), 'comm' (commute), or 'deform' (neither)."""
    if j == 0: return "zero"
    if j == k: return "axis"
    A = alg or Alg(N); Lj, Lk = A.Lmat(j), A.Lmat(k); c = Lj @ Lk
    if np.allclose(c, -(Lk @ Lj)): return "anti"
    if np.allclose(c, Lk @ Lj): return "comm"
    return "deform"


def relation_table(N, k=4):
    A = Alg(N); return [slot_relation(N, j, k, A) for j in range(N)]


def rule_class_shares(N, rule, primes, k=4, rng=None):
    tab = relation_table(N, k); slots = rule_slots(rule, primes, N, rng)
    names = ("zero", "axis", "anti", "comm", "deform"); c = {n: 0 for n in names}
    for s in slots: c[tab[int(s)]] += 1
    return {n: c[n] / len(slots) for n in names}


# ------------------------------------------------------------ the Klein group of the functional equation
def klein_check(points):
    """Max error of ξ(s) = ξ(1−s), ξ(s̄) = conj ξ(s), and the same for Λ(s,χ₋₄) = (π/4)^{−(s+1)/2} Γ((s+1)/2) L(s,χ₋₄); with mpmath.  points: complex numbers."""
    import mpmath as mp
    mp.mp.dps = 30
    xi = lambda s: 0.5 * s * (s - 1) * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)
    lam = lambda s: (mp.pi / 4) ** (-(s + 1) / 2) * mp.gamma((s + 1) / 2) * mp.dirichlet(s, [0, 1, 0, -1])
    out = {"xi_fe": 0, "xi_schwarz": 0, "L4_fe": 0, "L4_schwarz": 0}
    for s in points:
        s = mp.mpc(s)
        out["xi_fe"] = max(out["xi_fe"], float(abs(xi(s) - xi(1 - s)))); out["xi_schwarz"] = max(out["xi_schwarz"], float(abs(xi(mp.conj(s)) - mp.conj(xi(s)))))
        out["L4_fe"] = max(out["L4_fe"], float(abs(lam(s) - lam(1 - s)))); out["L4_schwarz"] = max(out["L4_schwarz"], float(abs(lam(mp.conj(s)) - mp.conj(lam(s)))))
    return out


def klein_orbit_size(rho):
    """number of distinct points among {ρ, 1−ρ, ρ̄, 1−ρ̄} (tolerance 1e-9): 2 on the line, 4 off it."""
    pts = [rho, 1 - rho, rho.conjugate(), 1 - rho.conjugate()]; u = []
    for p in pts:
        if not any(abs(p - q) < 1e-9 for q in u): u.append(p)
    return len(u)


# ------------------------------------------------------------ arithmetic J₂: fixed (inert) or swapped (split)
def galois_action(d, p):
    """The conjugation of ℚ(√d) on the primes above p: 'fixed' (inert, (d/p) = −1: p stays irreducible = a hole), 'swapped' (split, (d/p) = +1), 'ramified'."""
    k = kronecker(d, p); return "fixed" if k == -1 else "swapped" if k == 1 else "ramified"


def hole_shares(primes, discs=_DISC):
    """share of primes that are holes of each quadratic level (inert), and the share that are holes of NONE of the first m levels (Frobenius trivial: split completely)."""
    inert = {d: float(np.mean([galois_action(d, int(p)) == "fixed" for p in primes])) for d in discs}
    none = {}
    for m in range(1, len(discs) + 1):
        none[m] = float(np.mean([all(kronecker(d, int(p)) != -1 for d in discs[:m]) for p in primes]))
    return inert, none


def decomposition_check(primes, discs=_DISC[:3]):
    """ζ_K = ∏ over the 2^m characters of L(s,χ_S) for K = ℚ(√d₁,…,√d_m): at p with Frobenius b(p), the number of characters with χ(p) = +1 must be 2^m (b = 0, p splits completely)
    or 2^{m−1} (b ≠ 0, Frobenius of order 2): this is the Euler factor (1−x)^{−2^m} against (1−x²)^{−2^{m−1}}.  Returns (primes checked, failures); ramified primes skipped."""
    m = len(discs); bad = 0; n = 0; subsets = range(1 << m)
    for p in primes:
        p = int(p)
        ks = [kronecker(d, p) for d in discs]
        if 0 in ks: continue
        n += 1; plus = sum(1 for S in subsets if math.prod(ks[j] for j in range(m) if S >> j & 1) == 1)
        b = sum((1 << j) for j in range(m) if ks[j] == -1); want = 2 ** m if b == 0 else 2 ** (m - 1)
        bad += plus != want
    return n, bad
