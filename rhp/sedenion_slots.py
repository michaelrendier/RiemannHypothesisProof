# SPDX-License-Identifier: GPL-3.0-only
"""Sedenion slot-rule test: which prime sits on which slot of T_N (N = 16, 32, 64, 128)?   Everything :CALCULATED.

Machinery (exact integers): the Cayley–Dickson sign table  e_i e_j = S[i,j] e_{i⊕j}  built by doubling, convention (a₁,a₂)(b₁,b₂) = (a₁b₁ − b₂*a₂, b₂a₁ + a₂b₁*)
(the one in SedenionSpectralRelativity/emerger_spectrum.cd_mul; checked against it at N = 16).  L_k = left multiplication by e_k, a signed permutation with L_k² = −I.

The prime-slot idea (wiki 130): the Euler factor of p on σ = ½ is the circle  (1 − p^{-1/2} e^{-it ln p})^{-1}; put the circle's i on slot k, i.e. i → L_k.  Closed form
    F_p = (αI − βL_k)/(α²+β²),  α = 1 − p^{-1/2} cos θ,  β = p^{-1/2} sin θ,  θ = t ln p,
so the sedenion partial Euler product is v_K = F_{p_K} … F_{p_1} e₀, an O(N) update per prime.  All primes on ONE slot = the complex partial product exactly (control).
A slot that is the identity (k = 0) has no circle: that prime's factor is the real scalar (1 − p^{-1/2}(cos θ − sin θ))^{-1} (stated, not hidden).

Rules (a rule is a map  prime → slot ∈ [0, N)):
  same    every prime on slot 1                                 (control: complex ζ)
  seq0    slot = (index − 1) mod N,  index 1-based              (wiki 130 table: slot 0 ↔ 2, slot 15 ↔ 53; wraps past N)
  seq1    slot = ((index − 1) mod (N − 1)) + 1                  (k-th prime on slot k, k = 1..N−1; wraps)
  mod     slot = p mod N                                        (residue class; the slot 0 is never used)
  char    slot = Σ_j 2^j [ (d_j / p) = −1 ],  d_j = −4, 8, −3, 5, −7, −11, 13   (Kronecker symbols: slot(pq) = slot(p) ⊕ slot(q), a homomorphism; m = log₂N bits)
  random  uniform over slots 1..N−1 (the null; many draws)
"""
from __future__ import annotations
import math
import numpy as np

# ---------------------------------------------------------------- Cayley–Dickson sign table
_S_CACHE: dict[int, np.ndarray] = {}


def sign_table(N: int) -> np.ndarray:
    """S[i,j] = ±1 with e_i e_j = S[i,j] e_{i⊕j}, N a power of 2."""
    if N in _S_CACHE: return _S_CACHE[N]
    if N == 1: S = np.ones((1, 1), np.int8)
    else:
        n = N // 2; Sn = sign_table(n); S = np.empty((N, N), np.int8)
        cd = np.where(np.arange(n) == 0, 1, -1).astype(np.int8)               # sign of e_d* = cd[d]·e_d
        S[:n, :n] = Sn                                                          # (e_a,0)(e_c,0)
        S[:n, n:] = Sn.T                                                        # (e_a,0)(0,e_d) = (0, e_d e_a):    S[a, n+d] = Sn[d,a]
        S[n:, :n] = Sn * cd[None, :]                                            # (0,e_b)(e_c,0) = (0, e_b e_c*):   S[n+b, c] = Sn[b,c] cd[c]
        S[n:, n:] = -(Sn.T * cd[None, :])   # (0,e_b)(0,e_d) = (−e_d* e_b, 0): S[n+b,n+d] = −cd[d] Sn[d,b]
    _S_CACHE[N] = S
    return S


def basis_product(i: int, j: int, N: int):
    s = sign_table(N)[i, j]
    return int(s), i ^ j


class Alg:
    """Left multiplications of T_N as signed permutations (O(N) action)."""
    def __init__(self, N):
        self.N = N; self.S = sign_table(N); self.idx = np.arange(N)
        self.perm = [self.idx ^ k for k in range(N)]

    def L(self, k, v):
        out = np.empty_like(v); out[self.perm[k]] = self.S[k] * v; return out

    def Lmat(self, k):
        M = np.zeros((self.N, self.N)); M[self.perm[k], self.idx] = self.S[k]; return M


# ---------------------------------------------------------------- structure (which pairs deform)
def deform_matrix(N: int):
    """D[i,j] = True iff L_i L_j + L_j L_i ≠ 0 (imaginary slots i ≠ j, both ≥ 1); the 'Assessor'-type pairs (42 in T_16)."""
    A = Alg(N); Ls = [A.Lmat(k) for k in range(N)]
    D = np.zeros((N, N), bool)
    for i in range(1, N):
        for j in range(i + 1, N):
            ac = Ls[i] @ Ls[j] + Ls[j] @ Ls[i]
            D[i, j] = D[j, i] = bool(np.abs(ac).max() > 1e-9)
    return D


# ---------------------------------------------------------------- primes and rules
def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for p in range(2, int(n ** .5) + 1):
        if s[p]: s[p * p::p] = False
    return np.nonzero(s)[0]


_DISC = (-4, 8, -3, 5, -7, -11, 13)


def kronecker(d: int, p: int) -> int:
    """(d/p) for prime p: Legendre for odd p, Kronecker at 2, 0 when p | d."""
    if p == 2:
        if d % 2 == 0: return 0
        return 1 if d % 8 in (1, 7) else -1
    if d % p == 0: return 0
    r = pow(d % p, (p - 1) // 2, p)
    return 1 if r == 1 else -1


def rule_slots(rule: str, primes, N: int, rng=None):
    primes = np.asarray(primes); n = len(primes); idx = np.arange(1, n + 1)
    if rule == "same": return np.ones(n, int)
    if rule == "seq0": return (idx - 1) % N
    if rule == "seq1": return (idx - 1) % (N - 1) + 1
    if rule == "mod": return primes % N
    if rule == "char":
        m = int(math.log2(N)); out = np.zeros(n, int)
        for j, d in enumerate(_DISC[:m]):
            out += np.array([1 << j if kronecker(d, int(p)) == -1 else 0 for p in primes])
        return out
    if rule == "charx":
        m = int(math.log2(N)) - 1; out = np.zeros(n, int)
        for j, d in enumerate(_DISC[:m]):
            out += np.array([1 << j if kronecker(d, int(p)) == -1 else 0 for p in primes])
        out[out == 0] = N // 2
        return out
    if rule == "random": return rng.integers(1, N, n)
    raise ValueError(rule)


# ---------------------------------------------------------------- the operator Euler product
def euler_path(slots, primes, t, N, alg: Alg | None = None):
    """v_k = F_{p_k} v_{k−1}, v₀ = e₀ (sedenion-style partial Euler product, circle's i → L_slot).  Returns the array of v_k (K × N)."""
    alg = alg or Alg(N); K = len(primes)
    v = np.zeros(N); v[0] = 1.0; out = np.empty((K, N))
    for k, (p, s) in enumerate(zip(primes, slots)):
        a = p ** -0.5; th = t * math.log(p); c, sn = math.cos(th), math.sin(th)
        if s == 0:
            v = v / (1.0 - a * (c - sn))                                      # identity slot: real scalar factor (no circle)
        else:
            al, be = 1.0 - a * c, a * sn
            v = (al * v - be * alg.L(int(s), v)) / (al * al + be * be)
        out[k] = v
    return out


def complex_path(primes, t):
    """The ordinary complex partial Euler product ∏(1 − p^{-1/2−it})^{-1}, running."""
    z = np.cumprod(1.0 / (1.0 - np.asarray(primes, float) ** (-0.5) * np.exp(-1j * t * np.log(np.asarray(primes, float)))))
    return z


def stats(path, zc, N):
    """norm ratio |v_K|/|P_K|, scalar-shadow error |v₀ − Re P_K|/|P_K| (at the last K), occupancy entropy of the slot energies (normalised), and the mean over K of |log(|v|/|P|)|."""
    v = path[-1]; P = zc[-1]
    e = v * v / (v @ v); H = -(e[e > 0] * np.log(e[e > 0])).sum() / math.log(N)
    nr = np.linalg.norm(path, axis=1) / np.abs(zc)
    return dict(ratio=float(nr[-1]), shadow=float(abs(v[0] - P.real) / abs(P)), entropy=float(H), drift=float(np.mean(np.abs(np.log(nr)))))
