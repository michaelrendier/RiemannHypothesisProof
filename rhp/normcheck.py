# SPDX-License-Identifier: GPL-3.0-only
"""Norm check at the interfaces of the Cayley–Dickson tower T_1 ⊂ T_2 ⊂ T_4 ⊂ … ⊂ T_N.   Everything :CALCULATED (floating point) or exact (integers) as marked.

The composition law |xy| = |x||y| holds at N = 1, 2, 4, 8 (Hurwitz, ESTABLISHED) and fails from N = 16.  This module measures it at every doubling and splits the
failure by INTERFACE: a doubling T_{N/2} → T_N cuts T_N into a lower half L (indices < N/2, the old algebra) and an upper half U (indices ≥ N/2, the new directions).
  L×L  stays inside T_{N/2}: the old law, inherited exactly
  L×U, U×U  the products that only exist after the doubling
The defect  D(x, y) = |xy|² − |x|²|y|  (squared both sides) is zero where the law holds.  Zero divisors are the points where D = −|x|²|y|² (xy = 0).
Also exact (integer) checks on sums of two basis units, the objects of rhp.zdtree: |e_a ± e_b|² = 2, and |xy|² takes only a few values.
"""
from __future__ import annotations
import numpy as np
from .sedenion_slots import sign_table


def mul(x, y, N):
    """Product in T_N of batches x, y (shape (m, N) or (N,)), by the sign table: out[i⊕j] += S[i,j] x_i y_j."""
    x = np.atleast_2d(x); y = np.atleast_2d(y); S = sign_table(N); idx = np.arange(N)
    out = np.zeros(np.broadcast_shapes(x.shape, y.shape), dtype=np.result_type(x, y))
    for i in range(N):
        out[:, idx ^ i] += x[:, i:i + 1] * S[i][None, :] * y
    return out


def norm2(x):
    return (np.atleast_2d(x) ** 2).sum(axis=1)


def defect(x, y, N):
    """(|xy|² − |x|²|y|²) / (|x|²|y|²): zero where the composition law holds, −1 at a zero divisor."""
    p = norm2(mul(x, y, N)); q = norm2(x) * norm2(y)
    return p / q - 1.0


def interface_samples(N, m, rng):
    """m random pairs for each interface class of the doubling T_{N/2} → T_N: LL, LU, UL, UU (Gaussian entries, support restricted to the half)."""
    h = N // 2
    def half(upper):
        v = rng.standard_normal((m, N))
        if upper: v[:, :h] = 0
        else: v[:, h:] = 0
        return v
    return {"LL": (half(0), half(0)), "LU": (half(0), half(1)), "UL": (half(1), half(0)), "UU": (half(1), half(1))}


def level_report(N, m=2000, seed=0):
    """Defect statistics at T_N: all (full random vectors) and per interface class (N ≥ 2)."""
    rng = np.random.default_rng(seed); out = {}
    x = rng.standard_normal((m, N)); y = rng.standard_normal((m, N))
    d = defect(x, y, N); out["all"] = dict(max_abs=float(np.abs(d).max()), mean_abs=float(np.abs(d).mean()), min=float(d.min()), max=float(d.max()))
    if N >= 2:
        for k, (a, b) in interface_samples(N, m, rng).items():
            d = defect(a, b, N); out[k] = dict(max_abs=float(np.abs(d).max()), mean_abs=float(np.abs(d).mean()))
    return out


def basis_pair_sums(N):
    """Exact: for every pair (a<b) of imaginary units and every pair (c<d), the values of |(e_a ± e_b)(e_c ± e_d)|² (integers, |x|²|y|² = 4).
    Returns {value: count}.  A value 0 = a zero divisor; a value 4 = norm preserved; 8 = norm doubled (the law broken upward)."""
    S = sign_table(N); vals = {}
    idxs = [(a, b, s) for a in range(1, N) for b in range(a + 1, N) for s in (1, -1)]
    X = np.zeros((len(idxs), N), np.int64)
    for r, (a, b, s) in enumerate(idxs): X[r, a] = 1; X[r, b] = s
    for r in range(len(idxs)):
        P = mul(np.broadcast_to(X[r], X.shape).astype(np.int64), X, N)
        for v in norm2(P).tolist(): vals[int(v)] = vals.get(int(v), 0) + 1
    return vals


def zd_from_norm(N):
    """Number of unordered {±}-pairs (x, y) of basis-pair sums with xy = 0 (norm_2 = 0), from basis_pair_sums: each zero is a pair of assessor-type elements."""
    return basis_pair_sums(N).get(0, 0)


def zero_product_halves(N):
    """Exact: classify every zero product (e_a ± e_b)(e_c ± e_d) = 0 of T_N by which half (L = index < N/2, U = index ≥ N/2) each unit lies in.
    Returns {(halves of x, halves of y): count}.  At T_16 all 336 are (L,U | L,U): every zero divisor mixes the old and the new directions."""
    h = N // 2; units = [(a, b, s) for a in range(1, N) for b in range(a + 1, N) for s in (1, -1)]
    X = np.zeros((len(units), N), np.int64)
    for r, (a, b, s) in enumerate(units): X[r, a] = 1; X[r, b] = s
    cnt = {}
    for r, (a, b, s) in enumerate(units):
        z = np.nonzero(norm2(mul(np.broadcast_to(X[r], X.shape).astype(np.int64), X, N)) == 0)[0]
        hx = "".join(sorted("L" if i < h else "U" for i in (a, b)))
        for q in z:
            hy = "".join(sorted("L" if i < h else "U" for i in units[q][:2])); cnt[(hx, hy)] = cnt.get((hx, hy), 0) + 1
    return cnt
