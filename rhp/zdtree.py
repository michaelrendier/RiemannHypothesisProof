# SPDX-License-Identifier: GPL-3.0-only
"""The Zero-Divisor Tree of the Cayley–Dickson tower, rebuilt in RHP from its own sign tables (rhp.sedenion_slots) and checked against AbrikosovTree/render/zd_propagation.py.   :CALCULATED
Definition (De Marrais / Moreno assessors): for imaginary a < b with δ = a ⊕ b, the element e_a + e_b is a zero divisor iff some e_c ± e_{c⊕δ} annihilates it.
Descent (zd_propagation.py): the XOR differences δ are the nodes; the doubling sends δ → {δ, δ + 2^k}."""
from __future__ import annotations
import numpy as np
from .sedenion_slots import sign_table


def zd_incidence(k):
    """boolean matrix Z[a, b] (a, b imaginary slots of the 2^k-ion): e_a + e_b is a zero divisor.  Also returns per-δ counts (zd pairs, all pairs)."""
    n = 2 ** k; S = sign_table(n).astype(np.int64); Z = np.zeros((n, n), bool)
    per = np.zeros(n, int); tot = np.zeros(n, int); cs = np.arange(1, n)
    for d in range(1, n):
        cc = cs[cs != d]; dd = cc ^ d
        for a in range(1, n):
            b = a ^ d
            if not a < b: continue
            tot[d] += 1
            plus = (S[a, cc] == -S[b, dd]) & (S[b, cc] == -S[a, dd])
            minus = (S[a, cc] == S[b, dd]) & (S[b, cc] == S[a, dd])
            if (plus | minus).any(): Z[a, b] = Z[b, a] = True; per[d] += 1
    return Z, per, tot


def classes(per, tot):
    """(ZD-carrying δ, clean δ) — clean = a difference class with pairs but no zero divisor."""
    carry = [int(d) for d in range(1, len(per)) if per[d] > 0]; clean = [int(d) for d in range(1, len(per)) if tot[d] > 0 and per[d] == 0]
    return carry, clean
