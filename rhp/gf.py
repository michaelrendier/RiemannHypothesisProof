# SPDX-License-Identifier: GPL-3.0-only
"""Small finite fields GF(p^k) (q ≤ ~5000) and the exact count of  x^n + y^n = 1  over GF(q), by log tables.   :CALCULATED (exact integers)."""
from __future__ import annotations
import itertools
import numpy as np


def _polymulmod(a, b, mod, p):
    k = len(mod) - 1; res = [0] * (2 * k - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b): res[i + j] = (res[i + j] + x * y) % p
    for d in range(2 * k - 2, k - 1, -1):                                   # reduce by the monic modulus (list low → high, mod[k] = 1)
        c = res[d]
        if c:
            for i in range(k + 1): res[d - k + i] = (res[d - k + i] - c * mod[i]) % p
    return res[:k]


def _find_field(p, k):
    """(exp table as list of digit-tuples, log dict) for GF(p^k) with a primitive element."""
    q = p ** k
    for coeffs in itertools.product(range(p), repeat=k):                    # monic modulus x^k + c_{k−1}x^{k−1}+…+c_0
        mod = list(coeffs) + [1]
        if k > 1 and mod[0] == 0: continue
        g = [0] * k; g[1 % k] = 1 if k > 1 else 0
        if k == 1:
            for gg in range(2, p):                                          # prime field: search a primitive root
                x, seen = 1, 0
                for e in range(1, p):
                    x = x * gg % p
                    if x == 1: seen = e; break
                if seen == p - 1: return [(pow(gg, e, p),) for e in range(p - 1)], {(pow(gg, e, p),): e for e in range(p - 1)}
            return [(1,)], {(1,): 0}
        x = [1] + [0] * (k - 1); exp = []; seen = set(); ok = True
        for e in range(q - 1):
            t = tuple(x)
            if t in seen: ok = False; break
            seen.add(t); exp.append(t); x = _polymulmod(x, g, mod, p)
        if ok and tuple(x) == tuple([1] + [0] * (k - 1)) and len(exp) == q - 1: return exp, {t: i for i, t in enumerate(exp)}
    raise RuntimeError("no primitive polynomial found")


def count_fermat_affine(p, k, n):
    """#{(x, y) ∈ GF(p^k)² : x^n + y^n = 1}."""
    q = p ** k; exp, log = _find_field(p, k)
    digs = np.array([[(i // p ** j) % p for j in range(k)] for i in range(q)], np.int64); code = lambda d: (d * (p ** np.arange(k))).sum(axis=-1)
    idx = {tuple(d): i for i, d in enumerate(digs.tolist())}
    pw = np.zeros(q, np.int64)                                               # x^n as an element index (0 ↦ 0)
    for e, t in enumerate(exp): pw[idx[t]] = idx[exp[(n * e) % (q - 1)]]
    cnt = np.bincount(pw, minlength=q)
    one = np.zeros(k, np.int64); one[0] = 1; om = code((one[None, :] - digs) % p)   # 1 − u for every u
    return int((cnt * cnt[om]).sum())
