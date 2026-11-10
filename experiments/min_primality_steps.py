# SPDX-License-Identifier: GPL-3.0-only
"""min_primality_steps.py — how few steps does an exact primality test need, on the range the 65 wheels cover and beyond?   :CALCULATED
Exhaustive over odd n < 10⁷: composites passing Fermat base 2, strong (Miller–Rabin) base 2, bases 2+3, bases 2+3+5; against the wheels π(√N) a sieve needs.
Writes experiments/min_primality_steps.json."""
import json, os, time
import numpy as np
X = 10_000_000
sv = np.ones(X + 1, bool); sv[:2] = False
for p in range(2, int(X ** .5) + 1):
    if sv[p]: sv[p * p::p] = False
def sprp(n, a):
    if n % a == 0: return n == a
    d = n - 1; s = 0
    while d % 2 == 0: d //= 2; s += 1
    x = pow(a, d, n)
    if x == 1 or x == n - 1: return True
    for _ in range(s - 1):
        x = x * x % n
        if x == n - 1: return True
    return False
odd = [int(n) for n in range(3, X + 1, 2)]
out = {}
fermat = [n for n in odd if pow(2, n - 1, n) == 1 and not sv[n]]
s2 = [n for n in odd if not sv[n] and sprp(n, 2)]
s23 = [n for n in s2 if sprp(n, 3)]; s235 = [n for n in s23 if sprp(n, 5)]
out["composites_passing_below_1e7"] = dict(fermat_base2=len(fermat), strong_base2=len(s2), strong_2_3=len(s23), strong_2_3_5=len(s235))
out["below_100489"] = dict(fermat_base2=sum(n < 100489 for n in fermat), strong_base2=sum(n < 100489 for n in s2), strong_2_3=sum(n < 100489 for n in s23))
out["first_exceptions"] = dict(fermat_base2=fermat[:6], strong_base2=s2[:6], strong_2_3=s23[:4], strong_2_3_5=s235[:3])
out["strong_base2_below_100489"] = [n for n in s2 if n < 100489]
out["wheels_needed_pi_sqrtN"] = {str(N): int(sv[:int(N ** .5) + 1].sum()) for N in (10 ** 5, 10 ** 6, 1_373_653, 10 ** 7)}
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "min_primality_steps.json"), "w"), indent=1)
print(out["composites_passing_below_1e7"], out["below_100489"])
