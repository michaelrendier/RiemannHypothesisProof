# SPDX-License-Identifier: GPL-3.0-only
"""bitflip_test.py — the slot bits computed two ways (exponent = Euler's criterion; reciprocity = Jacobi) for every odd n ≤ 10⁷.   :CALCULATED
Writes experiments/bitflip_test.json.   Run: python3 experiments/bitflip_test.py"""
import itertools, json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, ".."))
from rhp import bitflip as B, nshape as S
X = 10 ** 7; t0 = time.time()
P = S.primes_upto(X); isp = np.zeros(X + 1, bool); isp[P] = True
ns = np.arange(3, X + 1, 2, dtype=np.int64); prime = isp[ns]
agree = np.zeros((len(B.DISC), len(ns)), bool); expbit = np.zeros_like(agree); jacbit = np.zeros_like(agree)
for k, d in enumerate(B.DISC):
    e, j, a = B.two_ways(ns, d); agree[k] = a; expbit[k] = (e == ns - 1); jacbit[k] = (j == -1)
print("computed", round(time.time() - t0)); sys.stdout.flush()
out = dict(N_max=X, bases=list(B.DISC), odd_n=int(len(ns)), primes=int(prime.sum()), composites=int((~prime).sum()))
out["primes_agree_all_7"] = int(agree[:, prime].all(axis=0).sum()); out["primes_total"] = int(prime.sum())
bad_primes = ns[prime & ~agree.all(axis=0)]; out["primes_not_agreeing_first"] = bad_primes[:10].tolist()
comp = ~prime; out["per_base_composites_agreeing"] = {str(d): int(agree[k][comp].sum()) for k, d in enumerate(B.DISC)}
allagree = agree.all(axis=0) & comp; out["composites_agreeing_all_7"] = int(allagree.sum()); out["first_composites_agreeing_all_7"] = ns[allagree][:8].tolist()
# slot-vector agreement: exponent-bit vector equals Jacobi-bit vector (where the residue is ±1) — for composites that agree on every base this is identity
mask = np.zeros(len(ns), np.int64)
for k in range(len(B.DISC)): mask |= (agree[k].astype(np.int64) << k)
cnt = np.bincount(mask[comp], minlength=128)                      # composites by agreement mask
subs = {}
for r in range(1, 8):
    best = None
    for S_ in itertools.combinations(range(7), r):
        sm = sum(1 << k for k in S_); surv = int(sum(cnt[m] for m in range(128) if (m & sm) == sm))
        subs.setdefault(r, []).append((surv, S_))
    subs[r].sort()
out["best_subset_by_size"] = {str(r): dict(survivors=v[0][0], bases=[B.DISC[k] for k in v[0][1]]) for r, v in subs.items()}
pref = []
for r in range(1, 8):
    sm = (1 << r) - 1; pref.append(int(sum(cnt[m] for m in range(128) if (m & sm) == sm)))
out["survivors_after_prefix_in_listed_order"] = pref
for lim in (10 ** 5, 10 ** 6):
    sel = comp & (ns <= lim); c2 = np.bincount(mask[sel], minlength=128); mn = None
    for r in range(1, 8):
        for S_ in itertools.combinations(range(7), r):
            sm = sum(1 << k for k in S_)
            if sum(c2[m] for m in range(128) if (m & sm) == sm) == 0: mn = (r, [B.DISC[k] for k in S_]); break
        if mn: break
    out[f"minimal_subset_zero_survivors_below_{lim}"] = mn
# immediate (non-Euler) failures: residue not ±1 at base 5 (a Fermat-type witness) among composites
e5 = B.modpow(5, (ns - 1) // 2, ns); out["composites_with_residue_not_pm1_base5"] = int((comp & ~((e5 == 1) | (e5 == ns - 1))).sum())
# complementarity: where the two routes DISAGREE, how often is the exponent-side bit 1 vs the Jacobi bit 1
dis = comp & ~agree[3]
out["base5_disagreements_composites"] = int(dis.sum())
json.dump(out, open(os.path.join(HERE, "bitflip_test.json"), "w"), indent=1)
print(out["primes_agree_all_7"], "/", out["primes_total"], "primes agree on all 7;", out["composites_agreeing_all_7"], "composites agree on all 7:", out["first_composites_agreeing_all_7"])
print("per base composites agreeing:", out["per_base_composites_agreeing"]); print("best subsets:", out["best_subset_by_size"]); print("prefix survivors:", pref)
print({k: v for k, v in out.items() if k.startswith("minimal")})
