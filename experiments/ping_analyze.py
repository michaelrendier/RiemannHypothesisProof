# SPDX-License-Identifier: GPL-3.0-only
"""ping_analyze.py — how well the zero-sum jump separates primes from non-prime-powers, by block of n, and stragglers against other composites.  :CALCULATED
Reads experiments/zeros25k/ping_K{3066,25000}_n200000.npy (experiments/ping_extend.py).  Writes ping_accuracy_n200000.json.
Balanced accuracy = best single threshold, mean of true-positive and true-negative rate (primes vs non-prime-powers); AUC = P(prime scores above a non-prime-power).
Stragglers: composites with no prime factor ≤ 313 (the 65 wheels), n ≥ 317² = 100489.  (An earlier fixed-threshold 'false alarm' rate ½ln2 overstated the loss at large n: at large n every integer reads ≈ 1, the smooth slope.)"""
import json, math, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "zeros25k"); NMAX = 200000
ns = np.arange(2, NMAX + 1); spf = np.zeros(NMAX + 1, int)
for p in range(2, NMAX + 1):
    if spf[p] == 0: spf[p::p][spf[p::p] == 0] = p
s = spf[2:]; lam = np.zeros(NMAX + 1)
for p in range(2, NMAX + 1):
    if spf[p] == p:
        q = p
        while q <= NMAX: lam[q] = math.log(p); q *= p
lam = lam[2:]; prime = s == ns; plain = lam == 0; strag = plain & (s > 313); other = plain & (s <= 313) & ~prime
def acc(rj, m):
    P = m & prime; Q = m & plain; x = np.concatenate([rj[P], rj[Q]]); y = np.concatenate([np.ones(P.sum()), np.zeros(Q.sum())]); o = np.argsort(x); ys = y[o]
    c = np.concatenate([[0], np.cumsum(ys)])[:-1]; npos = ys.sum(); nneg = len(ys) - npos
    tpr = (npos - c) / npos; fpr = (nneg - (np.arange(len(ys)) - c)) / nneg
    r = np.empty(len(x)); r[o] = np.arange(1, len(x) + 1)
    return float(np.max(.5 * (tpr + 1 - fpr))), float((r[y == 1].sum() - npos * (npos + 1) / 2) / (npos * nneg))
edges = [2, 100, 300, 1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 12000, 15000, 20000, 30000, 50000, 100000, 200001]
res = {}
for K in (3066, 25000):
    rj = np.load(os.path.join(OUT, f"ping_K{K}_n{NMAX}.npy")).astype(float); res[K] = dict(blocks={}, stragglers={})
    for a, b in zip(edges[:-1], edges[1:]):
        m = (ns >= a) & (ns < b); ba, auc = acc(rj, m); res[K]["blocks"][f"{a}-{b - 1}"] = dict(balanced_acc=ba, auc=auc, prime_mean=float(rj[m & prime].mean()), nonpp_mean=float(rj[m & plain].mean()), nonpp_sd=float(rj[m & plain].std()))
    for a, b in ((100489, 150000), (150000, 200001)):
        m = (ns >= a) & (ns < b); S = m & strag; O = m & other
        d = (rj[S].mean() - rj[O].mean()) / math.sqrt(.5 * (rj[S].var() + rj[O].var()))
        res[K]["stragglers"][f"{a}-{b - 1}"] = dict(n_strag=int(S.sum()), n_other=int(O.sum()), strag_mean=float(rj[S].mean()), other_mean=float(rj[O].mean()), cohen_d=float(d))
    mm = plain & ~prime; res[K]["corr_rj_ln_spf_composites"] = float(np.corrcoef(rj[mm], np.log(s[mm]))[0, 1])
json.dump(res, open(os.path.join(OUT, f"ping_accuracy_n{NMAX}.json"), "w"), indent=1)
for K, r in res.items():
    print("K", K, "corr(rj, ln spf) over composites:", round(r["corr_rj_ln_spf_composites"], 3), "| stragglers:", {k: (round(v["strag_mean"], 2), round(v["other_mean"], 2), round(v["cohen_d"], 2)) for k, v in r["stragglers"].items()})
    print("   ", {k: round(v["balanced_acc"], 3) for k, v in r["blocks"].items()})
