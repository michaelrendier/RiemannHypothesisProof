# SPDX-License-Identifier: GPL-3.0-only
"""ping_extend.py — the per-integer ping (ψ_K jump across [n−½, n+½] against Λ(n)) for n = 2 … NMAX with K = 3066 and 25,000 zero-tones.  :CALCULATED
Writes experiments/zeros25k/ping_K<K>_n<NMAX>.npy (float32).  Run:  python3 experiments/ping_extend.py [NMAX]"""
import os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import zeros25k_rerun as Z
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
g = np.array([float(l.split("\t")[1]) for l in open(Z.HP).read().split("\n")[1:] if l.strip()])
ns = np.arange(2, NMAX + 1)
for K in (3066, 25000):
    t0 = time.time(); rj = Z.jump_recovered(ns, g[:K], chunk=32)
    np.save(os.path.join(Z.OUT, f"ping_K{K}_n{NMAX}.npy"), rj.astype(np.float32))
    Z.log("ping_extend K", K, "n ≤", NMAX, "secs", round(time.time() - t0))
Z.log("ping_extend DONE")
