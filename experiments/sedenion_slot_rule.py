# SPDX-License-Identifier: GPL-3.0-only
"""Sedenion slot-rule test (rhp/sedenion_slots.py).  E1 structure vs arithmetic with a shuffle null; E2 operator Euler product with a random-slot null.  :CALCULATED
Run: python3 experiments/sedenion_slot_rule.py   → experiments/sedenion_slot_rule.json"""
import json, math, os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from rhp import sedenion_slots as S

HERE = os.path.dirname(os.path.abspath(__file__)); out = {}
P = S.primes_upto(5000)
RULES = ("seq0", "seq1", "mod", "char", "charx")
TS = (14.134725, 21.022040, 30.0, 50.0, 100.0)
rng = np.random.default_rng(20261110)

# capacity
cap = {}
for N in (16, 32, 64, 128):
    cap[N] = dict(seq0_last_slotted=int(P[N - 1]), seq0_first_unslotted=int(P[N]), seq1_last_slotted=int(P[N - 2]), seq1_first_unslotted=int(P[N - 1]))
out["capacity"] = cap

# E1  structure ↔ arithmetic: phi correlation of the deform graph with pair features, shuffle null
def pair_features(pr):
    a = pr[:, None]; b = pr[None, :]
    both3 = (a % 4 == 3) & (b % 4 == 3)
    same4 = (a % 4) == (b % 4)
    leg = np.array([[S.kronecker(int(x), int(y)) if x != y else 0 for y in pr] for x in pr]) == -1
    return dict(reciprocity_minus=both3, same_mod4=same4, legendre_minus=leg)
def phi(x, y):
    x = x.astype(float); y = y.astype(float)
    if x.std() == 0 or y.std() == 0: return float("nan")
    return float(np.corrcoef(x, y)[0, 1])
E1 = {}
pr = P[6:206]                                   # 200 primes 17 … (avoids the d_j that divide a prime)
feat = pair_features(pr); iu = np.triu_indices(len(pr), 1)
for N in (16, 32, 64, 128):
    D = S.deform_matrix(N); E1[N] = {}
    for rule in RULES:
        slots = S.rule_slots(rule, pr, N)
        d = D[slots[:, None], slots[None, :]][iu]
        row = {}
        for fname, F in feat.items():
            f = F[iu]; obs = phi(d, f)
            null = []
            for _ in range(1500):
                sh = rng.permutation(slots); null.append(phi(D[sh[:, None], sh[None, :]][iu], f))
            null = np.array(null); null = null[~np.isnan(null)]
            row[fname] = dict(phi=obs, null_mean=float(null.mean()), null_sd=float(null.std()), z=float((obs - null.mean()) / null.std()) if null.std() > 0 else None,
                              p_two_sided=float((np.sum(np.abs(null - null.mean()) >= abs(obs - null.mean())) + 1) / (len(null) + 1)))
        E1[N][rule] = row
out["E1"] = E1
print("E1 done"); sys.stdout.flush()

# E2  operator Euler product
K = 400; pk = P[:K]; E2 = {}
for N in (16, 32, 64, 128):
    A = S.Alg(N); zc = {t: S.complex_path(pk, t) for t in TS}; E2[N] = {}
    def run(slots):
        acc = {k: [] for k in ("ratio_log", "shadow", "entropy", "drift")}
        for t in TS:
            st = S.stats(S.euler_path(slots, pk, t, N, A), zc[t], N)
            acc["ratio_log"].append(abs(math.log(st["ratio"]))); acc["shadow"].append(st["shadow"]); acc["entropy"].append(st["entropy"]); acc["drift"].append(st["drift"])
        return {k: float(np.mean(v)) for k, v in acc.items()}
    ctl = run(S.rule_slots("same", pk, N)); nulls = [run(S.rule_slots("random", pk, N, rng)) for _ in range(300)]
    E2[N]["control_same_slot"] = ctl
    E2[N]["null"] = {k: dict(mean=float(np.mean([r[k] for r in nulls])), sd=float(np.std([r[k] for r in nulls]))) for k in ctl}
    for rule in RULES:
        r = run(S.rule_slots(rule, pk, N)); E2[N][rule] = {k: dict(value=r[k], z=float((r[k] - E2[N]["null"][k]["mean"]) / E2[N]["null"][k]["sd"]),
                                                                  pct=float(np.mean([n[k] <= r[k] for n in nulls]))) for k in r}
    print("E2 N", N, "done"); sys.stdout.flush()
out["E2"] = E2
out["slot0_primes"] = {N: {rule: int((S.rule_slots(rule, pk, N) == 0).sum()) for rule in RULES} for N in (16, 32, 64, 128)}
json.dump(out, open(os.path.join(HERE, "sedenion_slot_rule.json"), "w"), indent=1)
