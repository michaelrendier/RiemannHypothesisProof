"""Diagonal prime quadratics as phase offsets, and the prime-phase torus.  Writes diagonal_quadratics.json.

Expectations written before the run:
  Q1  for f(n) = n^2 + b n + c with discriminant D = b^2 - 4c, the number of roots mod an odd prime p not dividing D is w(p) = 1 + (D/p) (Legendre symbol).
  Q2  a shift n -> n + d (a phase offset in n) leaves D unchanged; the prime counts along f(n+d) for n in a fixed window agree with the counts along f(n) within sampling noise.
  Q3  the Bateman-Horn prediction C(f) * sum 1/ln f(n), with C(f) = prod_p (1 - w(p)/p)/(1 - 1/p), tracks the observed counts across many quadratics.
  Q4  the pair of phases (t ln 2, t ln 3) mod 2 pi, sampled along the critical-line path, covers a 32 x 32 grid of the torus (Kronecker: ln 2 / ln 3 irrational)."""
import json, math, pathlib, sys
import numpy as np
root = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(root.parent))
LIM = 40_000_000
S = np.ones(LIM + 1, bool); S[:2] = False
for i in range(2, int(LIM ** 0.5) + 1):
    if S[i]: S[i * i::i] = False
PR = np.nonzero(S)[0]
def legendre(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1
R = {}
# Q1
bad = 0; tested = 0; rng = np.random.default_rng(3)
for _ in range(60):
    b = int(rng.integers(-20, 21)); c = int(rng.integers(1, 80)); D = b * b - 4 * c
    for p in PR[1:40]:
        p = int(p)
        if D % p == 0: continue
        w = sum(1 for n in range(p) if (n * n + b * n + c) % p == 0); tested += 1; bad += (w != 1 + legendre(D, p))
R['Q1'] = dict(tested=tested, mismatches=bad)
# Q2: Euler's polynomial n^2+n+41 (D = -163) and its shifts
N = 5000
def count(b, c, d=0, N=N):
    n = np.arange(d, d + N, dtype=np.int64); v = n * n + b * n + c; v = v[(v > 1) & (v <= LIM)]; return int(S[v].sum())
shifts = []; base = count(1, 41)
for d in range(0, 60, 5):
    b2 = 1 + 2 * d * 0                                                   # f(n + d) = n^2 + (1 + 2d) n + (d^2 + d + 41), window n in [0, N)
    shifts.append(dict(d=d, b=1 + 2 * d, c=d * d + d + 41, D=(1 + 2 * d) ** 2 - 4 * (d * d + d + 41), count=count(1 + 2 * d, d * d + d + 41)))
cs = np.array([s['count'] for s in shifts]); R['Q2'] = dict(D_all_equal=len({s['D'] for s in shifts}) == 1, D=shifts[0]['D'], base_count=base, counts=cs.tolist(), mean=float(cs.mean()), sd=float(cs.std()), poisson_sd=float(math.sqrt(cs.mean())))
# Q3: Bateman-Horn
P_ = PR[PR < 200000]
def C_of(b, c):
    D = b * b - 4 * c; prod = 1.0
    for p in P_:
        p = int(p)
        if p == 2: w = sum(1 for n in range(2) if (n * n + b * n + c) % 2 == 0)
        elif D % p == 0: w = 1
        else: w = 1 + legendre(D, p)
        prod *= (1 - w / p) / (1 - 1 / p)
    return prod
obs = []; pred = []; rows = []
for _ in range(120):
    b = int(rng.integers(-30, 31)); c = int(rng.integers(2, 200))
    n = np.arange(0, N, dtype=np.int64); v = n * n + b * n + c; ok = v > 1
    if not ok.all() and ok.sum() < 0.5 * N: continue
    o = int(S[v[ok]].sum()); pr = C_of(b, c) * float(np.sum(1.0 / np.log(v[ok].astype(float)))); obs.append(o); pred.append(pr); rows.append((b, c, b * b - 4 * c, o, round(pr, 1)))
obs = np.array(obs); pred = np.array(pred); R['Q3'] = dict(n=len(obs), pearson_r=float(np.corrcoef(obs, pred)[0, 1]), n_zero_prediction=int((pred <= 0).sum()), zero_prediction_observed_primes=int(obs[pred <= 0].sum()), median_ratio=float(np.median(obs[pred > 0] / pred[pred > 0])), ratio_p5=float(np.percentile(obs[pred > 0] / pred[pred > 0], 5)), ratio_p95=float(np.percentile(obs[pred > 0] / pred[pred > 0], 95)),
     euler_41=dict(observed=count(1, 41), predicted=round(C_of(1, 41) * float(np.sum(1.0 / np.log((np.arange(N) ** 2 + np.arange(N) + 41).astype(float)))), 1)), examples=rows[:6])
# Q4
from ValaQuenta.modules.zeta_uns_path import maths as M
ts = M.sample_by_phase(10.0, 25000, 16, jitter=True); a = np.mod(ts * math.log(2), 2 * math.pi); b = np.mod(ts * math.log(3), 2 * math.pi)
for G in (16, 32, 64): H = np.histogram2d(a, b, bins=G, range=[[0, 2 * math.pi]] * 2)[0]; R[f'Q4_grid{G}'] = float((H > 0).mean())
json.dump(R, open(root / 'experiments/diagonal_quadratics.json', 'w'), indent=1); print(json.dumps({k: (v if k != 'Q3' else {kk: vv for kk, vv in v.items() if kk != 'examples'}) for k, v in R.items()}, indent=1))
