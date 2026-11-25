"""Block flow of the zeta path (100 zeros at a time) and the flow of the primes through the zeros.  Writes block_flow_and_prime_flow.json.
Expectations written before the run:
  B1  each block of 100 zeros winds the same number of turns (theta advances pi per zero), so its swept AREA does not halve; what can shrink block by block is the NEW area a block adds to what the earlier blocks already covered
      (a coverage that saturates), and the t-span of a block (about 100 * 2 pi / ln(t / 2 pi)).
  L1  Landau (1911): sum_{0 < gamma <= T} cos(gamma ln n) = -(T / 2 pi) Lambda(n) / sqrt(n) + O(log T) for prime powers n, i.e. the zero spectrum has a spike of weight -Lambda(n)/sqrt(n) per (T/2pi) at x = ln n.
  S1  the fluctuation of the zero count S(T) = N(T) - theta(T)/pi - 1 is a sum of prime waves:  S(T) ~ -(1/pi) sum_p sum_k sin(k T ln p) / (k p^(k/2));  regressing S on sin/cos(T ln n) recovers amplitudes -1/(pi k p^(k/2))."""
import json, math, pathlib, sys
import numpy as np
root = pathlib.Path(__file__).resolve().parents[1]; sys.path.insert(0, str(root.parent))
from ValaQuenta.modules.zeta_uns_path import maths as M
R = {}
g = np.array([float(l.split()[1]) for l in open(root / 'data/zeros_25000_hp.tsv') if l[0].isdigit()])
ts = M.sample_by_phase(10.0, 25000, 16, jitter=True); z = M.zeta_line(ts); dt = np.gradient(ts)
# B1
edges = np.concatenate([[10.0], g[99::100]])[:251]; B = []
GN = 40; seen = np.zeros((GN, GN), bool); RNG = [[-4, 6], [-5, 5]]
for k in range(250):
    lo, hi = (g[100 * k - 1] if k else 10.0), g[100 * (k + 1) - 1]; m = (ts >= lo) & (ts < hi); zz = z[m]; w = dt[m] / dt[m].sum()
    X = np.stack([zz.real, zz.imag], 1); mu = (X * w[:, None]).sum(0); C = ((X - mu).T * w) @ (X - mu); ev, evec = np.linalg.eigh(C)
    H, _, _ = np.histogram2d(zz.real, zz.imag, bins=GN, range=RNG); cell = H > 0; new = int((cell & ~seen).sum()); seen |= cell
    B.append(dict(k=k + 1, t0=float(lo), t1=float(hi), span=float(hi - lo), turns=float((M.theta(hi) - M.theta(max(lo, 10.0))) / (2 * math.pi)), mean_re=float(mu[0]), mean_im=float(mu[1]), eig_small=float(ev[0]), eig_big=float(ev[1]),
                  axis_deg=float(np.degrees(np.arctan2(evec[1, 1], evec[0, 1]))), re_neg=float(np.sum(w[zz.real < 0])), cells=int(cell.sum()), new_cells=new, cum_cells=int(seen.sum())))
R['B1_first10'] = B[:10]; R['B1_last3'] = B[-3:]
sp = np.array([b['span'] for b in B]); nw = np.array([b['new_cells'] for b in B]); R['B1_summary'] = dict(span_first=float(sp[0]), span_k2=float(sp[1]), span_last=float(sp[-1]), span_pred_last=float(100 * 2 * math.pi / math.log(g[-1] / (2 * math.pi))),
    turns_min=min(b['turns'] for b in B[1:]), turns_max=max(b['turns'] for b in B[1:]), new_cells_first10=nw[:10].tolist(), new_ratio_first8=[float(nw[i + 1] / max(nw[i], 1)) for i in range(8)], total_cells=int(seen.sum()), cells_per_block_median=float(np.median([b['cells'] for b in B])))
# L1
def Lam(n):
    for p in M.primes_upto(int(n) + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return math.log(p) if m == 1 else 0.0
    return 0.0
rows = []
for n in range(2, 60):
    if Lam(n) == 0: continue
    pred = -Lam(n) / math.sqrt(n); r = []
    for nz in (100, 1000, 5000, 25000):
        T = g[nz - 1]; S = float(np.sum(np.cos(g[:nz] * math.log(n)))); r.append(S / ((T / (2 * math.pi)) * pred))
    rows.append(dict(n=n, pred_weight=pred, ratio_at_100_1000_5000_25000=r))
R['L1'] = rows; R['L1_mean_abs_dev_from_1_at_25000'] = float(np.mean([abs(x['ratio_at_100_1000_5000_25000'][3] - 1) for x in rows])); R['L1_mean_abs_dev_from_1_at_100'] = float(np.mean([abs(x['ratio_at_100_1000_5000_25000'][0] - 1) for x in rows]))
# S1
mid = 0.5 * (g[1:] + g[:-1]); Nn = np.arange(1, len(g)); Sfl = Nn - M.theta(mid) / math.pi - 1.0
pp = [(p, k) for p in M.primes_upto(200) for k in range(1, 4) if p ** k <= 200]; cols = []; names = []
for p, k in pp: cols += [np.sin(mid * k * math.log(p)), np.cos(mid * k * math.log(p))]; names += [(int(p), k, 'sin'), (int(p), k, 'cos')]
A = np.stack(cols, 1); coef, *_ = np.linalg.lstsq(A, Sfl, rcond=None); res = Sfl - A @ coef
R['S1'] = dict(var_S=float(np.var(Sfl)), var_pred=float(math.log(math.log(g[-1])) / (2 * math.pi ** 2)), R2=float(1 - np.var(res) / np.var(Sfl)), n_terms=len(pp),
    fit=[dict(p=p, k=k, fitted_sin=float(coef[2 * i]), predicted_sin=float(-1 / (math.pi * k * p ** (k / 2))), fitted_cos=float(coef[2 * i + 1])) for i, (p, k) in enumerate(pp) if (p, k) in ((2, 1), (3, 1), (5, 1), (7, 1), (11, 1), (13, 1), (2, 2), (2, 3), (3, 2))])
json.dump(R, open(root / 'experiments/block_flow_and_prime_flow.json', 'w'), indent=1, default=float)
print(json.dumps(R['B1_summary'], indent=0)); print([ (b['k'], round(b['span'],1), round(b['turns'],2), b['cells'], b['new_cells'], round(b['re_neg'],3), round(b['eig_big'],2), round(b['axis_deg'],0)) for b in B[:8]])
print('L1 dev', R['L1_mean_abs_dev_from_1_at_100'], R['L1_mean_abs_dev_from_1_at_25000']); print([(x['n'], [round(v,3) for v in x['ratio_at_100_1000_5000_25000']]) for x in rows[:9]])
print('S1', R['S1']['var_S'], R['S1']['var_pred'], R['S1']['R2']); [print(f) for f in R['S1']['fit']]
